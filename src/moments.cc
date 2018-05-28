// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>

// use slightly modified versions of dune/modelling/boundary.hh, ./equation.hh and ./forwardmodel.hh
#include<dune/ftg/override/boundary.hh>
#include<dune/ftg/override/equation.hh>
#include<dune/ftg/override/forwardmodel.hh>

#include<dune/modelling/forwardmodel.hh>
#include<dune/modelling/forwardadjointmodel.hh>

#include<dune/ftg/modeltraits.hh>
#include<dune/ftg/groundwater.hh>
#include<dune/ftg/moments_c.hh>
#include<dune/ftg/transport.hh>
#include<dune/ftg/geoelectrics.hh>
#include<dune/ftg/moments_ERT.hh>
#include<dune/ftg/ftg.hh>
#include<dune/pdelab/function/callableadapter.hh>

using namespace Dune::Modelling;

void moments(int argc, char** argv)
{
  Dune::Timer totalTimer;
  totalTimer.start();

  // initialize MPI if available
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

  // read in configuration from .ini file
  Dune::ParameterTree config;
  Dune::ParameterTreeParser parser;
  parser.readINITree("modelling.ini",config);

  // use double for coordinates and values, dim = 3;
  using ModelTraits = ModelTraits<double,double,3>;
  ModelTraits   modelTraits(helper,config,true);  //this will also read in the electrode configuration file

  // define forward model
  using ForwardModelList = ForwardModelList<ModelTraits>;
  ForwardModelList   forwardModelList(modelTraits);
  
  // define forward models
  using GroundwaterModel  = ForwardModel<ModelTraits,ModelTypes::Groundwater, Formulation::Stationary>;
  using Moments_c_Model    = ForwardModel<ModelTraits,ModelTypes::Moments_c,Formulation::Stationary>;
  using GeoelectricsModel = ForwardModel<ModelTraits,ModelTypes::Geoelectrics,Formulation::Stationary>;
  using Moments_ERT_Model = ForwardModel<ModelTraits,ModelTypes::Moments_ERT,Formulation::Stationary>;

  forwardModelList.add<GroundwaterModel>("groundwaterFlow");
  unsigned int highest_moment = config.template get<unsigned int>("moments.highest");
  for (unsigned int i = 0; i <= highest_moment; i++)
  {
    std::string model_name = "momentsTransport_" + std::to_string(i); //name for the k-th model is moments_ck
    if (i == 0)
    {
      forwardModelList.add<Moments_c_Model,GroundwaterModel>(model_name,std::list<std::string>{"groundwaterFlow"}); //create a new transport moment model
    }
    else
    {
      std::string model_name_dependency = "momentsTransport_" +  std::to_string(i-1); // a moment model depends on the model below
      forwardModelList.add<Moments_c_Model,GroundwaterModel,Moments_c_Model>(model_name,std::list<std::string>{"groundwaterFlow",model_name_dependency});
    }  
  }
  set_electrodes<ModelTraits>(&modelTraits);
  set_wells<ModelTraits>(&modelTraits);

  //generate N geoelectrics models for N electrodes
  for (int i = 0; i < modelTraits.electrodeconfiguration.no_electrodes; i++)
  {
    std::string name = "ERT_" + std::to_string(i); //name for the n-th model is geoelectricsn
    forwardModelList.add<GeoelectricsModel>(name); //create a new geoelectrics model
  }

  //generate (k+1)*N ERT moment models for N electrodes and the highest moment k
  for (int i = 0; i < modelTraits.electrodeconfiguration.no_electrodes; i++)
  {
    for (unsigned int j = 0; j <= highest_moment; j++)
    {
      std::string model_name = "momentsERT_" + std::to_string(i) + "_" + std::to_string(j); // example: momentsERT_42_1 -> 1st moment, 42nd electrode
      std::string c_model_name = "momentsTransport_"    + std::to_string(j);
      std::string ERT_model_name = "ERT_" + std::to_string(i);
      forwardModelList.add<Moments_ERT_Model,Moments_c_Model,GeoelectricsModel>(model_name,std::list<std::string>{c_model_name,ERT_model_name});
    }
  }

  // print information about model list, avoid multiple outputs if run is parallel
  if (helper.rank() == 0)
    forwardModelList.report(std::cout);

  // define parameters (input) and measurements (output)
  using ParameterList   = ModelTraits::ParameterList;
  using MeasurementList = ModelTraits::MeasurementList;
  std::shared_ptr<ParameterList>   parameterList  (new ParameterList(config.template get<std::string>("fields.location")));
  std::shared_ptr<MeasurementList> measurementList(new MeasurementList(modelTraits));

  // perform forward run
  forwardModelList.solve(parameterList,measurementList);

  if (helper.rank()== 0 && modelTraits.config().template get<bool>("output.unify_parallel_results",false))
  {
    if (modelTraits.config().template get<bool>("output.writeERT",false))
      unify_ERT_results<ModelTraits>(&modelTraits);
    if (modelTraits.config().template get<bool>("output.writeMomentsTransport",false))
      unify_momentsTransport_results<ModelTraits>(&modelTraits);
    if (modelTraits.config().template get<bool>("output.writeMomentsERT",false))
      unify_momentsERT_results<ModelTraits>(&modelTraits);
  }

  // give output of total elapsed time  
  totalTimer.stop();
  if (helper.rank() == 0)
    std::cout << "Total time elapsed: " << totalTimer.elapsed() << std::endl;
}

int main(int argc, char** argv)
{
  try
  {
    if (argc==1)
      moments(argc,argv); // try to run the problem
    else
      std::cout << "Possible options: \n"
      << "      (no option) -> run moment model\n"
      << std::endl;
    return 0;
  }
  catch (Dune::Exception &e)
  {
    std::cerr << "Dune reported error: " << e << std::endl;
  }
}
