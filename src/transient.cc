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
#include<dune/ftg/transport.hh>
#include<dune/ftg/geoelectrics.hh>
#include<dune/ftg/ftg.hh>
#include<dune/pdelab/function/callableadapter.hh>

using namespace Dune::Modelling;

void transient(int argc, char** argv, bool evaluateBasePotentials)
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

  ModelTraits   modelTraits(helper,config,evaluateBasePotentials);  //this will also read in the electrode configuration file
  
  // define forward model
  using ForwardModelList = ForwardModelList<ModelTraits>;
  ForwardModelList   forwardModelList(modelTraits);
  
  // define forward models for groundwater and transport
  using GroundwaterModel  = ForwardModel<ModelTraits,ModelTypes::Groundwater, Formulation::Stationary>;
  using TransportModel    = ForwardModel<ModelTraits,ModelTypes::Transport,   Formulation::Transient>;
  using GeoelectricsModel = ForwardModel<ModelTraits,ModelTypes::Geoelectrics,Formulation::Stationary>;

  // insert groundwater and transport model into forward model list
  if (!modelTraits.basePotentialEvaluation)
  {  
    forwardModelList.add<GroundwaterModel>("groundwaterFlow");
    forwardModelList.add<TransportModel,GroundwaterModel>("soluteTransport",std::list<std::string>{"groundwaterFlow"});
  }

  set_electrodes<ModelTraits>(&modelTraits);
  set_wells<ModelTraits>(&modelTraits);
  
  //generate N geoelectrics models for N electrodes
  std::stringstream temp_ss;
  for (int i = 0; i < modelTraits.electrodeconfiguration.no_electrodes; i++)
  {
    std::string name = "ERT_" + std::to_string(i); //name for the n-th model is geoelectricsn
    if (!modelTraits.basePotentialEvaluation)
      forwardModelList.add<GeoelectricsModel,TransportModel>(name,std::list<std::string>{"soluteTransport"}); //create a new geoelectrics model
    else
      forwardModelList.add<GeoelectricsModel>(name);
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

    if (modelTraits.config().template get<bool>("output.writeTransport",false) && !modelTraits.basePotentialEvaluation)
      unify_transport_results<ModelTraits>(&modelTraits);
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
      transient(argc,argv,false); // try to run the problem
    else if (std::string(argv[1]) == "--basepotential")
      transient(argc,argv,true);
    else
      std::cout << "Possible options: \n"
      << "      (no option) -> run transient model\n"
      << "      --basepotential -> run transient model, but only the electrical base potential will be evaluated\n"
      << std::endl;
    return 0;
  }
  catch (Dune::Exception &e)
  {
    std::cerr << "Dune reported error: " << e << std::endl;
  }
}
