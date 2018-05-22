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

#include<dune/ftg/modeltraits_moments.hh>
#include<dune/ftg/groundwater_moments.hh>
#include<dune/ftg/moments_c.hh>
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
  ModelTraits   modelTraits(helper,config);  //this will also read in the electrode configuration file

  // define forward model
  using ForwardModelList = ForwardModelList<ModelTraits>;
  ForwardModelList   forwardModelList(modelTraits);
  
  // define forward models for groundwater and transport
  using GroundwaterModel  = ForwardModel<ModelTraits,ModelTypes::Groundwater, Formulation::Stationary>;
  using Moments_c_Model    = ForwardModel<ModelTraits,ModelTypes::Moments_c,Formulation::Stationary>;

  // insert groundwater and transport models into forward model list
  forwardModelList.add<GroundwaterModel>("groundwater");
  //forwardModelList.add<Moments_c_Model,GroundwaterModel>("moments_c",std::list<std::string>{"groundwater"});
  
  unsigned int highest_moment = config.template get<unsigned int>("moments.highest");
  for (unsigned int i = 0; i <= highest_moment; i++)
  {
    std::string model_name = "moments_c";
    model_name += std::to_string(i); //name for the k-th model is moments_ck
    if (i == 0)
    {
      forwardModelList.add<Moments_c_Model,GroundwaterModel>(model_name,std::list<std::string>{"groundwater"}); //create a new transport moment model
    }
    else
    {
      std::string model_name_dependency = "moments_c";
      model_name_dependency += std::to_string(i-1); // a moment model depends on the model below   
      forwardModelList.add<Moments_c_Model,GroundwaterModel,Moments_c_Model>(model_name,std::list<std::string>{"groundwater",model_name_dependency}); //create a new transport moment model
    }  
  }

  
  set_electrodes<ModelTraits>(&modelTraits);
  set_wells<ModelTraits>(&modelTraits);
  
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
  
  // give output of total elapsed time  
  totalTimer.stop();
  if (helper.rank() == 0)
    std::cout << "Total time elapsed: " << totalTimer.elapsed() << std::endl;
}

int main(int argc, char** argv)
{
  try
  {
    moments(argc,argv); // try to run the problem
    return 0;
  }
  catch (Dune::Exception &e)
  {
    std::cerr << "Dune reported error: " << e << std::endl;
  }
}
