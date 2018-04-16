// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>

#include<dune/ftg/boundary.hh>
#include<dune/modelling/forwardmodel.hh>
#include<dune/modelling/forwardadjointmodel.hh>

#include<dune/ftg/modeltraits.hh>
#include<dune/ftg/groundwater.hh>
#include<dune/ftg/transport.hh>
#include<dune/ftg/geoelectrics.hh>
#include<dune/ftg/ftg.hh>
#include<fstream> //jonas

using namespace Dune::Modelling;

void transientTransport(int argc, char** argv)
{
  // initialize MPI if available
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

  // read in configuration from .ini file
  Dune::ParameterTree config;
  Dune::ParameterTreeParser parser;
  parser.readINITree("modelling.ini",config);
   
  // use double for coordinates and values, dim = 3; instantiate modelTraits -> this will already read in the electrode configuration file, as the model traits contain all properties, which are common to all models
  using ModelTraits = ModelTraits<double,double,3>;
  ModelTraits   modelTraits(helper,config);  

  // define forward model
  using ForwardModelList = ForwardModelList<ModelTraits>;
  ForwardModelList   forwardModelList(modelTraits);
  
  // define forward models for groundwater and transport
  using GroundwaterModel  = ForwardModel<ModelTraits,ModelTypes::Groundwater, Formulation::Transient>;
  using TransportModel    = ForwardModel<ModelTraits,ModelTypes::Transport,   Formulation::Transient>;
  using GeoelectricsModel = ForwardModel<ModelTraits,ModelTypes::Geoelectrics,Formulation::Stationary>; //jonas

  // insert groundwater and transport model into forward model list
  forwardModelList.add<GroundwaterModel>("groundwater");
  forwardModelList.add<TransportModel,GroundwaterModel>("transport",std::list<std::string>{"groundwater"});
  
  setelectrodes<ModelTraits>(&modelTraits);
  setwells<ModelTraits>(&modelTraits);

  //Jonas; add many geoelectrics models
  std::stringstream temp_ss;
  std::string temp_name;
  int no_electrodes = modelTraits.electrodeconfiguration.no_electrodes; // total number of electrodes
  
  //generate N models for N electrodes
  for (int i = 0; i < no_electrodes; i++) //jonas
  {
    //name for the n-th model is geoelectricsn
    temp_ss.str(std::string());
    temp_ss.clear();
    temp_ss << "geoelectrics";
    temp_ss << i;
    temp_name = temp_ss.str();

    //create a new geoelectrics model
    forwardModelList.add<GeoelectricsModel,TransportModel>(temp_name,std::list<std::string>{"transport"});
  }
  
  // print information about model list, avoid multiple outputs if run is parallel
  if (helper.rank() == 0)
  forwardModelList.report(std::cout);

  // define parameters (input) and measurements (output)
  using ParameterList   = ModelTraits::ParameterList;
  using MeasurementList = ModelTraits::MeasurementList;
  std::shared_ptr<ParameterList>   parameterList  (new ParameterList("./fields/ftg"));
  //std::shared_ptr<ParameterList>   parameterList  (new ParameterList(config));
  std::shared_ptr<MeasurementList> measurementList(new MeasurementList());
  //(*parameterList).generate(); 

  // perform forward run
  forwardModelList.solve(parameterList,measurementList);
  // perform forward run, one model after the other
  //forwardModelList.solveConsecutively(parameterList,measurementList);
  // perform forward run, step all models at same time (default)
  //forwardModelList.solveSimultaneously(parameterList,measurementList);
}

int main(int argc, char** argv)
{
  try
  {
    transientTransport(argc,argv); // try to run the problem
    return 0;
  }
  catch (Dune::Exception &e)
  {
    std::cerr << "Dune reported error: " << e << std::endl;
  }
}
