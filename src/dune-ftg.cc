// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>

// use slightly modified versions of dune/modelling/boundary.hh and ./equation.hh
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
#include<fstream>
#include<algorithm>

using namespace Dune::Modelling;

void transientTransport(int argc, char** argv)
{
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
  using TransportModel    = ForwardModel<ModelTraits,ModelTypes::Transport,   Formulation::Transient>;
  using GeoelectricsModel = ForwardModel<ModelTraits,ModelTypes::Geoelectrics,Formulation::Stationary>;

  // insert groundwater and transport model into forward model list
  forwardModelList.add<GroundwaterModel>("groundwater");
  forwardModelList.add<TransportModel,GroundwaterModel>("transport",std::list<std::string>{"groundwater"});
  
  set_electrodes<ModelTraits>(&modelTraits);
  set_wells<ModelTraits>(&modelTraits);

  //generate N geoelectrics models for N electrodes
  std::stringstream temp_ss;
  for (int i = 0; i < modelTraits.electrodeconfiguration.no_electrodes; i++)
  {
    temp_ss.str(std::string());
    temp_ss.clear();
    temp_ss << "geoelectrics" << i; //name for the n-th model is geoelectricsn
    forwardModelList.add<GeoelectricsModel,TransportModel>(temp_ss.str(),std::list<std::string>{"transport"}); //create a new geoelectrics model
  }
  
  // print information about model list, avoid multiple outputs if run is parallel
  if (helper.rank() == 0)
    forwardModelList.report(std::cout);

  // define parameters (input) and measurements (output)
  using ParameterList   = ModelTraits::ParameterList;
  using MeasurementList = ModelTraits::MeasurementList;
  std::shared_ptr<ParameterList>   parameterList  (new ParameterList(config.template get<std::string>("fields.location")));
  //std::shared_ptr<ParameterList>   parameterList  (new ParameterList(config));
  std::shared_ptr<MeasurementList> measurementList(new MeasurementList());
  //(*parameterList).generate(); 

  // perform forward run
  forwardModelList.solve(parameterList,measurementList);
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
