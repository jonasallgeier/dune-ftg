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

void transientTransport(int argc, char** argv)
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
  std::shared_ptr<MeasurementList> measurementList(new MeasurementList());

  // perform forward run
  forwardModelList.solve(parameterList,measurementList);
  
  // give output of total elapsed time  
  totalTimer.stop();
  if (helper.rank() == 0)
    std::cout << "Total time elapsed: " << totalTimer.elapsed() << std::endl;
}


void printParameters(int argc, char** argv)
{
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

  // read in configuration from .ini file
  Dune::ParameterTree config;
  Dune::ParameterTreeParser parser;
  parser.readINITree("modelling.ini",config);
   
  // use double for coordinates and values, dim = 3;
  using Traits = ModelTraits<double,double,3>;
  using GridTraits = typename Traits::GridTraits;
  using GV    = typename GridTraits::Grid::LeafGridView;
  using DomainField = typename GridTraits::DomainField;
  using RangeField  = typename GridTraits::RangeField;
  using ParameterList   = Traits::ParameterList;
  using ParameterField  = typename ParameterList::SubRandomField;

  Traits   modelTraits(helper,config);  //this will also read in the electrode configuration file
  std::shared_ptr<ParameterList>  parameterList  (new ParameterList(config.template get<std::string>("fields.location")));
  std::shared_ptr<ParameterField> conductivityField;
  conductivityField = (*parameterList).get(argv[2]);

  const GV gv = modelTraits.grid().leafGridView();

  auto glambda = [&](const Traits::GridTraits::Domain& x)
  {
    Traits::GridTraits::Scalar s;
    (*conductivityField).evaluate(x,s);
    return s;
  };
  auto g = Dune::PDELab::makeGridFunctionFromCallable(gv,glambda);

  enum {dim = GridTraits::dim};

  using FEM = Dune::PDELab::P0LocalFiniteElementMap<DomainField,RangeField,dim>;
  using VBE = Dune::PDELab::istl::VectorBackend<Dune::PDELab::istl::Blocking::fixed,1>;
  using CON = Dune::PDELab::P0ParallelConstraints;
  using GFS = Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE>;
  using Z = Dune::PDELab::Backend::Vector<GFS,RangeField>; // A coefficient vector
  using ZDGF = Dune::PDELab::DiscreteGridFunction<GFS,Z>; // Make a grid function out of it
  using VTKF =  Dune::PDELab::VTKGridFunctionAdapter<ZDGF>;

  FEM fem(Dune::GeometryType(Dune::GeometryType::cube,dim));
  GFS gfs(gv,fem);
  gfs.name(argv[2]);
  Z z(gfs); // initial value
  ZDGF zdgf(gfs,z);
  Dune::PDELab::interpolate(g,gfs,z); // Fill the coefficient vector

  Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,1);
  vtkwriter.addVertexData(std::shared_ptr<VTKF>(new VTKF(zdgf,"data")));
  vtkwriter.pwrite(argv[2],"vtk","", Dune::VTK::appendedraw);
}


int main(int argc, char** argv)
{
  try
  {
    if (argc==1)
      transientTransport(argc,argv); // try to run the problem
    else if (std::string(argv[1]) == "print")
      printParameters(argc,argv);
    else
      std::cout << "Possible options: \n"
      << "      (no option) -> run dune-ftg model\n"
      << "      print [field] -> print parameter field as VTK\n"
      << std::endl;
    return 0;
  }
  catch (Dune::Exception &e)
  {
    std::cerr << "Dune reported error: " << e << std::endl;
  }
}
