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
void printParameters(int argc, char** argv,bool logprint)
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

  Traits   modelTraits(helper,config,false);  //this will also read in the electrode configuration file
  std::shared_ptr<ParameterList>  parameterList  (new ParameterList(config.template get<std::string>("fields.location")));
  std::shared_ptr<ParameterField> parameterField;
  if (logprint)
    parameterField = (*parameterList).get(argv[2]);
  else
    parameterField = (*parameterList).get(argv[1]);

  const GV gv = modelTraits.grid().leafGridView();

  auto glambda = [&](const Traits::GridTraits::Domain& x)
  {
    Traits::GridTraits::Scalar s;
    (*parameterField).evaluate(x,s);
    if (logprint)
      s = log10(s);
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
  if (logprint)
    gfs.name(argv[2]);
  else
    gfs.name(argv[1]);

  Z z(gfs); // initial value
  ZDGF zdgf(gfs,z);
  Dune::PDELab::interpolate(g,gfs,z); // Fill the coefficient vector

  Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,0);
  vtkwriter.addVertexData(std::shared_ptr<VTKF>(new VTKF(zdgf,"data")));
  if (logprint)  
    vtkwriter.pwrite(argv[2],"vtk","", Dune::VTK::appendedraw);
  else
    vtkwriter.pwrite(argv[1],"vtk","", Dune::VTK::appendedraw);
}


int main(int argc, char** argv)
{
  try
  {
    if (argc==1)
     std::cout << "Possible options: \n"
      << "      [field] -> print parameter field as VTK\n"
      << "      --logprint [field] -> print log of  parameter field as VTK\n"
      << std::endl;
    else if (std::string(argv[1]) == "--logprint")
      printParameters(argc,argv,true);
    else
      printParameters(argc,argv,false);
    return 0;
  }
  catch (Dune::Exception &e)
  {
    std::cerr << "Dune reported error: " << e << std::endl;
  }
}
