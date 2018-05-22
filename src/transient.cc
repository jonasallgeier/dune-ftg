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

#include<dune/ftg/modeltraits_transient.hh>
#include<dune/ftg/groundwater.hh>
#include<dune/ftg/transport.hh>
#include<dune/ftg/geoelectrics.hh>
#include<dune/ftg/ftg.hh>
#include<dune/pdelab/function/callableadapter.hh>

using namespace Dune::Modelling;

void transient(int argc, char** argv)
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
  
  bool evaluateBasePotentials;
  if (std::string(argv[1]) == "--basepotential")
    evaluateBasePotentials = true;
  else
    evaluateBasePotentials = false;

  ModelTraits   modelTraits(helper,config,evaluateBasePotentials);  //this will also read in the electrode configuration file
    
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
  std::shared_ptr<MeasurementList> measurementList(new MeasurementList(modelTraits));

  // perform forward run
  forwardModelList.solve(parameterList,measurementList);
  
  // give output of total elapsed time  
  totalTimer.stop();
  if (helper.rank() == 0)
    std::cout << "Total time elapsed: " << totalTimer.elapsed() << std::endl;
  
  if (helper.rank()== 0 && modelTraits.config().template get<bool>("output.unify_parallel_results",false))
  {
    std::string filenamebase = modelTraits.config().template get<std::string>("output.writeGeoelectricsFilename","results");
    
    std::string timefilename;
    timefilename.append(filenamebase);
    timefilename.append(".times");          
    
    std::string temp;

    std::ifstream timefile(timefilename);
    unsigned int no_electrodes;
    unsigned int no_processors;
    timefile >> temp;
    timefile >> no_electrodes;
    timefile >> temp;
    timefile >> no_processors;
    timefile >> temp;

    std::string time;
    std::vector<std::string> times;
    while ( !timefile.eof() )
    {
      timefile >> time;
      times.push_back(time);
    }
    using RF  = ModelTraits::GridTraits::RangeField;

    for (auto const& timeString: times) 
    {
      std::map<std::pair<unsigned int, unsigned int>, RF > complete_map; // <injection electrode, <measured electrode, potential> > 
    
      for( unsigned int proc = 0; proc < no_processors; ++proc ) 
      {
        std::string infilename;
        infilename.append(filenamebase);

        if (modelTraits.basePotentialEvaluation)
        {
          infilename.append("_base");
        } else
        {
          infilename.append("_");
          infilename.append(timeString);
        }

        infilename.append("_");
        std::stringstream ss;
        ss << proc;
        std::string rank(ss.str());
        infilename.append(rank);
        infilename.append(".data");


        std::ifstream infile(infilename);

        //skip header line
        infile >> temp;
        infile >> temp;
        infile >> temp;

        unsigned int injection_electrode;
        unsigned int measured_electrode;
        RF potential;

        while ( !infile.eof() )
        {
          infile >> injection_electrode;
          infile >> measured_electrode;
          infile >> potential;

          std::pair<unsigned int, unsigned int> tmp_pair;
          tmp_pair = std::make_pair(injection_electrode,measured_electrode);

          complete_map.insert(std::pair<std::pair<unsigned int, unsigned int>, RF >(tmp_pair, potential));
        }
        infile.close();
        remove(infilename.c_str());
        
        std::ofstream outfile;
        std::string outfilename;
        outfilename.append(filenamebase);
        if (modelTraits.basePotentialEvaluation)
        {
          outfilename.append("_base");
        } else
        { 
          outfilename.append("_");
          outfilename.append(timeString);
        }
        outfilename.append(".data");

        outfile.open(outfilename, std::ios::out | std::ios::trunc);
        outfile << "injected_el measured_el potential" << std::endl;

        for (auto const & current_entry : complete_map)
        {
          outfile << current_entry.first.first << " " << current_entry.first.second << " " << current_entry.second << std::endl;
        }
        outfile.close();
        remove(timefilename.c_str());
      }
    }
  }

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

  Traits   modelTraits(helper,config,false);  //this will also read in the electrode configuration file
  std::shared_ptr<ParameterList>  parameterList  (new ParameterList(config.template get<std::string>("fields.location")));
  std::shared_ptr<ParameterField> conductivityField;
  conductivityField = (*parameterList).get(argv[2]);

  const GV gv = modelTraits.grid().leafGridView();

  auto glambda = [&](const Traits::GridTraits::Domain& x)
  {
    Traits::GridTraits::Scalar s;
    (*conductivityField).evaluate(x,s);
    if (std::string(argv[1]) == "--logprint")
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
  gfs.name(argv[2]);
  Z z(gfs); // initial value
  ZDGF zdgf(gfs,z);
  Dune::PDELab::interpolate(g,gfs,z); // Fill the coefficient vector

  Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,0);
  vtkwriter.addVertexData(std::shared_ptr<VTKF>(new VTKF(zdgf,"data")));
  vtkwriter.pwrite(argv[2],"vtk","", Dune::VTK::appendedraw);
}


int main(int argc, char** argv)
{
  try
  {
    if (argc==1 || std::string(argv[1]) == "--basepotential")
      transient(argc,argv); // try to run the problem
    else if (std::string(argv[1]) == "--print")
      printParameters(argc,argv);
    else if (std::string(argv[1]) == "--logprint")
      printParameters(argc,argv);
    else
      std::cout << "Possible options: \n"
      << "      (no option) -> run dune-ftg model\n"
      << "      --basepotential -> run dune-ftg model, but only the electrical base potential will be evaluated; therefore adjust the modelling.ini to equal time steps\n"
      << "      --print [field] -> print parameter field as VTK\n"
      << "      --logprint [field] -> print log of  parameter field as VTK\n"
      << std::endl;
    return 0;
  }
  catch (Dune::Exception &e)
  {
    std::cerr << "Dune reported error: " << e << std::endl;
  }
}
