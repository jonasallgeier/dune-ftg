// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_FTG_MODELTRAITS_TRANSIENT_HH
#define DUNE_FTG_MODELTRAITS_TRANSIENT_HH

#include<dune/grid/yaspgrid.hh>
// use slightly modified version of dune/randomfield.hh
#include<dune/ftg/override/randomfield.hh>

// the following includes are necessary for ERT matrix storage
#include<dune/pdelab/common/referenceelements.hh>
#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/constraints/p0.hh>
#include<dune/pdelab/backend/istl.hh>

/**
 * @brief Helper class for grid generation
 */
template<typename GT>
class GridHelper
{
  enum {dim = GT::dim};

  typedef typename GT::DomainField DF;
  typedef typename GT::RangeField RF;

  int levels;
  std::vector<unsigned int> minCells, maxCells;
  std::array<std::vector<DF>,dim> coordinate_vectors;

  public:

  //std::vector<DF> maxExt;
  bool isequidistant;
  std::vector<DF> lowerleft, upperright;

  GridHelper(const Dune::ParameterTree& config)
  {
    isequidistant = config.get<bool>("grid.equidistant",1);
    levels      = config.get<int>("grid.levels",1);
    if (isequidistant)       
    {
      lowerleft   = config.get<std::vector<DF> >  ("grid.lowerleft");   // get coordinates of lower left point
      upperright  = config.get<std::vector<DF> >  ("grid.upperright");  // get coordinates of upper right point
      maxCells = config.get<std::vector<unsigned int> >("grid.cells");
      minCells = maxCells;

      for (int i = 0; i < levels - 1; i++)
        for (unsigned int j = 0; j < maxCells.size(); j++)
        {
          if (minCells[j]%2 != 0)
            DUNE_THROW(Dune::Exception,"cannot create enough levels for hierarchical grid, check number of cells");

          minCells[j] /= 2;
        }
    } else {
      for (unsigned int i = 0; i < dim; i++)
      {
        std::vector<DF> vec_i;
        if (i==0)
          vec_i = config.get<std::vector<DF> >("grid.vector_x");
        else if (i==1)
          vec_i = config.get<std::vector<DF> >("grid.vector_y");
        else if (i==2)
          vec_i = config.get<std::vector<DF> >("grid.vector_z");
        else
          DUNE_THROW(Dune::Exception,"nonequidistant grid with more than four dimensions is not implemented");
        lowerleft.push_back(vec_i[0]);
        upperright.push_back(vec_i.back());
        coordinate_vectors[i] = vec_i;
      }
    }
  }

  // lower left point of grid
  Dune::FieldVector<DF,dim> LL() const
  {
    Dune::FieldVector<DF,dim> LLvector;

    for (unsigned int i = 0; i < dim; i++)
      LLvector[i] = lowerleft[i];

    return LLvector;
  }

  // upper right point of grid
  Dune::FieldVector<DF,dim> UR() const
  {
    Dune::FieldVector<DF,dim> URvector;

    for (unsigned int i = 0; i < dim; i++)
      URvector[i] = upperright[i];

    return URvector;
  }

  // array, containing vectors of all nodal coordinates
  std::array<std::vector<DF>,dim> coords() const
  {
    std::array<std::vector<DF>,dim> output;
    if (isequidistant)
    {
      for (unsigned int i = 0; i < dim; i++)
      {
        DF delta = (upperright[i]-lowerleft[i])/minCells[i];
        std::vector<DF> vec_i;
        vec_i.reserve(minCells[i]);
        DF pos = lowerleft[i];
        for (unsigned int j = 0; j <= minCells[i]; j++)
        {
          vec_i.push_back(pos);
          pos += delta;
        }
        output[i] = vec_i;
      }
    } else
    {
      output = coordinate_vectors;
    }
    return output;
  }

  std::array<int,dim> N() const
  {
    std::array<int,dim> Nvector;

    for (unsigned int i = 0; i < dim; i++)
      Nvector[i] = minCells[i];

    return Nvector;
  }

  std::bitset<dim> B() const
  {
    return std::bitset<dim>(false);
  }
};

/**
 * @brief Traits class defining central constants and types
 */
template<typename DF, typename RF, unsigned int dimension>
class ModelTraits
{
  private:
  std::vector<unsigned int> electrode_cell_indices;
  std::map<unsigned int, std::pair<RF, bool>> well_cells;
  
  public:
    // traits relating to geometry and underlying types
    struct GridTraits
    {
      enum {dim = dimension};
      using Grid     = Dune::YaspGrid<dim,Dune::TensorProductCoordinates<DF,dim>>;
      using GridView = typename Grid::LevelGridView;
      using RangeField   = RF;
      using Scalar       = Dune::FieldVector<RF,1>;
      using Vector       = Dune::FieldVector<RF,dim>;
      using Tensor       = Dune::FieldMatrix<RF,dim,dim>;
      using DomainField  = DF;
      using Domain       = Dune::FieldVector<DF,dim>;
      using IDomain      = Dune::FieldVector<DF,dim-1>;
      using Element      = typename GridView::Traits::template Codim<0>::Entity;
      using Intersection = typename GridView::Intersection;
    };

    // use a list of random Gaussian fields as parameter list
    using ParameterList = Dune::RandomField::RandomFieldList<GridTraits>;

    const Dune::MPIHelper& helper;
    const Dune::ParameterTree& duneConfig;
    const GridHelper<GridTraits> gh;

    bool basePotentialEvaluation;

    typename GridTraits::Grid yaspGrid;

    ModelTraits(Dune::MPIHelper& helper_, Dune::ParameterTree& duneConfig_, bool basePotentialEvaluation_)
      : helper(helper_), duneConfig(duneConfig_), gh(duneConfig), basePotentialEvaluation(basePotentialEvaluation_),
      yaspGrid(gh.coords(),gh.B(),1)
    {
    }

    const std::vector<DF> & corners_ll() const
    {
      return gh.lowerleft;
    }

    const std::vector<DF> & corners_ur() const
    {
      return gh.upperright;
    }

    const typename GridTraits::Grid& grid() const
    {
      return yaspGrid;
    }

    const Dune::ParameterTree& config() const
    {
      return duneConfig;
    }

    int rank() const
    {
      return helper.rank();
    }

    const auto& comm() const
    {
      return yaspGrid.comm();
    }

    // provide a leaf grid view
    const typename GridTraits::Grid::LeafGridView& gv() const
    {
      return yaspGrid.leafGridView();
    }

    // provide the vector of grid indices that contain electrodes
    std::vector<unsigned int> read_electrode_cell_indices() const
    {
      return electrode_cell_indices;
    }

    void set_electrode_cell_indices(std::vector<unsigned int> cell_indices)
    {
      electrode_cell_indices = cell_indices;
    }

    void set_well_cells(std::map<unsigned int, std::pair<RF, bool>> in_well_cells)
    {
      well_cells = in_well_cells;
    }

    std::map<unsigned int, std::pair<RF, bool>> read_well_cells() const
    {
      return well_cells;
    }

    // define electrode configuration; it is read in from a file by the constructor of this struct
    struct ElectrodeConfiguration
    {
        int no_electrodes;          // number of electrodes
        std::vector<int> name_elec; // name vector
        std::vector<RF> x_elec;     // vector of x coordinates
        std::vector<RF> y_elec;     // vector of y coordinates
        std::vector<RF> z_elec;     // vector of z coordinates
        std::vector<int> surf_elec; // vector containing information if electrode is a surface electrode 

        // constructor -> reads in the data from file
        ElectrodeConfiguration (std::string electrodefilename, int myrank) 
        {
	      std::ifstream file_econf; // open the file
          file_econf.open(electrodefilename);
          if (myrank == 0)
            {std::cout << "Attempting to read electrode configuration file...";}          
          if (file_econf.is_open())
          {           
            file_econf >> no_electrodes;  // first line equals number of electrodes
   
            // declare temporary storage variables
            int tmp_name;
            RF tmp_x;
            RF tmp_y;
            RF tmp_z;
            int tmp_surf;
        
            // read data in lines {#,x,y,z,s} from file via temp to vector 
            for (int i = 1; i < no_electrodes+1; i++) 
            {
              file_econf >> tmp_name;
              name_elec.push_back(tmp_name);
              file_econf >> tmp_x;
              x_elec.push_back(tmp_x);
              file_econf >> tmp_y;
              y_elec.push_back(tmp_y);
              file_econf >> tmp_z;
              z_elec.push_back(tmp_z);
              file_econf >> tmp_surf;
              surf_elec.push_back(tmp_surf);
            }
            
            // close the file and give console output
            file_econf.close();
            if (myrank == 0)
              {std::cout << " done. Found " << z_elec.size() << " electrode(s)." << std::endl;}
          }
          else
          {
            DUNE_THROW(Dune::Exception,"Unable to read electrode configuration file. Please check file name.");
          }
        };
    };
    
    // call constructor to read in electrode configuration from file specified in .ini; data is stored in a struct object
    ElectrodeConfiguration electrodeconfiguration {duneConfig.get<std::string>("configfiles.electrodes"),rank()};

    // define well configuration; it is read in from a file by the constructor of this struct
    struct WellConfiguration
    {
        int no_wells;                     // number of wells
        std::vector<RF> x_well;           // vector of x coordinates
        std::vector<RF> y_well;           // vector of y coordinates
        std::vector<RF> z1_well;          // vector of z1 coordinates
        std::vector<RF> z2_well;          // vector of z2 coordinates
        std::vector<RF> q_well;           // vector of pumping rates
        std::vector<bool> injection_well; // vector of bools indicating whether this well participates the tracer injection 

        // constructor -> reads in the data from file
        WellConfiguration (std::string wellfilename, int myrank) 
        {
	      std::ifstream file_wconf; // open the file
          file_wconf.open(wellfilename);
          if (myrank == 0)
            {std::cout << "Attempting to read well configuration file...";}          
          if (file_wconf.is_open())
          {           
            file_wconf >> no_wells;  // first line equals number of wells
   
            // declare temporary storage variables
            RF tmp_x;
            RF tmp_y;
            RF tmp_z1;
            RF tmp_z2;
            RF tmp_q;
            bool tmp_in;
        
            // read data in lines {#,x,y,z,s} from file via temp to vector 
            for (int i = 1; i < no_wells+1; i++) 
            {
              file_wconf >> tmp_x;
              x_well.push_back(tmp_x);
              file_wconf >> tmp_y;
              y_well.push_back(tmp_y);
              file_wconf >> tmp_z1;
              z1_well.push_back(tmp_z1);
              file_wconf >> tmp_z2;
              z2_well.push_back(tmp_z2);
              file_wconf >> tmp_q;
              q_well.push_back(tmp_q);
              file_wconf >> tmp_in;
              injection_well.push_back(tmp_in);
            }
            
            // close the file and give console output
            file_wconf.close();
            if (myrank == 0)
              {std::cout << " done. Found " << q_well.size() << " well(s)." << std::endl;}
          }
          else
          {
            DUNE_THROW(Dune::Exception,"Unable to read well configuration file. Please check file name.");
          }
        };
    };
    
    // call constructor to read in well configuration from file specified in .ini; data is stored in a struct object
    WellConfiguration wellconfiguration {duneConfig.get<std::string>("configfiles.wells"),rank()};

    // this class stores the ERT matrix during a single ERT measurement
    class ERTMatrixContainer
    {
      private:
        const ModelTraits& traits;
        // this is bad! flexible code needed here! TODO
        using M = Dune::PDELab::ISTL::BCRSMatrix<Dune::PDELab::GridFunctionSpace<Dune::GridView<Dune::DefaultLevelGridViewTraits<const Dune::YaspGrid<3, Dune::TensorProductCoordinates<double, 3> > > >, Dune::PDELab::P0LocalFiniteElementMap<double, double, 3>, Dune::PDELab::P0ParallelConstraints, Dune::PDELab::ISTL::VectorBackend<(Dune::PDELab::ISTL::Blocking)2, 1ul>, Dune::PDELab::LeafOrderingTag<Dune::PDELab::EmptyParams> >, Dune::PDELab::GridFunctionSpace<Dune::GridView<Dune::DefaultLevelGridViewTraits<const Dune::YaspGrid<3, Dune::TensorProductCoordinates<double, 3> > > >, Dune::PDELab::P0LocalFiniteElementMap<double, double, 3>, Dune::PDELab::P0ParallelConstraints, Dune::PDELab::ISTL::VectorBackend<(Dune::PDELab::ISTL::Blocking)2, 1ul>, Dune::PDELab::LeafOrderingTag<Dune::PDELab::EmptyParams> >, Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>, std::allocator<Dune::FieldMatrix<double, 1, 1> > >, Dune::PDELab::ISTL::PatternStatistics<long unsigned int> >;
        M m;
      public:
      ERTMatrixContainer(const ModelTraits& traits_)
        : traits(traits_)
      {
      }
        void set_matrix(M m_in)
        {
          m = m_in;
        }
        M read_matrix()
        {
          return m;
        }
    };


    class MeasurementList
    {
      private:
        const ModelTraits& traits;

        // initial writing will clear the files, all other writings will append
        bool clearFiles_ERT               = true;
        bool clearFiles_Transport         = true;
        bool clearFiles_MomentsTransport  = true;
        bool clearFiles_MomentsERT        = true;

        std::map<unsigned int, std::map<unsigned int, RF> > output_ERT_all_electrodes;
   
        // < <moment k, inj el> <meas el, value> >
        std::map<std::pair<unsigned int, unsigned int>, std::map<unsigned int, RF> > output_MomentsERT_all_electrodes;  
        bool writeERT; 
        bool writeTransport; 
        bool writeMomentsTransport;
        bool writeMomentsERT;      

        auto& index_set() const
        {
          return traits.grid().leafGridView().indexSet();
        };

        template<typename Storage>
        unsigned int extraction_helper(std::map<unsigned int, RF> & output_this_electrode, const Storage& storage, const RF& time)
        {
          std::vector<unsigned int> electrode_cells = traits.read_electrode_cell_indices();
          const typename GridTraits::Domain x = {0.5, 0.5, 0.5}; // get value at cell center

          for (const auto & elem : elements (traits.grid().leafGridView()))
          {
            // check if this cell is a measurement cell/is multiple measurement cells --> get electrode numbers of current cell
            unsigned int cell_index = index_set().index(elem);
            std::vector<unsigned int> affected_electrodes; 
            
            auto lambda_compare = [cell_index](unsigned int n) { return (n==cell_index); };
            
            std::vector<unsigned int>::iterator it = std::find_if (electrode_cells.begin(), electrode_cells.end(), lambda_compare);
            while ((electrode_cells.end()) != it) 
            {
              affected_electrodes.push_back( distance(electrode_cells.begin(), it) );
              it++;
              it = std::find_if (it, electrode_cells.end(), lambda_compare);
            }

            // if there is one or more electrodes -> measure the potential & write output in output vector
            if (affected_electrodes.size() > 0)
            {
              typename GridTraits::Scalar output = 0.0;  
              (*storage).value(time,elem,x,output);
              for(std::vector<unsigned int>::iterator iter = affected_electrodes.begin(); iter != affected_electrodes.end(); ++iter)
              {
                output_this_electrode.insert ( std::pair<unsigned int,RF>(*iter+1,output) );
              }
            }
          }
          return electrode_cells.size();
        }

      public:
      MeasurementList(const ModelTraits& traits_)
        : traits(traits_)
      {
        writeERT              = traits.config().template get<bool>("output.writeERT",false); 
        writeTransport        = traits.config().template get<bool>("output.writeTransport",false); 
        writeMomentsTransport = traits.config().template get<bool>("output.writeMomentsTransport",false);
        writeMomentsERT       = traits.config().template get<bool>("output.writeMomentsERT",false);      
      }

        void setTimes(const RF& one, const RF& two)
        {}

        RF suggestPreviousTime(const RF& firstTime, const RF& lastTime) const
        {
          return firstTime;
        }

        RF suggestNextTime(const RF& firstTime, const RF& lastTime) const
        {
          return lastTime;
        }

        template<typename Storage>
        void extract(const Storage& storage, const RF& firstTime, const RF& time,
            const std::string & modelname, const unsigned int& modelNumber,const std::string& timeString, Dune::Timer&  printTimer)              
        {
          // determine which modeltype we have and if the output is desired; could probably be performed using templates (I don't know how)
          bool isERT              = (modelname.find("ERT_") == 0);
          bool isTransport        = (modelname.find("soluteTransport") == 0);
          bool isMomentsTransport = (modelname.find("momentsTransport_") == 0);
          bool isMomentsERT       = (modelname.find("momentsERT_") == 0);
          
          if (writeERT && isERT)
          {
            printTimer.start();
            
            // declare output map for this source electrode; [1:144]
            std::map<unsigned int, RF> current_output;
            unsigned int no_electrodes = extraction_helper<Storage>(current_output,storage,time); // get the results of this model

            // add the output for the active electrode to the map of all electrodes
            // increase modelNumber by one, as models are counted from 0, while electrode names start with 1
            output_ERT_all_electrodes.insert(std::pair<unsigned int,std::map<unsigned int, RF> >(modelNumber+1, current_output));

            // if this is the last ERT model print output to file
            if (modelNumber+1 == no_electrodes)
            {
              std::string filenamebase = traits.config().template get<std::string>("output.writeERTFilename","resultsERT");
              if (traits.basePotentialEvaluation)
              {
                filenamebase = traits.config().template get<std::string>("output.writeERTBaseFilename","resultsERTbase");
              }

              if (traits.rank() == 0)
              {
                std::cout << "printing ERT results" << std::endl;

                // print the temporal information to a timefile
                std::ofstream timefile;
                std::string timefilename = filenamebase + ".times";
                if (clearFiles_ERT)
                { 
                  timefile.open(timefilename, std::ios::out | std::ios::trunc); // clear the file
                  clearFiles_ERT = false;
                  timefile << "electrodes " << no_electrodes << std::endl;
                  timefile << "processors " << traits.helper.size() << std::endl;
                  timefile << "times";
                } else 
                {
                  timefile.open(timefilename, std::ios::app); // append to file
                }
                timefile << " " << time;
                timefile.close();
              }

              std::ofstream outfile;
              std::string filename = filenamebase;
              // if these are base potentials, don't name after time
              if (traits.basePotentialEvaluation)
              {
                filename += ("_base_" + std::to_string(traits.rank()) + ".data");
              } else
              {
                filename += ("_" + timeString + "_" + std::to_string(traits.rank()) + ".data");
              }

              outfile.open(filename, std::ios::out | std::ios::trunc);
              outfile << "injected_el measured_el potential" << std::endl;

              // go through the collected data set -> get measured potentials for one injection
              for (auto const & injection_electrode : output_ERT_all_electrodes)
              {
                for (auto const & measurement_electrode : injection_electrode.second)
                {
                  outfile<<injection_electrode.first<<" "<<measurement_electrode.first<<" "<<measurement_electrode.second<<std::endl;
                }
              }
              outfile.close();
              output_ERT_all_electrodes.clear(); // reset local storage of all ERT measurements
            }
            printTimer.stop();
          }

          // only save transport results at ERT measurement times          
          RF current_time = time;
          bool timeforTransport = (std::fmod(current_time, traits.config().template get<RF>("time.step_ERT"))  == 0);

          if (writeTransport && isTransport && timeforTransport)
          {
            printTimer.start();

            // declare output map for this source electrode; [1:144]
            std::map<unsigned int, RF> current_output;
            unsigned int no_electrodes = extraction_helper<Storage>(current_output,storage,time); // get the results of this model
            
            std::string filenamebase = traits.config().template get<std::string>("output.writeTransportFilename","resultsTransport");
            if (traits.rank() == 0)
            {
              std::cout << "printing concentration results" << std::endl;

              // print the temporal information to a timefile
              std::ofstream timefile;
              std::string timefilename = filenamebase + ".times";

              if (clearFiles_MomentsTransport)
              { 
                timefile.open(timefilename, std::ios::out | std::ios::trunc); // clear the file
                clearFiles_MomentsTransport = false;
                timefile << "electrodes " << no_electrodes << std::endl;
                timefile << "processors " << traits.helper.size() << std::endl;
                timefile << "times ";
              } else 
              {
                timefile.open(timefilename, std::ios::app); // append to file
              }
              timefile << " " << time;
              timefile.close();
            }

            std::ofstream outfile;
            std::string filename = filenamebase + "_" + timeString + "_" + std::to_string(traits.rank()) + ".data";

            outfile.open(filename, std::ios::out | std::ios::trunc);
            outfile << "electrode concentration" << std::endl;
            for (auto const & electrode : current_output)
            {
              outfile << electrode.first << " " << electrode.second << std::endl;
            }
            outfile.close();
          
            printTimer.stop();
          }

          if (writeMomentsTransport && isMomentsTransport)
          {
            printTimer.start();

            // declare output map for this source electrode; [1:144]
            std::map<unsigned int, RF> current_output;
            unsigned int no_electrodes = extraction_helper<Storage>(current_output,storage,time); // get the results of this model
            
            std::string filenamebase = traits.config().template get<std::string>("output.writeMomentsTransportFilename","resultsMomentsTransport");
            if (traits.rank() == 0)
            {
              std::cout << "printing transport moment(" << modelNumber << ") results" << std::endl;

              // print the temporal information to a timefile
              std::ofstream timefile;
              std::string timefilename = filenamebase + ".moments";

              if (clearFiles_Transport)
              { 
                timefile.open(timefilename, std::ios::out | std::ios::trunc); // clear the file
                clearFiles_Transport = false;
                timefile << "electrodes " << no_electrodes << std::endl;
                timefile << "processors " << traits.helper.size() << std::endl;
                timefile << "moments";
              } else 
              {
                timefile.open(timefilename, std::ios::app); // append to file
              }
              timefile << " " << std::to_string(modelNumber);
              timefile.close();
            }

            std::ofstream outfile;
            std::string filename = filenamebase + "_" + std::to_string(modelNumber) + "_" + std::to_string(traits.rank()) + ".data";

            outfile.open(filename, std::ios::out | std::ios::trunc);
            outfile << "electrode momentConcentration" << std::endl;
            for (auto const & electrode : current_output)
            {
              outfile << electrode.first << " " << electrode.second << std::endl;
            }
            outfile.close();
          
            printTimer.stop();
          }

          if (writeMomentsERT && isMomentsERT)
          {
            printTimer.start();
            
            // declare output map for this source electrode; [1:144]
            std::map<unsigned int, RF> current_output;
            unsigned int no_electrodes = extraction_helper<Storage>(current_output,storage,time); // get the results of this model

            unsigned int k; // moment number

            std::string str = modelname;
            std::string common_base = "momentsERT_";
            auto start_position_to_erase = str.find(common_base);
            str.erase(start_position_to_erase, common_base.size());      
            std::replace(str.begin(), str.end(), '_', ' ');


            std::stringstream ss; 
            ss << str;
            std::string temp;
            ss >> temp; // get the ERT model number
            ss >> temp; // get the moment number k
            std::stringstream(temp) >> k;

            // add the output for the active electrode to the map of all electrodes
            // increase modelNumber by one, as models are counted from 0, while electrode names start with 1

            std::pair<unsigned int, unsigned int> key;
            key = std::make_pair (k,modelNumber+1);
            output_MomentsERT_all_electrodes.insert(std::pair< std::pair<unsigned int, unsigned int>, std::map<unsigned int, RF> >(key,current_output) );

            unsigned int highest_moment = traits.config().template get<unsigned int>("moments.highest");
            // if this is the last ERT moments model, print output to file
            if (modelNumber+1 == no_electrodes && k == highest_moment)
            {
              for (unsigned int moment = 0; moment <= highest_moment ; ++moment)
              {
                std::string filenamebase = traits.config().template get<std::string>("output.writeMomentsERTFilename","resultsMomentsERT");
                if (traits.rank() == 0)
                {
                  std::cout << "printing ERT moment(" << moment <<") results" << std::endl;

                  // print the temporal information to a timefile
                  std::ofstream timefile;
                  std::string timefilename = filenamebase + ".moments";
                  if (clearFiles_MomentsERT)
                  { 
                    timefile.open(timefilename, std::ios::out | std::ios::trunc); // clear the file
                    clearFiles_MomentsERT = false;
                    timefile << "electrodes " << no_electrodes << std::endl;
                    timefile << "processors " << traits.helper.size() << std::endl;
                    timefile << "moments ";
                  } else 
                  {
                    timefile.open(timefilename, std::ios::app); // append to file
                  }
                  timefile << " " << std::to_string(moment);
                  timefile.close();
                }

                std::ofstream outfile;
                std::string filename = filenamebase + "_" + std::to_string(moment) + "_" + std::to_string(traits.rank()) + ".data";

                outfile.open(filename, std::ios::out | std::ios::trunc);
                outfile << "injected_el measured_el potential" << std::endl;

                // go through the collected data set -> get measured potentials for one injection
                for (auto const & entry : output_MomentsERT_all_electrodes)
                {
                  for (auto const & measurement_electrode : entry.second)
                  {
                    unsigned int injection_electrode = entry.first.second;
                    unsigned int meas_electrode = measurement_electrode.first;
                    RF value = measurement_electrode.second;
                    if (entry.first.first  == moment)
                      outfile << injection_electrode << " " << meas_electrode << " " << value << std::endl;
                  }
                }
                outfile.close();
              }
              output_MomentsERT_all_electrodes.clear(); // reset local storage of all ERT measurements
            }
            printTimer.stop();
          }
      }; 
    };
};

/**
 * @brief Tags for existing model types
 */
struct ModelTypes
{
  struct Groundwater {};
  struct Transport {};
  struct Geoelectrics {};
  struct Moments_c {};
  struct Moments_ERT {};
};


#endif // DUNE_FTG_MODELTRAITS_TRANSIENT_HH
