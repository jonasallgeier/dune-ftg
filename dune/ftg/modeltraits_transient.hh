// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_FTG_MODELTRAITS_TRANSIENT_HH
#define DUNE_FTG_MODELTRAITS_TRANSIENT_HH

#include<dune/grid/yaspgrid.hh>
// use slightly modified version of dune/randomfield.hh
#include<dune/ftg/override/randomfield.hh>
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
    isequidistant = config.get<bool>              ("grid.equidistant",1);
    levels      = config.get<int>                 ("grid.levels"    ,1);
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
          vec_i = config.get<std::vector<DF> >    ("grid.vector_x");
        else if (i==1)
          vec_i = config.get<std::vector<DF> >    ("grid.vector_y");
        else if (i==2)
          vec_i = config.get<std::vector<DF> >    ("grid.vector_z");
        else
          DUNE_THROW(Dune::Exception,"nonequidistant grid with more than four dimensions is not implemented");
        lowerleft.push_back(vec_i[0]);
        upperright.push_back(vec_i.back());
        coordinate_vectors[i] = vec_i;
      }
    }
  }

  //lower left point of grid
  Dune::FieldVector<DF,dim> LL() const
  {
    Dune::FieldVector<DF,dim> LLvector;

    for (unsigned int i = 0; i < dim; i++)
      LLvector[i] = lowerleft[i];

    return LLvector;
  }

  //upper right point of grid
  Dune::FieldVector<DF,dim> UR() const
  {
    Dune::FieldVector<DF,dim> URvector;

    for (unsigned int i = 0; i < dim; i++)
      URvector[i] = upperright[i];

    return URvector;
  }

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
      using Grid     = Dune::YaspGrid<dim,Dune::TensorProductCoordinates<DF, dim> >;
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



    // dummy implementation
    struct SensitivityList
    {

    };

    // use a list of random Gaussian fields as parameter list
    using ParameterList = Dune::RandomField::RandomFieldList<GridTraits>;

    const Dune::MPIHelper& helper;
    const Dune::ParameterTree& duneConfig;
    const GridHelper<GridTraits> gh;

    typename GridTraits::Grid yaspGrid;

    ModelTraits(Dune::MPIHelper& helper_, Dune::ParameterTree& duneConfig_)
      : helper(helper_), duneConfig(duneConfig_), gh(duneConfig),
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

    // mechanism to define the vector of grid indices that contain electrodes
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
        std::vector<RF> x_elec; // vector of x coordinates
        std::vector<RF> y_elec; // vector of y coordinates
        std::vector<RF> z_elec; // vector of z coordinates
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
            file_econf >> no_electrodes;  //first line equals number of electrodes
   
            //declare temporary storage variables
            int tmp_name;
            RF tmp_x;
            RF tmp_y;
            RF tmp_z;
            int tmp_surf;
        
            //read data in lines {#,x,y,z,s} from file via temp to vector 
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
            
            //close the file and give console output
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
    
    //call constructor to read in electrode configuration from file specified in .ini; data is stored in a struct object
    ElectrodeConfiguration electrodeconfiguration {duneConfig.get<std::string>("configfiles.electrodes"),rank()};

    // define well configuration; it is read in from a file by the constructor of this struct
    struct WellConfiguration
    {
        int no_wells;                     // number of wells
        std::vector<RF> x_well;       // vector of x coordinates
        std::vector<RF> y_well;       // vector of y coordinates
        std::vector<RF> z1_well;      // vector of z1 coordinates
        std::vector<RF> z2_well;      // vector of z2 coordinates
        std::vector<RF> q_well;       // vector of pumping rates
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
            file_wconf >> no_wells;  //first line equals number of wells
   
            //declare temporary storage variables
            RF tmp_x;
            RF tmp_y;
            RF tmp_z1;
            RF tmp_z2;
            RF tmp_q;
            bool tmp_in;
        
            //read data in lines {#,x,y,z,s} from file via temp to vector 
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
            
            //close the file and give console output
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
    
    //call constructor to read in well configuration from file specified in .ini; data is stored in a struct object
    WellConfiguration wellconfiguration {duneConfig.get<std::string>("configfiles.wells"),rank()};

    class MeasurementList
    {
      private:
        const ModelTraits& traits;
        bool clearFiles = true;
        std::map<unsigned int, std::map<unsigned int, RF> > output_all_electrodes;
        
        auto& index_set() const
        {
          return traits.grid().leafGridView().indexSet();
        };

      public:
      MeasurementList(const ModelTraits& traits_)
        : traits(traits_)
      {}

      //struct SubMeasurements
      //{

        
        //bool printToFile = false;

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
        void extract(const Storage& storage, const RF& firstTime, const RF& time, const unsigned int& modelNumber,const std::string& timeString, Dune::Timer&  printTimer)              
        {
          std::vector<unsigned int> electrode_cells = traits.read_electrode_cell_indices();
          const typename GridTraits::Domain x = {0.5, 0.5, 0.5}; // get value at cell center

          //declare output map for this source electrode; [1:144]
          std::map<unsigned int, RF> output_this_electrode;
      
          for (const auto & elem : elements (traits.grid().leafGridView()))
          {
            //check if this cell is a measurement cell/is multiple measurement cells --> get electrode numbers of current cell
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
          output_all_electrodes.insert(std::pair<unsigned int,std::map<unsigned int, RF> >(modelNumber+1,output_this_electrode));

          // if this is the last geoelectrics model print output to file
          if (modelNumber+1 == electrode_cells.size())
          {
            printTimer.start();
            std::string filenamebase = traits.config().template get<std::string>("output.writeGeoelectricsFilename","results");
            if (traits.rank() == 0)
            {
              std::cout << "printing ERT results" << std::endl;

              // print the temporal information to a timefile
              std::ofstream timefile;
              std::string timefilename;
              timefilename.append(filenamebase);
              timefilename.append(".times");

              if (clearFiles)
              { 
                timefile.open(timefilename, std::ios::out | std::ios::trunc); // clear the file
                clearFiles = false;
                timefile << "electrodes " << electrode_cells.size() << std::endl;
                timefile << "processors " << traits.helper.size() << std::endl;
                timefile << "times";
              } else 
              {
                timefile.open(timefilename, std::ios::app); // append to file
              }
              timefile << " " << time;
              timefile.close();
            }

            std::stringstream ss;
            ss << traits.rank();
            std::string rank(ss.str());

            std::ofstream outfile;
            std::string filename;
            filename.append(filenamebase);
            filename.append("_");
            filename.append(timeString);
            filename.append("_");
            filename.append(rank);
            filename.append(".data");

            outfile.open(filename, std::ios::out | std::ios::trunc);

            outfile << "injected_el measured_el potential" << std::endl;

            for (auto const & injection_electrode : output_all_electrodes)
            {
              for (auto const & measurement_electrode : injection_electrode.second)
              {
                outfile << injection_electrode.first << " " << measurement_electrode.first << " " << measurement_electrode.second << std::endl;
              }
            }
            outfile.close();
            output_all_electrodes.clear(); // reset local storage of all ERT measurements
          
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
};


#endif // DUNE_FTG_MODELTRAITS_TRANSIENT_HH