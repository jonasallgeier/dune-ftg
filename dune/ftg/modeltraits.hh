// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_MODELLING_EXAMPLE_HH
#define DUNE_MODELLING_EXAMPLE_HH

#include<dune/grid/yaspgrid.hh>
#include<dune/randomfield/randomfield.hh>

/**
 * @brief Helper class for grid generation
 */
template<typename GT>
class GridHelper
{
  enum {dim = GT::dim};

  typedef typename GT::DomainField DF;
  //typedef typename GT::RangeField RF;

  int levels;
  std::vector<DF> maxExt;
  std::vector<DF> lowerleft, upperright;
  std::vector<unsigned int> minCells, maxCells;

  public:

  GridHelper(const Dune::ParameterTree& config)
  {
    levels      = config.get<int>                       ("grid.levels"    ,1);
    maxExt      = config.get<std::vector<DF> >          ("grid.extensions");    
    //lowerleft   = config.get<std::vector<DF> >          ("grid.lowerleft");   // get coordinates of lower left point
    //upperright  = config.get<std::vector<DF> >          ("grid.upperright"); // get coordinates of upper right point
    
    //maxExt = std::vector<DF>(3);
    //for (int i=0; i <3; i++)
    //{
    //  maxExt[i] = upperright[i]-lowerleft[i];
    //}    
    
    maxCells = config.get<std::vector<unsigned int> >("grid.cells");

    if (maxExt.size() != maxCells.size())
      DUNE_THROW(Dune::Exception,"cell and extension vectors differ in size");

    minCells = maxCells;
    for (int i = 0; i < levels - 1; i++)
      for (unsigned int j = 0; j < maxCells.size(); j++)
      {
        if (minCells[j]%2 != 0)
          DUNE_THROW(Dune::Exception,"cannot create enough levels for hierarchical grid, check number of cells");

        minCells[j] /= 2;
      }
  }

  Dune::FieldVector<DF,dim> L() const
  {
    Dune::FieldVector<DF,dim> Lvector;

    for (unsigned int i = 0; i < dim; i++)
      Lvector[i] = maxExt[i];

    return Lvector;
  }
  
  /*
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
  */

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
  std::vector<unsigned int> well_cell_indices;
  std::vector<double> well_rates;
  std::vector<unsigned int> tracer_cell_indices;

  public:
    // traits relating to geometry and underlying types
    struct GridTraits
    {
      enum {dim = dimension};
      using Grid     = Dune::YaspGrid<dim>;
//      using Grid     = Dune::YaspGrid<dim,Dune::EquidistantOffsetCoordinates<double, dim> >;
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
    struct MeasurementList
    {
      struct SubMeasurements
      {
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
        void extract(const Storage& storage, const RF& firstTime, const RF& lastTime)
        {
          std::cout << "(would extract measurements " << firstTime << " to " << lastTime << ")" << std::endl;
        }
      };

      std::shared_ptr<SubMeasurements> sub;
      const std::shared_ptr<SubMeasurements>& get(const std::string& name) const
      {
        return sub;
      }
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
      yaspGrid(gh.L(),gh.N(),gh.B(),1)
    {
    }
//    ModelTraits(Dune::MPIHelper& helper_, Dune::ParameterTree& duneConfig_)
//      : helper(helper_), duneConfig(duneConfig_), gh(duneConfig),
//      yaspGrid(gh.LL(),gh.UR(),gh.N(),gh.B(),1)
//    {
//    }



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

    // provide the vector of grid indices that contain wells
    std::vector<unsigned int> read_well_cell_indices() const
    {
      return well_cell_indices;
    }

    // provide the vector of pumping rates
    std::vector<double> read_well_rates() const
    {
      return well_rates;
    }

    // provide the vector of cell indices of tracer injection cells
    std::vector<unsigned int> read_tracer_cell_indices() const
    {
      return tracer_cell_indices;
    }

    // mechanism to define the vector of grid indices that contain wells
    void set_well_cell_indices(std::vector<unsigned int> in_well_cell_indices, std::vector<double> in_well_rates, std::vector<unsigned int> in_tracer_cell_indices)
    {
      well_cell_indices = in_well_cell_indices;
      well_rates = in_well_rates;
      tracer_cell_indices = in_tracer_cell_indices;
    }  

    // define electrode configuration; it is read in from a file by the constructor of this struct
    struct ElectrodeConfiguration
    {
        int no_electrodes;          // number of electrodes
        std::vector<int> name_elec; // name vector
        std::vector<double> x_elec; // vector of x coordinates
        std::vector<double> y_elec; // vector of y coordinates
        std::vector<double> z_elec; // vector of z coordinates
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
            double tmp_x;
            double tmp_y;
            double tmp_z;
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

    // define electrode configuration; it is read in from a file by the constructor of this struct
    struct WellConfiguration
    {
        int no_wells;                     // number of wells
        std::vector<double> x_well;       // vector of x coordinates
        std::vector<double> y_well;       // vector of y coordinates
        std::vector<double> z1_well;      // vector of z1 coordinates
        std::vector<double> z2_well;      // vector of z2 coordinates
        std::vector<double> q_well;       // vector of pumping rates
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
            file_wconf >> no_wells;  //first line equals number of electrodes
   
            //declare temporary storage variables
            double tmp_x;
            double tmp_y;
            double tmp_z1;
            double tmp_z2;
            double tmp_q;
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
    
    //call constructor to read in electrode configuration from file specified in .ini; data is stored in a struct object
    WellConfiguration wellconfiguration {duneConfig.get<std::string>("configfiles.wells"),rank()};    
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


#endif // DUNE_MODELLING_EXAMPLE_HH
