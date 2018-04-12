#ifndef DUNE_FTG_HH
#define DUNE_FTG_HH
using namespace Dune::Modelling;

template<typename ModelTraits>
void setelectrodes(auto modelTraits)
{
  // electrode position was already obtained; now, the coordinates of the electrodes are assigned to the cells they are located in
  const typename ModelTraits::GridTraits::Grid::LeafGridView& mygridview = modelTraits->grid().leafGridView(); // get leaf grid view
  auto& index_set = mygridview.indexSet();                              //get index set
  typedef Dune::HierarchicSearch<typename ModelTraits::GridTraits::Grid, decltype(index_set)> HierarchicSearch;      // define how a hierarchic search object  looks like
  HierarchicSearch hsearch(modelTraits->grid(), index_set);              // define a hierarchic search object
  
  typename ModelTraits::GridTraits::Vector  electrode_coordinates;      // temporary vector containing coordinates of current electrode
  int no_electrodes = modelTraits->electrodeconfiguration.no_electrodes; // total number of electrodes
  auto electrode_cells = std::vector<int>(no_electrodes,-1);               // define a vector for electrode cell indices, is later given to modelTraits
  
  for (int i = 0; i < no_electrodes; i++) // go through all electrodes and assign them  
  { 
    electrode_coordinates[0] = modelTraits->electrodeconfiguration.x_elec[i];  // get current coordinates
    electrode_coordinates[1] = modelTraits->electrodeconfiguration.y_elec[i];  // get current coordinates
    electrode_coordinates[2] = modelTraits->electrodeconfiguration.z_elec[i];  // get current coordinates
    // throw exception, when a single electrode is outside of domain    
    try
    {
      const auto test_answer = hsearch.findEntity(electrode_coordinates);     // find out which cell contains the electrode coordinates
      electrode_cells[i] = index_set.index(test_answer);                      // assign the electrode coordinate  
      //std::cout << "electrode " << i+1 << " of " << no_electrodes << " electrodes: "<<  electrode_coordinates << " -> " << electrode_cells[i] << std::endl; // get the cell and give output  
    }
    catch (Dune::Exception &e)
    {
      //std::cerr << "There is an electrode outside of the grid. This is not allowed!" << std::endl;
      //throw;
      //if (mygridview.comm().rank()==0)
      //{ 
      //std::cerr << "electrode could not be located!" << std::endl;
      //}
    }
    // if there is an electrode that could not be localized, throw an exception!
    std::vector<int>::iterator electrode_check = std::min_element(std::begin(electrode_cells), std::end(electrode_cells));
    if (electrode_cells[std::distance(std::begin(electrode_cells), electrode_check)] == -1 ) { DUNE_THROW( Dune::GridError, "There is an electrode outside of the grid." );}
  }
  modelTraits->set_electrodecells(electrode_cells); // give electrode cell indices to modelTraits

}


template<typename ModelTraits>
void setwells(auto modelTraits)
{
  const typename ModelTraits::GridTraits::Grid::LeafGridView& mygridview = modelTraits->grid().leafGridView(); // get leaf grid view
  auto& index_set = mygridview.indexSet();                              //get index set
  
  typename ModelTraits::GridTraits::Vector  well_coordinates;       // temporary vector containing coordinates of current well
  int no_wells = modelTraits->wellconfiguration.no_wells;            // total number of wells
  std::vector<int> well_cells;         // define a vector for well cell indices, is later given to modelTraits
  std::vector<double> well_qs;
  std::vector<int> well_identities;  
  std::vector<double> well_x = modelTraits->wellconfiguration.x_well;  
  std::vector<double> well_y = modelTraits->wellconfiguration.y_well;
  std::vector<double> well_z1 = modelTraits->wellconfiguration.z1_well;
  std::vector<double> well_z2 = modelTraits->wellconfiguration.z2_well;
  std::vector<double> well_q = modelTraits->wellconfiguration.q_well;
  
  // go through all cells and determine the well contributions
  for (const auto & elem : elements (mygridview))
  {
    // get borders of cuboid cell from corner coordinates
    double x_min = elem.geometry().corner(0)[0];
    double x_max = elem.geometry().corner(1)[0];
    double y_min = elem.geometry().corner(0)[1];
    double y_max = elem.geometry().corner(2)[1];
    double z_min = elem.geometry().corner(0)[2];
    double z_max = elem.geometry().corner(4)[2];
    
    // go through all wells and check if they contribute to this cell
    for (int i = 0; i<no_wells; i++)
    {
      //double tmp_q = 0;
      if (well_x[i]>=x_min && well_x[i]<x_max && well_y[i]>=y_min && well_y[i]<y_max)
      {
        if (well_z2[i] <= z_min || well_z1[i] > z_max)
        {
          // well is outside of cell
        } else {
          // well is inside of cell
          well_cells.push_back(index_set.index(elem));
          well_qs.push_back(well_q[i]);
          well_identities.push_back(i);
        }
      }
    }
  }

  // distribute well inflow across all participating cells
  for (int i = 0; i<no_wells;i++)
  {
    double tmp = count(well_identities.begin(), well_identities.end(), i);
    for(unsigned int j = 0; j != well_qs.size(); j++) 
    {
      if (well_identities[j]==i) {well_qs[j] /= tmp;}
    }
  }
  
  // sum all distributions in a cell up
  auto well_cells_new = well_cells;
  std::sort( well_cells_new.begin(), well_cells_new.end() );
  well_cells_new.erase( std::unique( well_cells_new.begin(), well_cells_new.end() ), well_cells_new.end() );
  
  std::vector<double> well_qs_new = std::vector<double>(well_cells_new.size());
  
  for (unsigned int i = 0; i!=well_cells_new.size();i++)
  {
    well_qs_new[i]=0;
    for(unsigned int j = 0; j != well_qs.size(); j++) 
    {
      if (well_cells[j]==well_cells_new[i]) {well_qs_new[i] += well_qs[j];}
    }
  }
  modelTraits->set_wellcells(well_cells_new,well_qs_new); // give well cell indices to modelTraits

}





#endif // DUNE_FTG_HH
