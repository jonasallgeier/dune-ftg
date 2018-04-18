// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FTG_HH
#define DUNE_FTG_HH
using namespace Dune::Modelling;

// electrode position was already obtained; 
// now, the coordinates of the electrodes are assigned to the cells they are located in
template<typename ModelTraits>
void set_electrodes(auto modelTraits)
{
  // get gridview, index set and hierarchic search object
  const typename ModelTraits::GridTraits::Grid::LeafGridView& mygridview = modelTraits->grid().leafGridView();
  auto& index_set = mygridview.indexSet();
  typedef Dune::HierarchicSearch<typename ModelTraits::GridTraits::Grid, decltype(index_set)> HierarchicSearch;
  HierarchicSearch hsearch(modelTraits->grid(), index_set);
  
  typename ModelTraits::GridTraits::Vector  electrode_coordinates;        // temporary vector containing coordinates of current electrode
  int no_electrodes = modelTraits->electrodeconfiguration.no_electrodes; // total number of electrodes
  auto electrode_cell_indices = std::vector<int>(no_electrodes,-1);            // vector for electrode cell indices, is later given to modelTraits
  
  for (int i = 0; i < no_electrodes; i++) // go through all electrodes and assign them  
  { 
    electrode_coordinates[0] = modelTraits->electrodeconfiguration.x_elec[i];  // get current coordinates
    electrode_coordinates[1] = modelTraits->electrodeconfiguration.y_elec[i];  // get current coordinates
    electrode_coordinates[2] = modelTraits->electrodeconfiguration.z_elec[i];  // get current coordinates
        

    try // in parallel, not every processor can assign all electrodes (-> domain decomposition)
    {
      const auto elem = hsearch.findEntity(electrode_coordinates);     // find out which cell contains the electrode coordinates
      electrode_cell_indices[i] = index_set.index(elem);               // assign the electrode coordinate  
      
      std::cout << "electrode coord... " << electrode_coordinates[0] << "; "<< electrode_coordinates[1] << "; "<< electrode_coordinates[2] << "(cell #"<< index_set.index(elem)<<")" <<std::endl;

      std::cout << "element center...  " << elem.geometry().center()[0] << "; "<< elem.geometry().center()[1] << "; "<< elem.geometry().center()[2] << "(cell #"<< index_set.index(elem)<<")" <<std::endl;

    }
    catch (Dune::Exception &e) {}
  } 
  // if there is an electrode that could not be localized, throw an exception! (check via. ini -1 cell index)
  // TODO This does not work in parallel... maybe make use of another method...
  /*std::vector<int>::iterator electrode_check = std::min_element(std::begin(electrode_cell_indices), std::end(electrode_cell_indices));
  if (electrode_cell_indices[std::distance(std::begin(electrode_cell_indices), electrode_check)] == -1 ) 
    {
      for (auto value : electrode_cell_indices) {std::cout << value << "; ";}
      std::cout << std::endl;
      DUNE_THROW( Dune::GridError, "There is an electrode outside of the grid." );
    }
  */
  modelTraits->set_electrode_cell_indices(electrode_cell_indices); // give electrode cell indices to modelTraits

}

// well position was already obtained; 
// now, the coordinates of the wells are assigned to the cells they are located in
template<typename ModelTraits>
void set_wells(auto modelTraits)
{
  const typename ModelTraits::GridTraits::Grid::LeafGridView& mygridview = modelTraits->grid().leafGridView();
  auto& index_set = mygridview.indexSet();
  
  typename ModelTraits::GridTraits::Vector  well_coordinates;       // temporary vector containing coordinates of current well
  int no_wells = modelTraits->wellconfiguration.no_wells;           // total number of wells
  std::vector<int> well_cells;          // vector for well cell indices, is later given to modelTraits
  std::vector<double> well_qs;          // vector for well pumping rates, is later given to modelTraits
  std::vector<int> well_identities;     // int vector representing well identities (a cell can have contributions from two different wells)
  std::vector<bool> well_injectionindicator; // true if cell contains a tracer injection; zero else

  std::vector<double> well_x = modelTraits->wellconfiguration.x_well;  
  std::vector<double> well_y = modelTraits->wellconfiguration.y_well;
  std::vector<double> well_z1 = modelTraits->wellconfiguration.z1_well;
  std::vector<double> well_z2 = modelTraits->wellconfiguration.z2_well;
  std::vector<double> well_q = modelTraits->wellconfiguration.q_well;
  std::vector<bool> well_in = modelTraits->wellconfiguration.injection_well;
  
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
      if (well_x[i]>=x_min && well_x[i]<x_max && well_y[i]>=y_min && well_y[i]<y_max)
      {
        if (well_z2[i] <= z_min || well_z1[i] >= z_max)
        {
          // well is outside of cell
        } else {
          // well is (partly) inside of cell
          well_cells.push_back(index_set.index(elem));
          well_qs.push_back(well_q[i]);
          well_identities.push_back(i);
          well_injectionindicator.push_back(well_in[i]);
         }
      }
    }
  }

  // distribute well injection rates equally across all participating cells
  for (int i = 0; i<no_wells;i++)
  {
    double tmp = count(well_identities.begin(), well_identities.end(), i);
    for(unsigned int j = 0; j != well_qs.size(); j++) 
    {
      if (well_identities[j]==i) {well_qs[j] /= tmp;}
    }
  }
  
  // sum all injection rate contributions in a cell up, if a single contribution injects tracer -> this is a tracer injection cell
  auto well_cell_indices = well_cells;
  std::sort(well_cell_indices.begin(), well_cell_indices.end());
  well_cell_indices.erase( std::unique( well_cell_indices.begin(), well_cell_indices.end() ), well_cell_indices.end() ); //get unique well cells
  
  std::vector<double> well_rates = std::vector<double>(well_cell_indices.size());
  std::vector<int> tracer_cell_indices;
  
  for (unsigned int i = 0; i!=well_cell_indices.size();i++)
  {
    well_rates[i]=0;
    for(unsigned int j = 0; j != well_qs.size(); j++) 
    {
      if (well_cells[j]==well_cell_indices[i]) 
      {
        well_rates[i] += well_qs[j];
        if (well_injectionindicator[j] == true) {tracer_cell_indices.push_back(well_cell_indices[i]);}
      }
    }
  }
  
  // get unique injection well cells (-> remove duplicate information)
  std::sort( tracer_cell_indices.begin(), tracer_cell_indices.end() );
  tracer_cell_indices.erase( std::unique( tracer_cell_indices.begin(), tracer_cell_indices.end() ), tracer_cell_indices.end() );


  modelTraits->set_well_cell_indices(well_cell_indices,well_rates,tracer_cell_indices); // give well cell indices to modelTraits

}





#endif // DUNE_FTG_HH
