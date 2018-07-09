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
  auto electrode_cell_indices = std::vector<unsigned int>(no_electrodes,-1);            // vector for electrode cell indices, is later given to modelTraits
  
  for (int i = 0; i < no_electrodes; i++) // go through all electrodes and assign them  
  { 
    electrode_coordinates[0] = modelTraits->electrodeconfiguration.x_elec[i];  // get current coordinates
    electrode_coordinates[1] = modelTraits->electrodeconfiguration.y_elec[i];  // get current coordinates
    electrode_coordinates[2] = modelTraits->electrodeconfiguration.z_elec[i];  // get current coordinates
        

    try // in parallel, not every processor can assign all electrodes (-> domain decomposition)
    {
      const auto elem = hsearch.findEntity(electrode_coordinates);     // find out which cell contains the electrode coordinates
      electrode_cell_indices[i] = index_set.index(elem);               // assign the electrode coordinate
    }
    catch (Dune::Exception &e) {}
  }
  modelTraits->set_electrode_cell_indices(electrode_cell_indices); // give electrode cell indices to modelTraits

}

// well position was already obtained; 
// now, the coordinates of the wells are assigned to the cells they are located in
template<typename ModelTraits>
void set_wells(auto modelTraits)
{
  const typename ModelTraits::GridTraits::Grid::LeafGridView& mygridview = modelTraits->grid().leafGridView();
  using RF = typename ModelTraits::GridTraits::RangeField;
  auto& index_set = mygridview.indexSet();
  
  typename ModelTraits::GridTraits::Vector  well_coordinates;       // temporary vector containing coordinates of current well
  int no_wells = modelTraits->wellconfiguration.no_wells;           // total number of wells
  std::vector<unsigned int> well_cells_temp;          // vector for well cell indices, is later given to modelTraits
  std::vector<RF> well_qs;          // vector for well pumping rates, is later given to modelTraits
  std::vector<int> well_identities;     // int vector representing well identities (a cell can have contributions from two different wells)
  std::vector<bool> well_injectionindicator; // true if cell contains a tracer injection; zero else

  std::vector<RF> well_x = modelTraits->wellconfiguration.x_well;  
  std::vector<RF> well_y = modelTraits->wellconfiguration.y_well;
  std::vector<RF> well_z1 = modelTraits->wellconfiguration.z1_well;
  std::vector<RF> well_z2 = modelTraits->wellconfiguration.z2_well;
  std::vector<RF> well_q = modelTraits->wellconfiguration.q_well;
  std::vector<bool> well_in = modelTraits->wellconfiguration.injection_well;
  
  // go through all cells and determine the well contributions
  for (const auto & elem : elements (mygridview))
  {
    // get borders of cuboid cell from corner coordinates
    RF x_min = elem.geometry().corner(0)[0];
    RF x_max = elem.geometry().corner(1)[0];
    RF y_min = elem.geometry().corner(0)[1];
    RF y_max = elem.geometry().corner(2)[1];
    RF z_min = elem.geometry().corner(0)[2];
    RF z_max = elem.geometry().corner(4)[2];
    
    // go through all wells and check if they contribute to this cell
    for (int i = 0; i<no_wells; i++)
    {
      if (well_x[i]>=x_min && well_x[i]<x_max && well_y[i]>=y_min && well_y[i]<y_max)
      {
        if (!((well_z2[i] <= z_min) || (well_z1[i] >= z_max)))
        {
          // well is (at least partly) inside of cell
          well_cells_temp.push_back(index_set.index(elem));
          well_identities.push_back(i);
          well_injectionindicator.push_back(well_in[i]);
          RF local_q;
          RF V_well = (x_max-x_min)*(y_max-y_min)*(well_z2[i]-well_z1[i]);

          RF local_max = std::min(z_max,well_z2[i]);
          RF local_min = std::max(z_min,well_z1[i]);
          RF h = local_max-local_min;
          RF V_cell = (x_max-x_min)*(y_max-y_min)*h;

          local_q = well_q[i]*V_cell/V_well;
          well_qs.push_back(local_q);
        }
      }
    }
  }

/*
  // distribute well injection rates equally across all participating cells
  for (int i = 0; i<no_wells;i++)
  {
    RF tmp = count(well_identities.begin(), well_identities.end(), i);
    for(unsigned int j = 0; j != well_qs.size(); j++) 
    {
      if (well_identities[j]==i) {well_qs[j] /= tmp;}
    }
  }
*/  
  // sum all injection rate contributions in a cell up, if a single contribution injects tracer -> this is a tracer injection cell
  std::vector<unsigned int> well_cell_indices = well_cells_temp;
  
  std::sort(well_cell_indices.begin(), well_cell_indices.end());
  well_cell_indices.erase( std::unique( well_cell_indices.begin(), well_cell_indices.end() ), well_cell_indices.end() ); //get unique well cells
  
  std::vector<RF> well_rates = std::vector<RF>(well_cell_indices.size());
  std::map<unsigned int, std::pair<RF, bool>> well_cells; // cell_indices & pumping rate

  for (unsigned int i = 0; i!=well_cell_indices.size();i++)
  {
    well_rates[i]=0;
    bool injection_cell = false;
    for(unsigned int j = 0; j != well_qs.size(); j++) 
    {
      if (well_cells_temp[j]==well_cell_indices[i]) 
      {
        well_rates[i] += well_qs[j];
        if (well_injectionindicator[j] == true) {injection_cell = true;}
      }
    }
    well_cells.insert(std::make_pair(well_cell_indices[i],std::make_pair(well_rates[i], injection_cell) ) );
  }
  modelTraits->set_well_cells(well_cells); // give well cell indices to modelTraits
}


template<typename ModelTraits>
void unify_ERT_results(auto modelTraits)
{
  std::string filenamebase = modelTraits->config().template get<std::string>("output.writeERTFilename","resultsERT");
  if (modelTraits->basePotentialEvaluation)
  {
    filenamebase = modelTraits->config().template get<std::string>("output.writeERTBaseFilename","resultsERTBase");
  }

  std::string timefilename = filenamebase + ".times";
  
  std::string temp;
  std::ifstream timefile(timefilename);
  unsigned int no_electrodes;
  unsigned int no_processors;
  timefile >> temp; // skip label
  timefile >> no_electrodes;
  timefile >> temp; // skip label
  timefile >> no_processors;
  timefile >> temp; // skip label

  std::string time;
  std::vector<std::string> times;
  while ( !timefile.eof() )
  {
    timefile >> time;
    times.push_back(time);
  }
  using RF  = typename ModelTraits::GridTraits::RangeField;

  for (auto const& timeString: times) 
  {
    std::map<std::pair<unsigned int, unsigned int>, RF > complete_map; // <injection electrode, <measured electrode, potential> > 
  
    for( unsigned int proc = 0; proc < no_processors; ++proc ) 
    {
      std::string infilename = filenamebase;

      if (modelTraits->basePotentialEvaluation)
      {
        infilename += ("_base_" + std::to_string(proc) + ".data");
      } else
      {
        infilename += ("_" + timeString + "_" + std::to_string(proc) + ".data");
      }

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

        std::pair<unsigned int, unsigned int> key;
        key = std::make_pair(injection_electrode,measured_electrode);

        complete_map.insert(std::pair<std::pair<unsigned int, unsigned int>, RF >(key, potential));
      }
      infile.close();
      remove(infilename.c_str());
      
      std::ofstream outfile;
      std::string outfilename = filenamebase;

      if (modelTraits->basePotentialEvaluation)
      {
        outfilename += (".data");
      } else
      {
        outfilename += ("_" + timeString + ".data");
      }

      outfile.open(outfilename, std::ios::out | std::ios::trunc);
      outfile << "injected_el measured_el potential" << std::endl;
      for (auto const & current_entry : complete_map)
      {
        outfile << current_entry.first.first << " " << current_entry.first.second << " " << current_entry.second << std::endl;
      }
      outfile.close();
    }
  }
}

template<typename ModelTraits>
void unify_transport_results(auto modelTraits)
{
  std::string filenamebase = modelTraits->config().template get<std::string>("output.writeTransportFilename","resultsTransport");
  
  std::string timefilename = filenamebase + ".times";
 
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
  using  RF  = typename ModelTraits::GridTraits::RangeField;

  for (auto const& timeString: times) 
  {
    std::map<unsigned int, RF > complete_map; // < electrode, concentration > 
  
    for( unsigned int proc = 0; proc < no_processors; ++proc ) 
    {
      std::string infilename = filenamebase + "_" + timeString + "_" + std::to_string(proc) + ".data";
      std::ifstream infile(infilename);

      // skip header line
      infile >> temp;
      infile >> temp;

      while ( !infile.eof() )
      {
        unsigned int electrode;
        RF concentration;
        infile >> electrode;
        infile >> concentration;
        complete_map.insert(std::pair<unsigned int, RF>(electrode, concentration));
      }
      infile.close();
      remove(infilename.c_str());
      
      std::ofstream outfile;
      std::string outfilename = filenamebase + "_" + timeString + ".data";
      outfile.open(outfilename, std::ios::out | std::ios::trunc);

      outfile << "electrode concentration" << std::endl;
      for (auto const & current_entry : complete_map)
      {
        outfile << current_entry.first << " " << current_entry.second << std::endl;
      }
      outfile.close();
    }
  }
}

template<typename ModelTraits>
void unify_momentsTransport_results(auto modelTraits)
{
  std::string filenamebase = modelTraits->config().template get<std::string>("output.writeMomentsTransportFilename","resultsMomentsTransport");
  std::string timefilename = filenamebase + ".moments";
 
  std::string temp;
  std::ifstream timefile(timefilename);
  unsigned int no_electrodes;
  unsigned int no_processors;
  timefile >> temp;
  timefile >> no_electrodes;
  timefile >> temp;
  timefile >> no_processors;
  timefile >> temp;

  std::string moment;
  std::vector<std::string> moments;
  while ( !timefile.eof() )
  {
    timefile >> moment;
    moments.push_back(moment);
  }
  using  RF  = typename ModelTraits::GridTraits::RangeField;

  for (auto const& momentString: moments) 
  {
    std::map<unsigned int, RF > complete_map; // < electrode, concentration > 
  
    for( unsigned int proc = 0; proc < no_processors; ++proc ) 
    {
      std::string infilename = filenamebase + "_" + momentString + "_" + std::to_string(proc) + ".data";
      std::ifstream infile(infilename);

      // skip header line
      infile >> temp;
      infile >> temp;

      while ( !infile.eof() )
      {
        unsigned int electrode;
        RF momentTransport;
        infile >> electrode;
        infile >> momentTransport;
        complete_map.insert(std::pair<unsigned int, RF>(electrode, momentTransport));
      }
      infile.close();
      remove(infilename.c_str());
      
      std::ofstream outfile;
      std::string outfilename = filenamebase + "_" + momentString + ".data";
      outfile.open(outfilename, std::ios::out | std::ios::trunc);

      outfile << "electrode momentTransport" << std::endl;
      for (auto const & current_entry : complete_map)
      {
        outfile << current_entry.first << " " << current_entry.second << std::endl;
      }
      outfile.close();
    }
  }
}

template<typename ModelTraits>
void unify_momentsERT_results(auto modelTraits)
{
  std::string filenamebase = modelTraits->config().template get<std::string>("output.writeMomentsERTFilename","resultsMomentsERT");
  
  std::string timefilename = filenamebase + ".moments";
  
  std::string temp;
  std::ifstream timefile(timefilename);
  unsigned int no_electrodes;
  unsigned int no_processors;
  timefile >> temp; // skip label
  timefile >> no_electrodes;
  timefile >> temp; // skip label
  timefile >> no_processors;
  timefile >> temp; // skip label

  std::string moment;
  std::vector<std::string> moments;
  while ( !timefile.eof() )
  {
    timefile >> moment;
    moments.push_back(moment);
  }
  using RF  = typename ModelTraits::GridTraits::RangeField;

  for (auto const& momentString: moments) 
  {
    std::map<std::pair<unsigned int, unsigned int>, RF > complete_map; // <injection electrode, <measured electrode, potential> > 
  
    for( unsigned int proc = 0; proc < no_processors; ++proc ) 
    {
      std::string infilename = filenamebase;

      infilename += ("_" + momentString + "_" + std::to_string(proc) + ".data");

      std::ifstream infile(infilename);
      //skip header line
      infile >> temp;
      infile >> temp;
      infile >> temp;

      unsigned int injection_electrode;
      unsigned int measured_electrode;
      RF value;

      while ( !infile.eof() )
      {
        infile >> injection_electrode;
        infile >> measured_electrode;
        infile >> value;

        std::pair<unsigned int, unsigned int> key;
        key = std::make_pair(injection_electrode,measured_electrode);

        complete_map.insert(std::pair<std::pair<unsigned int, unsigned int>, RF >(key, value));
      }
      infile.close();
      remove(infilename.c_str());
      
      std::ofstream outfile;
      std::string outfilename = filenamebase;

      outfilename += ("_" + momentString + ".data");
      outfile.open(outfilename, std::ios::out | std::ios::trunc);
      outfile << "injected_el measured_el potential" << std::endl;
      for (auto const & current_entry : complete_map)
      {
        outfile << current_entry.first.first << " " << current_entry.first.second << " " << current_entry.second << std::endl;
      }
      outfile.close();
    }
  }
}

#endif // DUNE_FTG_HH
