#include <iostream>
#include <unordered_map>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <iterator>
#include <iomanip>
#include <limits>
#include <chrono>

#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/JetsWithoutJets.hh"
#include "fastjet/contrib/VariableRPlugin.hh"

#include "fastjet/PseudoJet.hh"

#include "../interface/Event.h"
#include "../interface/Property.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/tools/Filter.hh"


using namespace std;
using namespace fastjet;
using namespace contrib;



void analyze_event(MOD::Event & event_being_read, ofstream & output_file, int & event_serial_number, unordered_map<int, pair <int, double>> &lumi_info_by_lumi_block, unordered_map<int, pair <int, double>> &lumi_info_by_lumi_block2);

double angularity_lambda(PseudoJet jet, double jet_radius, float k, float beta);
double pT_D(PseudoJet jet);
unordered_map<int, pair <int, double>> get_lumi_info_lumi_block(string & input_file, int number_of_events_to_process);
unordered_map<int, pair <int, double>> get_lumi_info_lumi_block2(string & input_file, int number_of_events_to_process);
double highest_this_pT;
double lowest_this_pT;

int main(int argc, char * argv[]) {

    auto start = std::chrono::steady_clock::now();
    
    int number_of_events_to_process;
    
    if (argc <= 2) {
        std::cerr << "ERROR: You need to supply five arguments- first, path to the input data; second, path to the output file; third, number of events to process. The path has to be either absolute or relative to the bin directory:" << std::endl << std::endl << "./analysis (input_file.dat) (output_file.dat) [optional Nev]" << std::endl;
        return 1;
    }
    else if (argc == 3) {
        // Fifth argument is missing, process everything.
        number_of_events_to_process = std::numeric_limits<int>::max();
    }
    else {
        // Fifth argument gives the number of events to process.
        number_of_events_to_process = stoi(argv[3]);
    }
    
    ifstream data_file(argv[1]);
    ofstream output_file(argv[2], ios::out | ios::app);
    
    
    cout << endl << endl << "Starting analysis with the following given arguments: " << endl;
    cout << "Input file: " << argv[1] << endl;
    cout << "Output file: " << argv[2] << endl;
    cout << "Number of events2: ";
    
    
    
    if(argc == 3)
        cout << "ALL" << "  " << number_of_events_to_process << endl << endl;
    else
        cout << number_of_events_to_process << endl << endl;
    
    MOD::Event event_being_read;
    
    int event_serial_number = 1;
    unordered_map<int, pair <int, double>> lumi_info_by_lumi_block;
    unordered_map<int, pair <int, double>> lumi_info_by_lumi_block2;

    
    string input_file(argv[1]);
    
    lumi_info_by_lumi_block = get_lumi_info_lumi_block(input_file, number_of_events_to_process);
    
    cout << sizeof(lumi_info_by_lumi_block) << endl;

    
    lumi_info_by_lumi_block2 = get_lumi_info_lumi_block2(input_file, number_of_events_to_process);
    
    cout << sizeof(lumi_info_by_lumi_block2) << endl;
    
    
    while( event_being_read.read_event(data_file) && ( event_serial_number <= number_of_events_to_process ) ) {
        
        if( (event_serial_number % 5000) == 0 )
            cout << "Processing event number " << event_serial_number << endl;
        
        analyze_event(event_being_read, output_file, event_serial_number,lumi_info_by_lumi_block, lumi_info_by_lumi_block2 );
        
        
        event_being_read = MOD::Event();
        
        event_serial_number++;
        
    }
    


   auto finish = std::chrono::steady_clock::now();
   double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start).count();
   cout << "Finished processing " << (event_serial_number - 1) << " events in " << elapsed_seconds << " seconds!" << endl;

   return 0;
}


void analyze_event(MOD::Event & event_being_read, ofstream & output_file, int & event_serial_number, unordered_map<int, pair <int, double>> &lumi_info_by_lumi_block, unordered_map<int, pair <int, double>> & lumi_info_by_lumi_block2) {

   
   
  // cout << "Event number is " << event_being_read.event_number() << endl;
    double jet_radius = 0.5;
    
    JetDefinition jet_def_cambridge(cambridge_algorithm, fastjet::JetDefinition::max_allowable_R);
    
    Selector pT_1_GeV_selector = SelectorPtMin(1.0);
    Selector pT_0_5_GeV_selector = SelectorPtMin(0.5);
    Selector pT_015_020_GeV_selector = SelectorPtRange(0.15, 0.2);
    

    
   

   vector<PseudoJet> hardest_jet_constituents = pT_1_GeV_selector(event_being_read.hardest_jet().constituents());
   vector<PseudoJet> filtered_hardest_jet_constituents = pT_015_020_GeV_selector(event_being_read.hardest_jet().constituents());

   // vector<PseudoJet> hardest_jet_constituents = event_being_read.hardest_jet().constituents();
    

   ClusterSequence cs = ClusterSequence(hardest_jet_constituents, jet_def_cambridge);
   
   if (cs.inclusive_jets().size() == 0) {
      return;
   }

   // This is the hardest FastJet jet. This does not have JEC applied to it. We need to go through this mess because we want to apply the pT selector to the hardest jet's constituents and then use the corresponding hardest jet.
   // But since JEC does not apply to individual PFCs, we can't directly use event_being_read.hardest_jet().constituents().

   ClusterSequence cs_uncorrected_jet_no_softkiller = ClusterSequence(event_being_read.hardest_jet().constituents(), jet_def_cambridge);
   PseudoJet uncorrected_hardest_jet_no_softkiller = cs_uncorrected_jet_no_softkiller.inclusive_jets()[0];


   PseudoJet uncorrected_hardest_jet_with_softkiller = cs.inclusive_jets()[0];
   

   // PseudoJet hardest_jet = uncorrected_hardest_jet_with_softkiller * jec;
   PseudoJet hardest_jet = uncorrected_hardest_jet_with_softkiller;
   
    
   double beta = 1.5;

   SoftDrop soft_drop(beta, 0.1);
   PseudoJet soft_drop_jet = soft_drop(hardest_jet);
    
    
    
    
    
    int lumi_block = event_being_read.condition().lumi_block();
    int run_number = event_being_read.condition().run_number();
    int hash_value = run_number*1000+lumi_block;
    
    pair <int, double> number_of_events_and_lumi;
    pair <int, double> number_of_events_and_lumi2;


    number_of_events_and_lumi = lumi_info_by_lumi_block[hash_value];
    
    number_of_events_and_lumi2 = lumi_info_by_lumi_block2[hash_value];
    
    
    //cout << "lumi block " << event_being_read.condition() << endl;
    
    //cout << "hash_value" << hash_value << " " << number_of_events_and_lumi2.first << endl;
   // cout << number_of_events_and_lumi.first << " " << number_of_events_and_lumi.second << "numebr of events" << endl;
    
   vector<MOD::Property> properties;

    double jec = 1.0;
    
    if ((event_being_read.data_source() == 3) || (event_being_read.data_source() == 0))
        jec = event_being_read.get_hardest_jet_jec();

    

   properties.push_back(MOD::Property("# Entry", "  Entry"));
    
   //cout << "second" << number_of_events_and_lumi.second << "and integrated" << event_being_read.condition().


   properties.push_back(MOD::Property("prescale", event_being_read.weight()));
    
   
   // multiplicity, jet_multiplicity
   properties.push_back(MOD::Property("multiplicity", (int) event_being_read.particles().size()));
    
    
   properties.push_back(MOD::Property("crosssection", number_of_events_and_lumi2.first*1.0/number_of_events_and_lumi.second));
    
   properties.push_back(MOD::Property("trigger_fired", event_being_read.assigned_trigger_name()));

    
   //cout << number_of_events_and_lumi.second << " " << event_being_read.condition().integrated_recorded_lumi() << endl;
    

   properties.push_back(MOD::Property("events_being_read", event_being_read.weight()));
   properties.push_back(MOD::Property("hardest_pT", jec * event_being_read.hardest_jet().pt()));
    
    if(event_being_read.assigned_trigger_name() == "HLT_Jet150"){
    
   cout << jec * event_being_read.hardest_jet().pt() << endl;
   //cout << "trigger: " << event_being_read.assigned_trigger_name() << endl;
    }

   properties.push_back(MOD::Property("jec", jec));
    
   // cout << " jec " << jec << endl;

   properties.push_back( MOD::Property("hardest_eta", event_being_read.hardest_jet().eta()) );


    
 
   
   properties.push_back( MOD::Property("softkill_pT_loss", (uncorrected_hardest_jet_no_softkiller.pt() - uncorrected_hardest_jet_with_softkiller.pt() ) / uncorrected_hardest_jet_no_softkiller.pt() ) );
   properties.push_back( MOD::Property("frac_pT_loss", (hardest_jet.pt() - soft_drop( hardest_jet ).pt() ) / hardest_jet.pt() ) );
   properties.push_back( MOD::Property("hardest_eta", event_being_read.hardest_jet().eta()) );   



   properties.push_back( MOD::Property("hardest_phi", event_being_read.hardest_jet().phi()) );

   if ((event_being_read.data_source() == 3) || (event_being_read.data_source() == 0))
      properties.push_back( MOD::Property("hardest_area", event_being_read.get_hardest_jet_area()) );
   else
      properties.push_back( MOD::Property("hardest_area", 0.0) );


   vector<pair<string, double>> zg_cuts { make_pair("05", 0.05), make_pair("10", 0.1), make_pair("20", 0.2) };
    
/*
   // zg, dr, and mu for zg_cuts of 0.05, 0.1 and 0.2.

   for (unsigned i = 0; i < zg_cuts.size(); i++) {

      string label = zg_cuts[i].first;
      double zg_cut = zg_cuts[i].second;

      SoftDrop soft_drop(beta, zg_cut);
   
      PseudoJet soft_drop_jet = soft_drop(hardest_jet);

      double zg = soft_drop_jet.structure_of<SoftDrop>().symmetry();
      double rg = soft_drop_jet.structure_of<SoftDrop>().delta_R() / 0.5;

      properties.push_back(MOD::Property("zg_" + label, zg));
      properties.push_back(MOD::Property("mu_" + label, soft_drop_jet.structure_of<SoftDrop>().mu()));

      properties.push_back(MOD::Property("rg_" + label, rg));
      properties.push_back(MOD::Property("e1_" + label, rg * zg));
      properties.push_back(MOD::Property("e2_" + label, pow(rg, 2) * zg ));
      properties.push_back(MOD::Property("e05_" + label, sqrt(rg) * zg ));
   }
    
    */
    // H_T(summed scalar transverse momentum) test
    
    // ShapeJetMultiplicity Nj(jet_radius, 200.0);
    //cout << "N_jet=" << Nj(hardest_jet_constituents) << endl;

    

   // Analysis related to the effects of SoftDrop- observables before and after SoftDrop.
   
   properties.push_back(MOD::Property("pT_after_SD", soft_drop_jet.pt()));

   properties.push_back( MOD::Property("mul_pre_SD", (int) hardest_jet_constituents.size()) );
   properties.push_back( MOD::Property("mul_post_SD", (int) soft_drop(hardest_jet).constituents().size() ) );
    
   properties.push_back( MOD::Property("mul_filtered_SD", (int) filtered_hardest_jet_constituents.size()) );


   properties.push_back( MOD::Property("mass_pre_SD", hardest_jet.m()) );
   properties.push_back( MOD::Property("mass_post_SD", soft_drop(hardest_jet).m()) );


   properties.push_back( MOD::Property("pT_D_pre_SD", pT_D(hardest_jet)) );
   properties.push_back( MOD::Property("pT_D_post_SD", pT_D(soft_drop(hardest_jet)) ));

   properties.push_back( MOD::Property("LHA_pre_SD", angularity_lambda(hardest_jet, jet_radius, 1, 0.5)) );
   properties.push_back( MOD::Property("LHA_post_SD", angularity_lambda(soft_drop(hardest_jet), jet_radius, 1, 0.5)) );

   properties.push_back( MOD::Property("width_pre_SD", angularity_lambda(hardest_jet, jet_radius, 1, 1)) );
   properties.push_back( MOD::Property("width_post_SD", angularity_lambda(soft_drop(hardest_jet), jet_radius, 1, 1)) );

   properties.push_back( MOD::Property("thrust_pre_SD", angularity_lambda(hardest_jet, jet_radius, 1, 2)) );
   properties.push_back( MOD::Property("thrust_post_SD", angularity_lambda(soft_drop(hardest_jet), jet_radius, 1, 2)) );


    
    //----------------------------------------------------------
    // illustrate how this SubjetCounting contrib works
    
    

   // ================================================================ Track Based Analysis ================================================================

/*


   // Get all charged particles with 0.5 GeV particles removed.
   vector<fastjet::PseudoJet> track_constituents = MOD::filter_charged(pT_0_5_GeV_selector(event_being_read.hardest_jet().constituents()));

   // Cluster them using Cambridge/Alachen with infinite radius. This makes sure that we get the same jets as "regular" ak5 jets except now with just charged particles.
   ClusterSequence cs_track(track_constituents, jet_def_cambridge);

   if (cs_track.inclusive_jets().size() > 0 ) {

      PseudoJet hardest_track_jet = cs_track.inclusive_jets()[0];

      for (unsigned i = 0; i < zg_cuts.size(); i++) {
         string label = zg_cuts[i].first;
         double zg_cut = zg_cuts[i].second;

         SoftDrop soft_drop_track(beta, zg_cut);
         PseudoJet soft_drop_jet_track = soft_drop_track(hardest_track_jet);

         double zg = soft_drop_jet_track.structure_of<SoftDrop>().symmetry();
         double rg = soft_drop_jet_track.structure_of<SoftDrop>().delta_R() / 0.5;
       

         properties.push_back(MOD::Property("track_zg_" + label, zg));
         properties.push_back(MOD::Property("track_mu_" + label, soft_drop_jet_track.structure_of<SoftDrop>().mu()));

         properties.push_back(MOD::Property("track_rg_" + label, rg));
         properties.push_back(MOD::Property("track_e1_" + label, rg * zg));
         properties.push_back(MOD::Property("track_e2_" + label, pow(rg, 2) * zg ));
         properties.push_back(MOD::Property("track_e05_" + label, sqrt(rg) * zg ));
      }
 

      properties.push_back( MOD::Property("track_mul_pre_SD", (int) hardest_track_jet.constituents().size()) );
      properties.push_back( MOD::Property("track_mul_post_SD", (int) soft_drop(hardest_track_jet).constituents().size()) );

      properties.push_back( MOD::Property("track_mass_pre_SD", hardest_track_jet.m()) );
      properties.push_back( MOD::Property("track_mass_post_SD", soft_drop(hardest_track_jet).m()) );


      properties.push_back( MOD::Property("track_pT_D_pre_SD", pT_D(hardest_track_jet)) );
      properties.push_back( MOD::Property("track_pT_D_post_SD", pT_D(soft_drop(hardest_track_jet))) );

      properties.push_back( MOD::Property("track_LHA_pre_SD", angularity_lambda(hardest_track_jet, jet_radius, 1, 0.5)) );
      properties.push_back( MOD::Property("track_LHA_post_SD", angularity_lambda(soft_drop(hardest_track_jet), jet_radius, 1, 0.5)) );

      properties.push_back( MOD::Property("track_width_pre_SD", angularity_lambda(hardest_track_jet, jet_radius, 1, 1)) );
      properties.push_back( MOD::Property("track_width_post_SD", angularity_lambda(soft_drop(hardest_track_jet), jet_radius, 1, 1)) );

      properties.push_back( MOD::Property("track_thrust_pre_SD", angularity_lambda(hardest_track_jet, jet_radius, 1, 2)) );
      properties.push_back( MOD::Property("track_thrust_post_SD", angularity_lambda(soft_drop(hardest_track_jet), jet_radius, 1, 2)) );

   }
   else {

      for (unsigned i = 0; i < zg_cuts.size(); i++) {
         string label = zg_cuts[i].first;
         properties.push_back(MOD::Property("track_zg_" + label, -1.0));
         properties.push_back(MOD::Property("track_mu_" + label, -1.0));

         properties.push_back(MOD::Property("track_rg_" + label, -1.));
         properties.push_back(MOD::Property("track_e1_" + label, -1.));
         properties.push_back(MOD::Property("track_e2_" + label, -1.));
         properties.push_back(MOD::Property("track_e05_" + label, -1.));
      }

      properties.push_back( MOD::Property("track_mul_pre_SD", -1. ));
      properties.push_back( MOD::Property("track_mul_post_SD", -1. ));

      properties.push_back( MOD::Property("track_mass_pre_SD", -1. ));
      properties.push_back( MOD::Property("track_mass_post_SD", -1. ));


      properties.push_back( MOD::Property("track_pT_D_pre_SD", -1.) );
      properties.push_back( MOD::Property("track_pT_D_post_SD", -1.) );

      properties.push_back( MOD::Property("track_LHA_pre_SD", -1.) );
      properties.push_back( MOD::Property("track_LHA_post_SD", -1.) );

      properties.push_back( MOD::Property("track_width_pre_SD", -1.) );
      properties.push_back( MOD::Property("track_width_post_SD", -1.) );

      properties.push_back( MOD::Property("track_thrust_pre_SD", -1.) );
      properties.push_back( MOD::Property("track_thrust_post_SD", -1.) );


   }
    */
   
 
   // Now that we've calculated all observables, write them out.

   string name;
   
   int padding = 35;

   if (event_serial_number == 1) {
      for (unsigned p = 0; p < properties.size(); p++) {
         if (p > 0)
            output_file << setw(padding);
         
         output_file << properties[p].name();
      }

      output_file << endl;
   }

   for (unsigned q = 0; q < properties.size(); q++) {
      if (q > 0)
         output_file << setw(padding);
      output_file << fixed << setprecision(8) << properties[q];
   }

   output_file << endl;
   
   
}





double angularity_lambda(PseudoJet jet, double jet_radius, float k, float beta) {
      

   JetDefinition jet_def_cambridge(cambridge_algorithm, fastjet::JetDefinition::max_allowable_R, WTA_pt_scheme);
   ClusterSequence cs = ClusterSequence(jet.constituents(), jet_def_cambridge);
   

   if (cs.inclusive_jets().size() == 0) {
      return 0;
   }

   PseudoJet recombined_jet = cs.inclusive_jets()[0];

   double lambda = 0.0;

   double R = jet_radius;   // Jet Radius.

   double total_pT = 0.0;
   for (unsigned j = 0; j < recombined_jet.constituents().size(); j++) {
      total_pT += recombined_jet.constituents()[j].pt();
   }

   for (unsigned i = 0; i < recombined_jet.constituents().size(); i++) {
      
      PseudoJet constituent = recombined_jet.constituents()[i];

      double z_i = constituent.pt() / total_pT;
      
      double delta_R = constituent.delta_R(jet);

      double theta_i = delta_R / R;

      lambda += pow(z_i, k) * pow(theta_i, beta);

   }

   return lambda;

}







double pT_D(PseudoJet jet) {
      
   if (jet.constituents().size() == 0) {
      return 0;
   }

   double total_pT = 0.0;
   for (unsigned j = 0; j < jet.constituents().size(); j++) {
      total_pT += jet.constituents()[j].pt();
   }

   double numerator = 0.0;
   for (unsigned i = 0; i < jet.constituents().size(); i++) {
      numerator += jet.constituents()[i].pt() * jet.constituents()[i].pt();
   }

   return pow(numerator, 0.5) / total_pT;

}

unordered_map<int, pair <int, double>> get_lumi_info_lumi_block(string & input_file, int number_of_events_to_process) {
    
    unordered_map<int, pair <int, double>> lumi_block_lumi_info;

    ifstream data_file("/Users/prekshanaik/Documents/MEng_Project/11MODAnalyzer-master/2011lumibyls.csv");
    vector <vector <string> > data;
    
    while (data_file)
    {
        string s;
        if (!getline( data_file, s )) break;
        
        istringstream ss( s );
        vector <string> record;
        
        while (ss)
        {
            string s;
            if (!getline( ss, s, ',' )) break;
            record.push_back( s );
        }
        
        data.push_back( record );
    
    }
    
    for (int i=0; i< data.size(); i++) {
        
        if (data[i][0].find('#') == string::npos) {
            
            replace(data[i][0].begin(), data[i][0].end(), ':', ' ');
            vector<int> array;
            stringstream ss(data[i][0]);
            int temp;
            while (ss >> temp)
                array.push_back(temp);
                
            int run_number = array[0];
            
            replace(data[i][1].begin(), data[i][1].end(), ':', ' ');
            vector<int> array2;
            stringstream ss2(data[i][1]);
            int temp2;
            while (ss2 >> temp2)
                array2.push_back(temp2);

            int lumi_block = array2[0];

            double recorded_lumi = atof(data[i][6].c_str());
        
            //cout << "this is the lumi_data" << lumi_block << " " << run_number << " " << recorded_lumi << endl;
            
            int hash_value = run_number*1000+lumi_block;
            
            auto it = lumi_block_lumi_info.find(hash_value);
            if(it != lumi_block_lumi_info.end()){
                pair <int, double> current_pair;
                current_pair = it->second;
                current_pair.first += 1;
                current_pair.second += recorded_lumi;
                it->second = current_pair;
            }
            else{
                pair <int, double> new_pair;
                new_pair = make_pair(1, recorded_lumi);
                lumi_block_lumi_info[hash_value] = new_pair;
            }
            
            //cout << hash_value << endl;

        
        }
    
    }
    
    int hash_test = 160578000;
    
    cout << lumi_block_lumi_info[hash_test].first << endl;
    
    return lumi_block_lumi_info;
    
    

}



