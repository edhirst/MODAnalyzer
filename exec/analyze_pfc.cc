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


#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/SoftDrop.hh"

#include "../interface/Event.h"
#include "../interface/Property.h"

using namespace std;
using namespace fastjet;
using namespace contrib;


void analyze_pfc(MOD::Event & event_being_read, ofstream & output_file, int & event_serial_number,string & data_set_to_process);
unordered_map<int, double> get_lumi_info_by_block(string lumi_data_file_name);
string get_just_file_name(string & input_file_name);
string replace_string(string subject, const string& search, const string& replace);

//void get_all_files_to_process(std::vector<string> & all_files, vector<boost::filesystem::path> input_paths);

double angularity_lambda(PseudoJet jet, float k, float beta);

int main(int argc, char * argv[]) {
    
    auto start = std::chrono::steady_clock::now();
    
    // Default arguments - process all events, lumi file located in same directory
    int number_of_events_to_process = numeric_limits<int>::max();
    string lumi2011_file_to_process = "2011lumibyls.csv";
    
    
    if (argc <= 3) {
        std::cerr << "ERROR: You need to supply at least three arguments- first, path to the input data; second, path to the output file; third, whether this is 2011 or simulated data; The paths should be either absolute or relative to the bin directory:" << std::endl << std::endl << "./analysis (input_file.dat) (output_file.dat) (sim or 2011)[optional Nev]" << std::endl;
        return 1;
    }
    
    string data_set_to_process(argv[3]);
    if ((data_set_to_process != "sim") and (data_set_to_process != "2011")){
        std::cerr << "You have given an invalid flag: " << argv[3] << ". Please specify either '2011' or 'sim' " << std::endl;
        return 1;
    }
    
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-2011_lumi_file") == 0){
            string new_lumi_file_location(argv[i+1]);
            lumi2011_file_to_process = new_lumi_file_location;
        }
        if (strcmp(argv[i], "-number_events") == 0){
            number_of_events_to_process = stoi(argv[i+1]);
        }
    }
    
    cout << endl << endl << "Starting analysis with the following given arguments: " << endl;
    cout << "Input file: " << argv[1] << endl;
    cout << "Output file: " << argv[2] << endl;
    cout << "Number of events: ";
    if(number_of_events_to_process == numeric_limits<int>::max()){    cout << "ALL " << endl;} else{cout << number_of_events_to_process << endl;}
    if (data_set_to_process == "sim") {cout << "Processing Simulated Data" << endl << endl;} else {cout << "Processing 2011 Data with file: " << lumi2011_file_to_process << endl << endl;}
    
    
    int event_serial_number_2011 = 1;
    int event_serial_number_analyzed_2011 = 1;
    if (data_set_to_process == "2011") {
        
        
        ifstream data_file(argv[1]);
        string output_file_name(argv[2]);
        ofstream output_file(output_file_name, ios::out);
        
        unordered_map<int, double> lumi_block_lumi_info = get_lumi_info_by_block(lumi2011_file_to_process);
        
        string possible_triggers[9] = { "HLT_Jet30", "HLT_Jet60", "HLT_Jet80", "HLT_Jet110", "HLT_Jet150", "HLT_Jet190", "HLT_Jet240", "HLT_Jet300", "HLT_Jet370" };
        unordered_map<string, double> recorded_luminosities;
        unordered_map<string, int> checked_lumi_blocks;
        
        MOD::Event event_being_read;
        
        while( event_being_read.read_event(data_file) && ( event_serial_number_2011 <= number_of_events_to_process ) ) {
            
            if( (event_serial_number_2011 % 500) == 0 ){
                cout << "Processing event number " << event_serial_number_2011 << endl;
            }
            
            // cout << "Processing event number " << event_serial_number_2011 << endl;
            
            // Compute 'hash_value' for lumi_block. Since number of lumi_block less than <1000, this is okay as unique id
            int lumi_block = event_being_read.condition().lumi_block();
            int run_number = event_being_read.condition().run_number();
            int hash_value = run_number*1000+lumi_block;
            
            // If the event being read has the assigned trigger fired and event is good, analyze and write it out
            if (event_being_read.assigned_trigger_fired() and (lumi_block_lumi_info.find(hash_value) != lumi_block_lumi_info.end())) {
                analyze_pfc(event_being_read, output_file, event_serial_number_analyzed_2011, data_set_to_process);
                event_serial_number_analyzed_2011++;
            }
            
            // Go through all possible triggers (regardless of whether assigned trigger fired), check if trigger exists and event is good, and keep cumulative total of effective luminosities
            // Checks to make sure for each trigger, each lumi_block is only counted once
            for (int i = 0; i < 10; i++ ){
                string trigger_name = possible_triggers[i];
                if (event_being_read.trigger_exists(trigger_name) and (lumi_block_lumi_info.find(hash_value) != lumi_block_lumi_info.end()) ){
                    MOD::Trigger trigger_current;
                    trigger_current = event_being_read.trigger_by_short_name(trigger_name);
                    int trigger_current_prescale;
                    trigger_current_prescale = trigger_current.prescale();
                    if (checked_lumi_blocks.find(to_string(hash_value)+trigger_name) == checked_lumi_blocks.end()){
                        recorded_luminosities[trigger_name] += lumi_block_lumi_info[hash_value]/trigger_current_prescale*1.0;
                        checked_lumi_blocks[to_string(hash_value)+trigger_name] = 1;
                    }
                    
                    
                }
            }
            
            event_serial_number_2011++;
            event_being_read = MOD::Event();
            
        }
        
        ofstream effective_luminosity_output("effective_luminosity_by_trigger.csv", ios::out | ios::app);
        
        unordered_map<string, double>::iterator it = recorded_luminosities.begin();
        
        // For each trigger, write out the mod_file it came from and the trigger_name. We keep track of the mod_file to make sure values are not double counted and
        // also as a check to make sure analysis completed.
        while (it != recorded_luminosities.end())
        {
            string trigger_name = it->first;
            
            double luminosity = it->second;
            
            string input_file_name(argv[1]);
            effective_luminosity_output << get_just_file_name(input_file_name) << "," << trigger_name << "," << luminosity;
            effective_luminosity_output << endl;
            
            it++;
        }
        effective_luminosity_output << endl;
        
        cout << "finished writing to luminosity by trigger" << endl;
        
    }
    int event_serial_number_pfc = 1;
    int event_serial_number_gen = 1;
    
    if (data_set_to_process == "sim"){
        
        string output_file_name(argv[2]);
        
        
        // First, analyze and write out the pfc_data
        MOD::Event event_being_read;
        string output_file_pfc_name = replace_string(output_file_name, ".dat", "_sim_pfc.dat");
        ofstream output_file_pfc(output_file_pfc_name, ios::out);
        ifstream data_file_pfc(argv[1]);
        
        while( event_being_read.read_event(data_file_pfc, "sim_pfc") && ( event_serial_number_pfc <= number_of_events_to_process ) ) {
            
            if( (event_serial_number_pfc % 500) == 0 ){
                cout << "Processing event number " << event_serial_number_pfc << endl;
            }
            
            analyze_pfc(event_being_read, output_file_pfc, event_serial_number_pfc, data_set_to_process);
            event_being_read = MOD::Event();
            event_serial_number_pfc++;
        }
        
        cout << event_serial_number_pfc << " " << event_serial_number_gen << " " << number_of_events_to_process << endl;
        
        // Next, analyze and write out the gen_data
        string output_file_gen_name = replace_string(output_file_name, ".dat", "_sim_gen.dat");
        ofstream output_file_gen(output_file_gen_name, ios::out);
        ifstream data_file_gen(argv[1]);
        
        while( event_being_read.read_event(data_file_gen, "sim_gen") && ( event_serial_number_gen <= number_of_events_to_process ) ) {
            
            if( (event_serial_number_gen % 500) == 0 ){
                cout << "Processing event number " << event_serial_number_gen << endl;
            }
            
            analyze_pfc(event_being_read, output_file_gen, event_serial_number_gen, data_set_to_process);
            event_being_read = MOD::Event();
            event_serial_number_gen++;
        }
        
        ofstream output_event_counts("event_count_by_pythia_and_mod.csv", ios::out | ios::app);
        output_event_counts << get_just_file_name(output_file_pfc_name) << "," << event_serial_number_pfc -1  << endl;
        output_event_counts << get_just_file_name(output_file_gen_name) << "," << event_serial_number_gen -1  << endl;
        
        
    }
    
    auto finish = std::chrono::steady_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start).count();
    if (data_set_to_process == "2011"){cout << "Finished processing " << (event_serial_number_2011 - 1) << " events in " << elapsed_seconds << " seconds!" << endl;}
    if (data_set_to_process == "sim"){cout << "Finished processing " << (event_serial_number_pfc - 1) << " PFC events in " << elapsed_seconds << " seconds!" << endl;
        cout << "Finished processing " << (event_serial_number_gen - 1) << " GEN events in " << elapsed_seconds << " seconds!" << endl;}
    
    return 0;
    
}


void analyze_pfc(MOD::Event & event_being_read, ofstream & output_file, int & event_serial_number, string & data_set_to_process) {


   JetDefinition jet_def_cambridge(cambridge_algorithm, fastjet::JetDefinition::max_allowable_R);
   PseudoJet hardest_jet = event_being_read.hardest_jet();


   if ( ! (hardest_jet.E() > 0.0)) {
      return;
   }


   double jec;
   if ((event_being_read.data_source() == 0) || (event_being_read.data_source() == 3)) {	// EXPERIMENT or PRISTINE
   	jec = event_being_read.get_hardest_jet_jec();	
   }
   else {
   	jec = 1.0;	
   }
   
   
 
   vector<PseudoJet> hardest_jet_pfcs = hardest_jet.constituents();
      
   for (unsigned i = 0; i < hardest_jet_pfcs.size(); i++) {

   		vector<MOD::Property> properties;

   		properties.push_back(MOD::Property("# Entry", "  Entry"));

        properties.push_back(MOD::Property("event_number", event_being_read.event_number()));
        properties.push_back(MOD::Property("prescale", event_being_read.weight()));
	   	properties.push_back(MOD::Property("hardest_pT", jec * hardest_jet.pt()));
	   	properties.push_back(MOD::Property("jet_eta", hardest_jet.eta()));
	   	properties.push_back(MOD::Property("pfc_pT", hardest_jet_pfcs[i].pt()));
        // cout << "hardest jet pfc: " << hardest_jet_pfcs[i];
	   	properties.push_back(MOD::Property("pfc_eta", hardest_jet_pfcs[i].eta()));
       
	   	properties.push_back(MOD::Property("pfc_pdgId", hardest_jet_pfcs[i].user_info<MOD::InfoPFC>().pdgId()));
       if (data_set_to_process == "2011"){   properties.push_back(MOD::Property("trigger_fired", event_being_read.assigned_trigger_name()));}
       
       if (data_set_to_process == "sim"){   properties.push_back(MOD::Property("trigger_fired", "no_trigger"));}

   		// Now that we've calculated all observables, write them out.

	   string name;
	   
	   int padding = 40;

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
	      output_file << properties[q];
	   }

	   output_file << endl;


   }

   

   
   
   
}


unordered_map<int, double> get_lumi_info_by_block(string lumi_data_file_name){
    
    unordered_map<int, double> lumi_info_by_block;
    
    ifstream lumi_data_file(lumi_data_file_name);
    vector <vector <string> > data;
    
    while (lumi_data_file)
    {
        string s;
        if (!getline( lumi_data_file, s )) break;
        
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
            
            
            int hash_value = run_number*1000+lumi_block;
            
            auto it = lumi_info_by_block.find(hash_value);
            if(it != lumi_info_by_block.end()){
                it->second += recorded_lumi;
            }
            else{
                lumi_info_by_block[hash_value] = recorded_lumi;
            }
        }
    }
    
    
    return lumi_info_by_block;
}


string get_just_file_name(string & input_file_name){
    
    char c = '/';
    string buff{""};
    vector<string> v;
    
    for(auto n:input_file_name)
    {
        if(n != c) buff+=n; else
            if(n == c && buff != "") { v.push_back(buff); buff = ""; }
    }
    if(buff != "") v.push_back(buff);
    
    return v.back();
    
}




double angularity_lambda(PseudoJet jet, float k, float beta) {
   
   double lambda = 0.0;

   double R = 0.5;   // Jet Radius.

   double total_pT = 0.0;
   for (unsigned j = 0; j < jet.constituents().size(); j++) {
      total_pT += jet.constituents()[j].pt();
   }

   for (unsigned i = 0; i < jet.constituents().size(); i++) {
      
      PseudoJet constituent = jet.constituents()[i];

      double z_i = constituent.pt() / total_pT;
      
      double delta_R = constituent.delta_R(jet);

      double theta_i = delta_R / R;

      lambda += pow(z_i, k) * pow(theta_i, beta);

   }

   return lambda;

}

string replace_string(string subject, const string& search, const string& replace) {
    size_t pos = 0;
    while ((pos = subject.find(search, pos)) != std::string::npos) {
        subject.replace(pos, search.length(), replace);
        pos += replace.length();
    }
    return subject;
}

