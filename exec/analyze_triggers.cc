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

void analyze_event(MOD::Event & event_being_read, ofstream & output_file, int & event_serial_number);
unordered_map<int, double> get_lumi_info_by_block(string lumi_data_file_name);


int main(int argc, char * argv[]) {

   auto start = std::chrono::steady_clock::now();

   int number_of_events_to_process;

   if (argc <= 2) {
        std::cerr << "ERROR: You need to supply three arguments- first, path to the input data; second, path to the output file; third, number of events to process. The path has to be either absolute or relative to the bin directory:" << std::endl << std::endl << "./analysis (input_file.dat) (output_file.dat) [optional Nev]" << std::endl;
        return 1;
   }
   else if (argc == 3) {
      // Third argument is missing, process everything.
      number_of_events_to_process = std::numeric_limits<int>::max();
   }
   else {
      // Third argument gives the number of events to process.
      number_of_events_to_process = stoi(argv[3]);
   }

   ifstream data_file(argv[1]);
   ofstream output_file(argv[2], ios::out | ios::app);

   
   cout << endl << endl << "Starting analysis with the following given arguments: " << endl;
   cout << "Input file: " << argv[1] << endl;
   cout << "Output file: " << argv[2] << endl;
   cout << "Number of events: ";
   if(argc == 3)
      cout << "ALL" << endl << endl;
   else
      cout << number_of_events_to_process << endl << endl;

   //vector<double> cone_radii = {0.3, 0.5, 0.7};
   //svector<double> pt_cuts = {50.0, 80.0, 110.0};

   MOD::Event event_being_read;
   string lumi2011_file_to_process = "2011lumibyls.csv";

    
   unordered_map<int, double> lumi_block_lumi_info = get_lumi_info_by_block(lumi2011_file_to_process);
    
    

   int event_serial_number = 1;
   while( event_being_read.read_event(data_file) && ( event_serial_number <= number_of_events_to_process ) ) {
      
      if( (event_serial_number % 5000) == 0 )
         cout << "Processing event number " << event_serial_number << endl;
       
       // Compute 'hash_value' for lumi_block. Since number of lumi_block less than <1000, this is okay as unique id
       int lumi_block = event_being_read.condition().lumi_block();
       int run_number = event_being_read.condition().run_number();
       int hash_value = run_number*1000+lumi_block;


       if(lumi_block_lumi_info.find(hash_value) != lumi_block_lumi_info.end()){
           analyze_event(event_being_read, output_file, event_serial_number);
       }
      
      event_being_read = MOD::Event();
      event_serial_number++;
   }

   auto finish = std::chrono::steady_clock::now();
   double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start).count();
   cout << "Finished processing " << (event_serial_number - 1) << " events in " << elapsed_seconds << " seconds!" << endl;

   return 0;
}


void analyze_event(MOD::Event & event_being_read, ofstream & output_file, int & event_serial_number) {

   fastjet::PseudoJet trigger_jet = event_being_read.trigger_jet();
   fastjet::PseudoJet const closest_jet_to_trigger_jet = event_being_read.closest_fastjet_jet_to_trigger_jet();

   vector<MOD::Property> properties;
   
   if ( (event_being_read.cms_jets().size() == 0) or (event_being_read.jets().size() == 0) or ( ! trigger_jet.has_user_info())) {
      return;
   }

   vector<MOD::Trigger> triggers = event_being_read.triggers();
  

   
   try {
      for (unsigned i = 0; i < triggers.size(); i++) {
         
         if (triggers[i].fired()) {

            fastjet::PseudoJet trigger_jet = event_being_read.trigger_jet();

            properties.push_back(MOD::Property("# Entry", "  Entry"));

            properties.push_back(MOD::Property("event_number", event_being_read.event_number()));
            properties.push_back(MOD::Property("run_number", event_being_read.run_number()));

            properties.push_back(MOD::Property("trig_jet_matched", (int) event_being_read.is_trigger_jet_matched())); 
            properties.push_back(MOD::Property("jet_quality", trigger_jet.user_info<MOD::InfoCalibratedJet>().jet_quality())); 
   
            properties.push_back(MOD::Property("hardest_pT", trigger_jet.pt()));  
            properties.push_back(MOD::Property("corr_hardest_pT", trigger_jet.pt() * trigger_jet.user_info<MOD::InfoCalibratedJet>().JEC()));  
            properties.push_back(MOD::Property("hardest_eta", trigger_jet.eta()));  
            properties.push_back(MOD::Property("prescale", triggers[i].prescale()));
            properties.push_back(MOD::Property("trigger_name", triggers[i].name()));
       
            string name;
   
            int padding = 36;

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

            properties.clear();

         }
      }
   }
   catch (exception& e) {

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

   
