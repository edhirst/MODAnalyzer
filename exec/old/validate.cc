#include <iostream>
#include <unordered_map>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <iterator>
#include <iomanip>

#include "fastjet/ClusterSequence.hh"
#include "../interface/event.h"
#include "../interface/fractional_jet_multiplicity.h"

using namespace std;
using namespace fastjet;

void validate_events(MOD::Event & event_being_read, ofstream & output_file);
bool jets_match(MOD::Event & event_being_read);
bool pseudojets_compare(PseudoJet a, PseudoJet b);
bool write_jets(MOD::Event & event_being_read, ofstream cms_output_file, ofstream fastjet_output_file);


int main(int argc, char * argv[]) {
   
   int number_of_events_to_process;

   if (argc <= 2) {
        std::cerr << "ERROR: You need to supply three arguments- first, path to the input data; second, path to the output file; third, number of events to process. The path has to be either absolute or relative to the bin directory:" << std::endl << std::endl << "./duplication (input_file.dat) (output_file.dat) [optional Nev]" << std::endl;
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
   ofstream output_file(argv[2], ios::out);

   ofstream cms_output_file("data/cms.dat", ios::out);
   ofstream fastjet_output_file("data/fastjet.dat", ios::out);

   cout << endl << endl << "Starting validation with the following given arguments: " << endl;
   cout << "Input file: " << argv[1] << endl;
   cout << "Output file: " << argv[2] << endl;
   cout << "Number of events: ";
   if(argc == 3)
      cout << "ALL" << endl << endl;
   else
      cout << number_of_events_to_process << endl << endl;
   
   int events_with_mismatched_jets = 0;

   MOD::Event event_being_read;

   int event_serial_number = 1;
   while( event_being_read.read_event(data_file) && ( event_serial_number <= number_of_events_to_process ) ) {

      if( (event_serial_number % 100) == 0 )
         cout << "Validating event number " << event_serial_number << endl;

      output_file << event_being_read;

      if( ! jets_match(event_being_read)) {
         events_with_mismatched_jets++;
         cout << "AK5 Jets don't match for event: " << event_being_read.event_number() << endl << endl;
      }


      // write_jets(event_being_read, 0.5, cms_output_file, fastjet_output_file);

      
      event_being_read = MOD::Event();
      event_serial_number++;
   }
   
   cout << "Validation complete!" << endl << endl;
   cout << events_with_mismatched_jets << " events have mismatched AK5 jets!" << endl;
   

}

bool jets_match(MOD::Event & event_being_read) {
   
   double pt_cut = 3.00;
   double cone_radius = 0.5;
   
   vector<PseudoJet> cms_jets = event_being_read.CMS_pseudojets();
   vector<PseudoJet> pfcandidates = event_being_read.pseudojets();

   // Cluster the pfcandidates using Fastjet.
   JetDefinition jet_def(antikt_algorithm, cone_radius);
   ClusterSequence cs(pfcandidates, jet_def);
   vector<PseudoJet> fastjet_jets = cs.inclusive_jets(pt_cut);

   // Compare the number of jets first.
   if (cms_jets.size() != fastjet_jets.size()) {
      cout << "Different jet size for CMS vs. FastJet; " << cms_jets.size() << " vs " << fastjet_jets.size() << "." << endl;
      return false;
   }

   // Next, compare if fastjet_jets with cms_jets or not upto 10e-4 precision.
   
   unsigned max_number_of_jets =  max(fastjet_jets.size(), cms_jets.size());

  
   // First sort both vectors by px so that we can compare them one by one.
   sort(fastjet_jets.begin(), fastjet_jets.end(), pseudojets_compare);
   sort(cms_jets.begin(), cms_jets.end(), pseudojets_compare);

   double tolerance = pow(10, -3);
   for (unsigned i = 0; i < max_number_of_jets; i++) {
      
      
      if ( ( abs(fastjet_jets[i].px() - cms_jets[i].px()) > tolerance ) || ( abs(fastjet_jets[i].py() - cms_jets[i].py()) > tolerance ) || ( abs(fastjet_jets[i].pz() - cms_jets[i].pz()) > tolerance ) || ( abs(fastjet_jets[i].E() - cms_jets[i].E()) > tolerance ) ) {
         
         cout << fixed << setprecision(5) << "FastJet: " <<  fastjet_jets[i].px() << "   " << fastjet_jets[i].py() << "   " << fastjet_jets[i].pz() << "   " << fastjet_jets[i].E() << endl;
         cout << fixed << setprecision(5) << "CMS:     " << cms_jets[i].px() << "   " << cms_jets[i].py() << "   " << cms_jets[i].pz() << "   " << cms_jets[i].E() << endl << endl;
         
         return false;

      }
   }

   return true;  
}

bool pseudojets_compare(PseudoJet a, PseudoJet b) {
   if (a.pt() > b.pt())
      return true;
   return false;
}

bool write_jets(MOD::Event & event_being_read, ofstream cms_output_file, ofstream fastjet_output_file) {

   double pt_cut = 3.00;
   double cone_radius = 0.5;
    
   vector<PseudoJet> cms_jets = event_being_read.CMS_pseudojets();
   
   // First sort the jets by px so that we can compare them one by one.
   sort(cms_jets.begin(), cms_jets.end(), pseudojets_compare);

   // Write out CMS jets.
   cms_output_file << "            px            py            pz        energy" << endl;
   for (unsigned i = 0; i < cms_jets.size(); i++) {
      cms_output_file 
                     << setw(14) << fixed << setprecision(8) << cms_jets[i].px()
                     << setw(14) << fixed << setprecision(8) << cms_jets[i].py()
                     << setw(14) << fixed << setprecision(8) << cms_jets[i].pz()
                     << setw(14) << fixed << setprecision(8) << cms_jets[i].E()
                     << endl;
   }

   vector<PseudoJet> pfcandidates = event_being_read.pseudojets();

   // Cluster the pfcandidates using Fastjet.
   JetDefinition jet_def(antikt_algorithm, cone_radius);
   ClusterSequence cs(pfcandidates, jet_def);
   vector<PseudoJet> fastjet_jets = cs.inclusive_jets(pt_cut);

   // First sort the jets by px so that we can compare them one by one.
   sort(fastjet_jets.begin(), fastjet_jets.end(), pseudojets_compare);

   // Write out FastJet jets.
   fastjet_output_file << "            px            py            pz        energy       pt" << endl;
   for (unsigned i = 0; i < fastjet_jets.size(); i++) {
      fastjet_output_file 
                     << setw(14) << fixed << setprecision(8) << fastjet_jets[i].px()
                     << setw(14) << fixed << setprecision(8) << fastjet_jets[i].py()
                     << setw(14) << fixed << setprecision(8) << fastjet_jets[i].pz()
                     << setw(14) << fixed << setprecision(8) << fastjet_jets[i].E()
                     << setw(14) << fixed << setprecision(8) << fastjet_jets[i].pt()
                     << endl;
   }

   return true;  
}