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
#include <utility>
#include <chrono>
#include <algorithm>

#include <boost/unordered_map.hpp>
#include <boost/filesystem.hpp>

#include "../interface/Event.h"
#include "../interface/Property.h"

using namespace std;
using namespace boost::filesystem;

void analyze_event(MOD::Event & event_being_read, string output_path, unordered_map<string, string> & registry_info, ofstream & completed_events_file_output, boost::unordered_map<string, int> & completed_events, int & number_of_events_written);
string find_correct_file_for_event(int run_number, int event_number, unordered_map<string, string> & registry_info);
void load_registry_information(string registry_filename, unordered_map<string, string> & registry_info);
void get_all_files_to_process(std::vector<string> & all_files, boost::filesystem::path input_path);

int main(int argc, char * argv[]) {

   auto start = std::chrono::steady_clock::now();


   if (argc <= 3) {
        std::cerr << "ERROR: You need to supply five arguments- first, path to the input data; second, path to the output file; third, number of events to process. The path has to be either absolute or relative to the bin directory:" << std::endl << std::endl << "./analysis (input_file.dat) (output_file.dat) [optional Nev]" << std::endl;
        return 1;
   }


   string input_path(argv[1]);
   string output_path(argv[2]);
   string registry_file(argv[3]);
   string completed_log_filename(argv[4]);

   
   cout << endl << endl << "Starting process with the following given arguments: " << endl;
   cout << "Input Path: " << argv[1] << endl;
   cout << "Output Path   : " << argv[2] << endl;
   cout << "Registry File : " << argv[3] << endl;
   cout << "Completed Log : " << argv[4] << endl << endl << endl;

   unordered_map<string, string> registry_info;


   load_registry_information(registry_file, registry_info);

   cout << "Finished loading registry to memory." << endl;


   boost::unordered_map<string, int> completed_events;

   // Load completed events to the vector "completed_events_"
   ifstream completed_file(completed_log_filename.c_str());
   
   string line;
   int line_number = 1;
   while(getline(completed_file, line)) {
      
      if (line_number % 100000 == 0)
         cout << "On line number " << line_number << endl;

      istringstream iss(line);
      int event_number, run_number;
      iss >> run_number >> event_number;
      
      completed_events.insert(make_pair(to_string(run_number) + "_" + to_string(event_number), 1));

      line_number++;
   }

   cout << "Finished loading \"completed\" events to memory." << endl;

   ofstream completed_events_file_output;
   completed_events_file_output.open(completed_log_filename.c_str(), ios::out | ios::app );


   ofstream log;
   log.open("./log_0000.dat", ios::out | ios::app);

   vector<string> all_filenames;
   get_all_files_to_process(all_filenames, input_path);

   int file_counter = 0;
   for (unsigned i=0; i < all_filenames.size(); i++) {
      ifstream input_mod_file(all_filenames[i]);

      file_counter++;

      // if ((file_counter % 100) == 0)
      //    cout << "Processing file number " << file_counter << " / " << all_filenames.size() << endl;

      cout << "Processing file " << (i + 1) << " / " << all_filenames.size() << ": " << all_filenames[i] << endl;
      
      MOD::Event event_being_read;

      int event_serial_number = 1;
      int number_of_events_written = 0;
      while ( event_being_read.read_event(input_mod_file) ) {
         
         // if( (event_serial_number % 5000) == 0 )
         //    cout << "Processing event number " << event_serial_number << endl;

         
         analyze_event(event_being_read, output_path, registry_info, completed_events_file_output, completed_events, number_of_events_written);
         
         event_being_read = MOD::Event();
         event_serial_number++;
      }

      log  << "Written " << number_of_events_written << " / " << (event_serial_number - 1) << " events for " << all_filenames[i] << endl;
   }


   auto finish = std::chrono::steady_clock::now();
   double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start).count();
   cout << "Finished processing " << all_filenames.size() << " files in " << elapsed_seconds << " seconds!" << endl;

   return 0;
}


void analyze_event(MOD::Event & event_being_read, string output_path, unordered_map<string, string> & registry_info, ofstream & completed_events_file_output, boost::unordered_map<string, int> & completed_events, int & number_of_events_written) {

   ofstream log_file(output_path + "/error.log", ios::out | ios::app);

   if (completed_events.find(to_string(event_being_read.run_number()) + "_" + to_string(event_being_read.event_number())) == completed_events.end()) {
      
      // Event not processed already.

      string correct_filename = find_correct_file_for_event(event_being_read.run_number(), event_being_read.event_number(), registry_info);

      // cout << "Correct filename is " << output_path << "/" << correct_filename << ".mod" << endl;

      if (correct_filename == "") {
         log_file << "Could not find the correct filename for the following event: " << event_being_read.event_number() << " " << event_being_read.run_number() << endl;
         return;
      }
      else {
         ofstream output_file(output_path + "/" + correct_filename + ".mod", ios::out | ios::app);

         // See if the event has already been added to the list. 
         
         // cout << event_being_read << endl;
         // cout << "Writing file " << output_path << "/" << correct_filename << ".mod" << endl;
         output_file << event_being_read;

         number_of_events_written++;
   
      }

      // Add it to the map completed_events and to the file.
      completed_events.insert(make_pair(to_string(event_being_read.run_number()) + "_" + to_string(event_being_read.event_number()), 1));

      completed_events_file_output << event_being_read.run_number() << "\t" << event_being_read.event_number() << endl; 

   }
   else {
      cout << "Event " << event_being_read.run_number() << "\t" << event_being_read.event_number() << " already processed so skipping." << endl;
      log_file << "Event " << event_being_read.run_number() << "\t" << event_being_read.event_number() << " is a duplicate so skipping." << endl;
   }
   
}

string find_correct_file_for_event(int run_number, int event_number, unordered_map<string, string> & registry_info) {

   auto search = registry_info.find(to_string(run_number) + "_" + to_string(event_number));
   
   if(search != registry_info.end()) {
      // std::cout << "Found " << search->first << " " << search->second << '\n';
      return search->second;
   }
   else {
      // cout << "Oops didn't find it." << endl;
   }

   return "";
}





void get_all_files_to_process(std::vector<string> & all_files, boost::filesystem::path input_path) {
   
   directory_iterator end_itr;

   for (directory_iterator itr(input_path); itr != end_itr; ++itr) {
      
      if (is_regular_file(itr->path())) {
         string current_file = itr->path().string();
         
         if (current_file.substr( current_file.length() - 3, current_file.length()) == "mod") {
            all_files.push_back(current_file);   
         }

      }
      else {
         // cout << itr->path().string() << endl;
         get_all_files_to_process(all_files, itr->path());
      }
   }

}






void load_registry_information(string registry_filename, unordered_map<string, string> & registry_info) {

   cout << "Loading registry to memory." << endl;

   ifstream registry(registry_filename);

   int line_number = 1;

   string line;
   // while ((getline(registry, line)) && (line_number < 100000000)) {
   while (getline(registry, line)) {

      if (line_number % 100000 == 0)
         cout << "On line number " << line_number << endl;

      istringstream iss(line);

      int registry_event_number, registry_run_number;
      string root_filename;

      iss >> registry_event_number >> registry_run_number >> root_filename;

      registry_info.emplace(to_string(registry_run_number) + "_" + to_string(registry_event_number), root_filename.substr(0, 36));  // 36 because the filenames without extensions are 37 characters long.

      line_number++;
   }
}