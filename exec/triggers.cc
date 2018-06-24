#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <iterator>
#include <iomanip>
#include <limits>
#include <chrono>




#include <iostream>

#include <unordered_map>




#include <algorithm>



#include <map>




using namespace std;


// void analyze_event(MOD::Event & event_being_read, ofstream & output_file, int & event_serial_number);
void read_file(istream & data_stream, unordered_map<string, float> & trigger_prescales, unordered_map<string, int> & trigger_numbers, string output_filename);
void write_average_prescales(string output_filename, unordered_map<string, float> & trigger_prescales, unordered_map<string, int> & trigger_numbers);

int main(int argc, char * argv[]) {
    
    auto start = std::chrono::steady_clock::now();
    
    
    if (argc <= 2) {
        std::cerr << "ERROR: You need to supply three arguments- first, path to the input data; second, path to the output file; third, number of events to process. The path has to be either absolute or relative to the bin directory:" << std::endl << std::endl << "./analysis (input_file.dat) (output_file.dat) [optional Nev]" << std::endl;
        return 1;
    }
    
    
    ifstream data_file(argv[1]);
    // ofstream output_file(argv[2], ios::out);
    string output_filename = argv[2];
    
    
    cout << endl << endl << "Starting analysis with the following given arguments: " << endl;
    cout << "Input file: " << argv[1] << endl;
    cout << "Output file: " << argv[2] << endl;
    
    
    unordered_map<string, float> trigger_prescales;
    unordered_map<string, int> trigger_numbers;
    
    read_file(data_file, trigger_prescales, trigger_numbers, output_filename);
    write_average_prescales(output_filename, trigger_prescales, trigger_numbers);
    
    auto finish = std::chrono::steady_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start).count();
    
    cout << "Finished processing in " << elapsed_seconds << " seconds!" << endl;
    
    return 0;
}


void write_average_prescales(string output_filename, unordered_map<string, float> & trigger_prescales, unordered_map<string, int> & trigger_numbers) {
    
    ofstream output_stream(output_filename, ios::out);
    
    for ( auto it = trigger_prescales.begin(); it != trigger_prescales.end(); ++it ) {
        output_stream << it-> first << "\t" << it->second << "\t" << trigger_numbers[it->first] << endl;
    }
    
}

void read_file(istream & data_stream, unordered_map<string, float> & trigger_prescales, unordered_map<string, int> & trigger_numbers, string output_filename) {
    string line;
    
    int line_number = 0;
    
    while(getline(data_stream, line)) {
        
        line_number++;
        
        if (line_number == 1000000) {
            cout << "Finished 1M lines." << endl;
            write_average_prescales(output_filename, trigger_prescales, trigger_numbers);
            line_number = 0;
        }
        
        istringstream iss(line);
        
        string tag;
        int event_number, run_number;
        int trig_jet_matched, jet_quality;
        float hardest_pT, corr_hardest_pT, hardest_eta;
        int prescale;
        string trigger_name;
        int trigger_fired;
        
        iss >> tag;
        istringstream stream(line);
        
        
        if ((line.empty()) || (tag == "#")) {
            continue;
        }
        else {
            
            stream >> tag >> event_number >> run_number >> trig_jet_matched >> jet_quality >> hardest_pT >> corr_hardest_pT >> hardest_eta >> prescale >> trigger_name >> trigger_fired;
            
            if (trig_jet_matched == 0) {
                cout << "Trig. not matched!" << endl;
            }
            else {
                
                    unordered_map<string, float>::const_iterator got = trigger_prescales.find(trigger_name);
                    
                    if ( got == trigger_prescales.end() ) {
                        // Trigger name is not in the hashmap yet.
                        trigger_prescales.emplace(trigger_name, prescale);
                        trigger_numbers.emplace(trigger_name, 1);
                    }
                    else {
                        // Trigger name already in the hashmap.
                        // Get the total number of prescales already in the hashmap.
                        int n = trigger_numbers[trigger_name];
                        float summation_x = n * trigger_prescales[trigger_name];
                        float new_mean = (summation_x + prescale) / (n + 1);
                        
                        trigger_prescales[trigger_name] = new_mean;
                        trigger_numbers[trigger_name]++;
                    }
            }
            
            
            
        }
    }
}
