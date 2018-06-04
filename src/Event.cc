#include "Event.h"


using namespace std;
using namespace fastjet;

MOD::Event::Event(int run_number, int event_number, int lumi_block, double inst_lumi) : _run_number(run_number), _event_number(event_number) {}


MOD::Event::Event() {}

const int MOD::Event::event_number() const {
   return _condition.event_number();
}

const int MOD::Event::run_number() const {
   return _condition.run_number();
}

const MOD::Condition MOD::Event::condition() const {
   return _condition;
}

const int MOD::Event::version() const {
   return _version;
}

const pair<string, string> MOD::Event::data_type() const {
   return _data_type;
}


void MOD::Event::set_version(int version) {
   _version = version;
}

void MOD::Event::set_data_type(string a, string b) {
   _data_type = make_pair(a, b);;
}

const vector<PseudoJet> & MOD::Event::particles() const {
   return _particles;
}


const vector<PseudoJet> & MOD::Event::cms_jets() const {
   return _cms_jets;
}


const vector<PseudoJet> & MOD::Event::jets() const {
   return _jets;
}


void MOD::Event::add_particle(istringstream & input_stream) {
   string tag;
   double px, py, pz, energy;
   int pdgId;

   input_stream >> tag >> px >> py >> pz >> energy >> pdgId;

  //cout << "reached here" << endl;
   PseudoJet new_particle = PseudoJet(px, py, pz, energy);
   new_particle.set_user_info(new InfoPFC(pdgId, tag));
    

   _particles.push_back(new_particle);
}

void MOD::Event::add_cms_jet(istringstream & input_stream) {
   string tag;
   double px, py, pz, energy, JEC, area;

   double neutral_hadron_fraction, neutral_em_fraction, charged_hadron_fraction, charged_em_fraction;
   int number_of_constituents, charged_multiplicity;

   input_stream >> tag >> px >> py >> pz >> energy >> JEC >> area >> number_of_constituents >> charged_multiplicity >> neutral_hadron_fraction >> neutral_em_fraction >> charged_hadron_fraction >> charged_em_fraction;
   
   PseudoJet new_jet = PseudoJet(px, py, pz, energy);
   double eta = new_jet.eta();

   new_jet.set_user_info(new InfoCalibratedJet(tag, JEC, area, number_of_constituents, charged_multiplicity, neutral_hadron_fraction, neutral_em_fraction, charged_hadron_fraction, charged_em_fraction, eta));

   // We never read other jets directly so the jets that reach here will always be cms jets.
   _cms_jets.push_back(new_jet);
}


void MOD::Event::add_condition(istringstream & input_stream) {
   MOD::Condition new_condition(input_stream);
   _condition = new_condition;
}



void MOD::Event::add_trigger(istringstream & input_stream) {
   _triggers.push_back(MOD::Trigger(input_stream));
}

const MOD::Event::data_source_t MOD::Event::data_source() const {
   return _data_source;
}

void MOD::Event::set_data_source(int data_source) {
   // 0 => Experiment, 1 => MC_TRUTH, 2 => MC_RECO.
   _data_source = static_cast<data_source_t>(data_source);
}

const MOD::Trigger MOD::Event::trigger_by_name(string name) const {
   for (unsigned i = 0; i < triggers().size(); i++) {
      const MOD::Trigger& current_trigger = triggers().at(i);

      if (current_trigger.name() == name) 
         return current_trigger;
   }

   return Trigger();
}

const MOD::Trigger MOD::Event::trigger_by_short_name(string short_name) const {
   for (unsigned i = 0; i < triggers().size(); i++) {
      const MOD::Trigger& current_trigger = triggers().at(i);

      if (current_trigger.short_name() == short_name)
         return current_trigger;
   }

   return Trigger();
}

// You can give a trigger's full name or short name here.
const bool MOD::Event::trigger_exists(string trigger_name) const {
   return ( trigger_by_name(trigger_name).is_valid() || trigger_by_short_name(trigger_name).is_valid() );
}

const vector<MOD::Trigger> & MOD::Event::triggers() const {
   return _triggers;
}


const double MOD::Event::weight() const {
   //return _weight;
    // cout <<"assigned trigger" << _assigned_trigger.prescale() << "   " << _weight << endl;
   return _assigned_trigger.prescale();
}


void MOD::Event::set_weight(double weight) {
   _weight = weight;
}

void MOD::Event::set_multiplicity(int multiplicity) {
    _multiplicity = multiplicity;
}


const string MOD::Event::stringify_jet(PseudoJet jet) const {
   stringstream ss;

   if (jet.has_user_info() && jet.user_info<MOD::InfoCalibratedJet>().tag() == "AK5") {
      // This is a CMS jet.
      ss << "     " << jet.user_info<MOD::InfoCalibratedJet>().tag()
        << setw(16) << fixed << setprecision(8) << jet.px()
        << setw(16) << fixed << setprecision(8) << jet.py()
        << setw(16) << fixed << setprecision(8) << jet.pz()
        << setw(16) << fixed << setprecision(8) << jet.E()
        << setw(16) << fixed << setprecision(8) << jet.user_info<MOD::InfoCalibratedJet>().JEC()
        << setw(16) << fixed << setprecision(8) << jet.user_info<MOD::InfoCalibratedJet>().area()
        << setw(16) << fixed << setprecision(8) << jet.user_info<MOD::InfoCalibratedJet>().number_of_constituents()   
        << setw(16) << fixed << setprecision(8) << jet.user_info<MOD::InfoCalibratedJet>().charged_multiplicity()
        << setw(16) << fixed << setprecision(8) << jet.user_info<MOD::InfoCalibratedJet>().neutral_hadron_fraction()
        << setw(16) << fixed << setprecision(8) << jet.user_info<MOD::InfoCalibratedJet>().neutral_em_fraction()
        << setw(16) << fixed << setprecision(8) << jet.user_info<MOD::InfoCalibratedJet>().charged_hadron_fraction()   
        << setw(16) << fixed << setprecision(8) << jet.user_info<MOD::InfoCalibratedJet>().charged_em_fraction() 
        << endl;
   }
   

   return ss.str();
}


const string MOD::Event::stringify_pfc(PseudoJet particle) const {
   stringstream ss;
   ss << "     " << particle.user_info<MOD::InfoPFC>().tag()
        << setw(16) << fixed << setprecision(8) << particle.px()
        << setw(16) << fixed << setprecision(8) << particle.py()
        << setw(16) << fixed << setprecision(8) << particle.pz()
        << setw(16) << fixed << setprecision(8) << particle.E()
        << setw(16) << noshowpos << particle.user_info<MOD::InfoPFC>().pdgId()
        << endl;

   return ss.str();
}


string MOD::Event::make_string() const {
   stringstream file_to_write;
   
 
   file_to_write << "BeginEvent Version " << _version << " " << _data_type.first << " " << _data_type.second;  // Don't put an endl here because for "pristine", we'll put a "Hardest_Jet_Selection" here. 


   if (_data_source == EXPERIMENT) { // Data is from experiment.
      
      file_to_write << endl;

      // First, write out conditions.
      file_to_write << _condition.make_header_string();
       
     //  cout << "condition " << _condition << endl;
      file_to_write << _condition;   

      // Then, write out all triggers. 
      for(unsigned int i = 0; i < _triggers.size(); i++) {
         if (i == 0)
            file_to_write << _triggers[i].make_header_string();
         file_to_write << _triggers[i];
      }

      // Next, write out AK5 calibrated jets.
      for(unsigned int i = 0; i < _cms_jets.size(); i++) {
         if (i == 0)
            file_to_write << _cms_jets[i].user_info<MOD::InfoCalibratedJet>().header();
         file_to_write << stringify_jet(_cms_jets[i]);
      }
      
      // Finally, write out all particles.
      for (unsigned int i = 0; i < _particles.size(); i++) {
         if (i == 0)
            file_to_write << _particles[i].user_info<MOD::InfoPFC>().header();
         file_to_write << stringify_pfc(_particles[i]);
      }

   }
   else if ( (_data_source == MC_TRUTH) || (_data_source == MC_RECO) ) {
      
      file_to_write << endl;

      for (unsigned i = 0; i < _particles.size(); i++) {
         if (i == 0)
            file_to_write << _particles[i].user_info<MOD::InfoPFC>().header();
         file_to_write << stringify_pfc(_particles[i]);
      }

   }
   else if (_data_source == PRISTINE) {
      

      file_to_write << " Hardest_Jet_Selection" << endl;

      // This output is for pristine form. Data in "pristine form" is used as-it-is i.e. without any consideration for things like jet quality levels and triggers.

      // First, write out conditions.
      file_to_write << _condition.make_header_string();
      file_to_write << _condition;   
      
      // We only output the hardest jet and its constituents.
      // The reason for that is, we want to do our analysis on FastJet-clustered jets instead of on cms jets. However, we need to apply JEC to our AK5 jets. 
      // Those JEC factors are known for cms jets but not for FastJet-clustered jets. This means, we need to figure out a correspondance between cms and FastJet-clustered jets.
      // This is hard to do for jets other than the hardest one.
      // So we select the hardest cms jet, use delta R to find the corresponding FastJet-clustered jet and then just output constituents of that jet.

      // While it's possible to do that for all jets, we output just the hardest jet's constituents because all analyses are going to be performed on the hardest jet anyway.

      file_to_write << "#  1JET" << "              px              py              pz          energy             jec            area          weight" << endl;
      file_to_write  << "   1JET"
                     << setw(16) << fixed << setprecision(8) << _jets[0].px()
                     << setw(16) << fixed << setprecision(8) << _jets[0].py()
                     << setw(16) << fixed << setprecision(8) << _jets[0].pz()
                     << setw(16) << fixed << setprecision(8) << _jets[0].E()
                     << setw(16) << fixed << setprecision(8) << _jets[0].user_info<MOD::InfoCalibratedJet>().JEC()
                     << setw(16) << fixed << setprecision(8) << _jets[0].user_info<MOD::InfoCalibratedJet>().area()
                     << setw(16) << _weight
                     << endl;

      file_to_write << "# PDPFC" << "              px              py              pz          energy   pdgId" << endl;

      for (unsigned i = 0; i < _particles.size(); i++) {
         
         file_to_write << "  PDPFC"
                       << setw(16) << fixed << setprecision(8) << _particles[i].px()
                       << setw(16) << fixed << setprecision(8) << _particles[i].py()
                       << setw(16) << fixed << setprecision(8) << _particles[i].pz()
                       << setw(16) << fixed << setprecision(8) << _particles[i].E()
                       << setw(8) << noshowpos << _particles[i].user_info<MOD::InfoPFC>().pdgId()
                       << endl;
      }

   }
   else {
      throw runtime_error("Invalid data source!");
   }

   
   file_to_write << "EndEvent" << endl;

   return file_to_write.str();
}



const PseudoJet & MOD::Event::hardest_jet() const {
   
   // EXPERIMENT = 0, MC_TRUTH = 1, MC_RECO = 2, PRISTINE = 3 

   if (_data_source == EXPERIMENT)
      return _closest_fastjet_jet_to_trigger_jet;

   return _hardest_jet;
}


void MOD::Event::convert_to_one_jet() {

   PseudoJet jet = _closest_fastjet_jet_to_trigger_jet;


   jet.set_user_info(new MOD::InfoCalibratedJet("1JET", _trigger_jet.user_info<MOD::InfoCalibratedJet>().JEC(), _trigger_jet.user_info<MOD::InfoCalibratedJet>().area()));

   vector<PseudoJet> particles = jet.constituents();

   int pdgId;
   for (unsigned i = 0; i < particles.size(); i++) {
      pdgId = particles[i].user_info<MOD::InfoPFC>().pdgId();

      particles[i].set_user_info( new MOD::InfoPFC(pdgId, "PDPFC") );
   }

   vector<PseudoJet> just_one_jet{jet};

   _particles = particles;
   _jets = just_one_jet;


   _data_source = PRISTINE;   // Set the data source to "Pristine".
   _weight = _assigned_trigger.prescale();



   // cout << _jets[0].px() << endl;

   // cout << _jets[0].px() << endl;

   


}


void MOD::Event::establish_properties() {

   JetDefinition jet_def(antikt_algorithm, 0.5);
   ClusterSequence * cs = new ClusterSequence(_particles, jet_def);
   vector<PseudoJet> ak5_jets = sorted_by_pt(cs->inclusive_jets(0.0));
   _jets = ak5_jets;

   cs->delete_self_when_unused();

   // cout << data_source() << endl;

   if (data_source() == EXPERIMENT) {  // Experiment 

      // cout << "Setting trigger jet!" << endl;
      set_trigger_jet();

      // Next, find out the specific FastJet that's closest to _trigger_jet.
      // cout << "Closest Jet" << endl;
      set_closest_fastjet_jet_to_trigger_jet();

      // cout << "Trigger jet matched!" << endl;
      // set_trigger_jet_is_matched();

      // cout << "Assigned trigger!" << endl;
      set_assigned_trigger();

     // cout << "Assigned trigger is " << _assigned_trigger.name() << endl;

   }
   else if (data_source() == PRISTINE) {

   }



  
   set_hardest_jet();

   // cout << _hardest_jet.pt() << endl;
}




bool MOD::Event::read_event(istream & data_stream, std::string type_of_data) {

   string line;
   // cout << "in read event " << endl;
   while(getline(data_stream, line)) {
      istringstream iss(line);

      int version;
      double weight;
      double area;
      string tag, version_keyword, a, b, c;
      double px, py, pz, energy, jec;
       
      int cur_count;

      iss >> tag;      
      istringstream stream(line);

      if (line.empty()) {
         continue;
      }
      else {
         

         if (tag == "BeginEvent") {

            // stream >> tag >> version_keyword >> version >> a >> b;               // For MOD files.
            stream >> tag >> version_keyword >> version >> a >> b >> c >> weight;   // For pristine, I think. Or maybe skim.

            set_version(version);
            set_data_type(a, b);


            set_weight(weight);
             
            cur_count = 0;
            

         }
         else if (tag == "1JET") {
            try {
               set_data_source(PRISTINE);
               
               stream >> tag >> px >> py >> pz >> energy >> jec >> area >> weight;

               _weight = weight;

               PseudoJet jet = PseudoJet(px, py, pz, energy);
               jet.set_user_info( new MOD::InfoCalibratedJet("1JET", jec, area) );
               vector<PseudoJet> jets{jet};
               
               _cms_jets = jets;
               _jets = jets;

            }
            catch (exception& e) {
               throw runtime_error("Invalid file format 1JET! Something's wrong the way 1JETs have been written!");
            }
         }
         else if (tag == "PFC" and type_of_data == "2011") {
            try {
               set_data_source(EXPERIMENT);
               add_particle(stream);
            }
            catch (exception& e) {
               throw runtime_error("Invalid file format PFC! Something's wrong with the way PFCs have been written.");
            }
         }
         else if (tag == "PFC" and type_of_data == "sim_pfc") {
              try {
                  set_data_source(MC_RECO);
                  add_particle(stream);
              }
              catch (exception& e) {
                  throw runtime_error("Invalid file format PFC! Something's wrong with the way PFCs have been written.");
              }
         }
         else if ((tag == "GEN" or tag == "Gen") and type_of_data == "sim_gen") {
              try {
                  set_data_source(MC_TRUTH);
                  add_particle(stream);
              }
              catch (exception& e) {
                  throw runtime_error("Invalid file format PFC! Something's wrong with the way PFCs have been written.");
              }
         }
         else if (tag == "AK5") {
            try {
               set_data_source(EXPERIMENT);
               // cout << "Adding AK5" << endl;
               add_cms_jet(stream);
            }
        
            catch (exception& e) {
               throw runtime_error("Invalid file format AK5! Something's wrong with the way jets have been written.");
            }
         }
         else if (tag == "Trig") {
            try {
               // cout << "Trig" << endl;
               add_trigger(stream);
            }
            catch (exception& e) {
               throw runtime_error("Invalid file format! Something's wrong with the way triggers have been written.");
            }
         }
         else if (tag == "Cond") {
            try {
               // cout << "Cond" << endl;
               add_condition(stream);
            }
            catch (exception& e) {
               throw runtime_error("Invalid file format! Something's wrong with the way conditions have been written.");
            }
         }
         else if (tag == "TRUTH") {
            try {
               set_data_source(MC_TRUTH);
               add_particle(stream);
            }
            catch (exception& e) {
               throw runtime_error("Invalid file format! Something's wrong with the way Truth has been written. ;)");
            }
         }
         else if (tag == "RPFC") {
            try {
               set_data_source(MC_RECO);
               add_particle(stream);
            }
            catch (exception& e) {
               throw runtime_error("Invalid file format! Something's wrong with the way RPFC has been written. ;)");
            }
         }
         else if (tag == "PDPFC") {
            try {
               set_data_source(PRISTINE);
               add_particle(stream);
               cur_count++;
            }
            catch (exception& e) {
               throw runtime_error("Invalid file format! Something's wrong with the way RPFC has been written. ;)");
            }
         }
         else if (tag == "EndEvent") {
            
            // cout << "EndEvent" << endl;

            set_multiplicity(cur_count);
            establish_properties();
            return true;
         }
         else if (tag == "#" or ((tag == "GEN" or tag == "Gen") and type_of_data == "sim_pfc") or (tag == "PFC" and type_of_data == "sim_gen") ) {
            // Either this line in the data file represents a comment
            // Or this line is a sim_gen_data line and we are processing sim_pfc_data
            // Or this line is a sim_pfc_data line and we are processing sim_gen_data
            // Just ignore it.
         }
         else {
            throw runtime_error("Invalid file format! Unrecognized tag: " + tag + "!");
         }
      }
      
   }
    
   // cout << "processed error " << endl;

   return false;
}

string MOD::Event::assigned_trigger_name() const {
   return _assigned_trigger_name;
}

const MOD::Trigger MOD::Event::assigned_trigger() const {
   return _assigned_trigger;
}

bool MOD::Event::assigned_trigger_fired() const {
   // cout << _assigned_trigger.is_valid() << ", " << _assigned_trigger.fired() << endl;
   bool fired = _assigned_trigger.is_valid() && _assigned_trigger.fired();

   // cout <<  _assigned_trigger.is_valid() << _assigned_trigger.fired() << endl;
   return fired;
}

int MOD::Event::assigned_trigger_prescale() const {
   return _assigned_trigger.prescale();
}

void MOD::Event::set_assigned_trigger() {
   
   string trigger_to_use = "";
   string trigger = "";

   // First, figure out the hardest pT.
   PseudoJet trigger_jet = _trigger_jet;
   
   //cout << "JEC is: " << _hardest_jet_jec << endl;

   trigger_jet *= _hardest_jet_jec;
    
  

   if (trigger_jet.E() == 0.0) {
      _assigned_trigger_name = "";
      _assigned_trigger = Trigger();
      // cout << "energy is zero" << endl;
   }
   else {

      double hardest_pT = trigger_jet.pt();

      // cout << "hardest_pT" << hardest_pT << endl;

      // Next, lookup which trigger to use based on the pT value of the hardest jet.

      if ( (hardest_pT > 420) && trigger_exists("HLT_Jet370") ) {
         trigger = "HLT_Jet370";
      }
      else if ( (hardest_pT > 350) && trigger_exists("HLT_Jet300") ) {
         trigger = "HLT_Jet300";
      }
      else if ( (hardest_pT > 275) && trigger_exists("HLT_Jet240") ) {
         trigger = "HLT_Jet240";
      }
      else if ( (hardest_pT > 220) && trigger_exists("HLT_Jet190") ) {
         trigger = "HLT_Jet190";
      }
      else if ( (hardest_pT > 190) && trigger_exists("HLT_Jet150") ) {
         trigger = "HLT_Jet150";
      }
      else if ( (hardest_pT > 130) && trigger_exists("HLT_Jet110") ) {
           trigger = "HLT_Jet110";
       }
       else if ( (hardest_pT > 110) && trigger_exists("HLT_Jet80") ) {
           trigger = "HLT_Jet80";
       }
       else if ( (hardest_pT > 70) && trigger_exists("HLT_Jet60") ) {
           trigger = "HLT_Jet60";
       }
       else if (trigger_exists("HLT_Jet30")) {
           trigger = "HLT_Jet30";
       }

      trigger_to_use = trigger;

      if (trigger_to_use != "") {
         _assigned_trigger_name = trigger_to_use;
         _assigned_trigger = trigger_by_short_name(trigger_to_use);   
      }
      else {
         _assigned_trigger_name = "";
         _assigned_trigger = Trigger();
      }
   }

   
   // cout << "Here: " << _assigned_trigger.name() << endl;

}


double MOD::Event::get_hardest_jet_jec() const {
   // cout << typeid(sorted_by_pt(_jets)[0]).name() << endl;
   //  return sorted_by_pt(_jets)[0].user_info<InfoCalibratedJet>().JEC();
   return _hardest_jet_jec;
}

double MOD::Event::get_hardest_jet_area() const {
    // TODO: Area needs to be fixed.
    // return apply_jet_energy_corrections(_cms_jets).area();
    // return sorted_by_pt(_jets)[0].user_info<InfoCalibratedJet>().area();
    return sorted_by_pt(_cms_jets)[0].user_info<InfoCalibratedJet>().area();
}

void MOD::Event::set_trigger_jet() {
	
   /*
   for (unsigned i=0; i < _cms_jets.size(); i++) {
	cout << _cms_jets[i].pt() << ", ";
   }
	
   cout << endl << "========================================" << endl;
   */

   // TODO: THis selects trigger jet without JEC.

   // Get the hardest jet, apply JEC, and then eta cut.

   vector<PseudoJet> processed_jets = apply_jet_energy_corrections(_cms_jets);
   // Selector rapidity_selector = SelectorAbsRapMax(2.4);
   // processed_jets = rapidity_selector(processed_jets);

   /*
   for (unsigned i=0; i < processed_jets.size(); i++) {
	cout << processed_jets[i].pt() << ", ";
   }
	
   cout << endl << "========================================" << endl;
   

   for (unsigned i=0; i < processed_jets.size(); i++) {
	cout << processed_jets[i].user_info<InfoCalibratedJet>().JEC() << ", ";
   }
	
   cout << endl << "========================================" << endl;



   cout << processed_jets[0].pt() << endl;
   */

   // cout << processed_jets.size() << endl;

   // Then, sort the jets and store the hardest one as _trigger_jet.
   if (processed_jets.size() > 0) {
      vector<PseudoJet> sorted_processed_jets = sorted_by_pt(processed_jets);
      
      //cout << sorted_processed_jets[0].pt() << endl;

      int index = find(processed_jets.begin(), processed_jets.end(), sorted_processed_jets[0]) - processed_jets.begin();
      //auto index = distance(processed_jets.begin(), it);
      _trigger_jet = _cms_jets[index];

      //cout << "index is " << index << endl;

      //cout << "in processed jet it is, " << processed_jets[index].pt() << endl;

      //cout << "Trigger jet is " << _trigger_jet.pt() << endl;
   }
   else {
      _trigger_jet = PseudoJet();
   }
}


const fastjet::PseudoJet & MOD::Event::closest_fastjet_jet_to_trigger_jet() const {
   return _closest_fastjet_jet_to_trigger_jet;
}

void MOD::Event::set_closest_fastjet_jet_to_trigger_jet() {

   if (_trigger_jet.E() > 0) {   // This ensures that trigger jet is a valid jet and not an empty PseudoJet().

      // Next, check if the trigger jet has a rapidity of < 2.4.
      
      // if ( abs(_trigger_jet.eta()) < 2.4 ) {
         vector<PseudoJet> fastjet_jets = _jets;

         // Loop through all FastJet pseudojets, calculating delta R for each one.
         vector<double> delta_Rs;
         for (unsigned i = 0; i < fastjet_jets.size(); i++) {
            delta_Rs.push_back( _trigger_jet.delta_R(fastjet_jets[i]) );
         }
         
         // Find the index of the fastjet jet that has the lowest delta R. This will be the jet that's "closest" to the CMS jet. 
         int index = -1;
         double delta_R = numeric_limits<double>::max();
         for (unsigned i = 0; i < delta_Rs.size(); i++) {
            if (delta_Rs[i] <= delta_R) {
               delta_R = delta_Rs[i];
               index = i;
            }
         }

         if (index >= 0) {

            // We now have the corresponding "hardest" FastJet jet.
            _closest_fastjet_jet_to_trigger_jet = fastjet_jets[index];
            return;   
         }

      // }

   }

   _closest_fastjet_jet_to_trigger_jet = PseudoJet();
   return;
}


const PseudoJet & MOD::Event::trigger_jet() const {
   return _trigger_jet;
}

// void MOD::Event::set_trigger_jet_is_matched() {


//    // cout << "Energy: " << _trigger_jet.E() << endl;

//    if (_trigger_jet.E() == 0.0) {
//       _trigger_jet = PseudoJet();
//       _trigger_jet_is_matched = false;
//       cout << "Not matched " << _event_number << endl;
//       return;
//    }


//    // Compare the number of constituents first.
//    if ((unsigned) _trigger_jet.user_info<MOD::InfoCalibratedJet>().number_of_constituents() != (unsigned) _closest_fastjet_jet_to_trigger_jet.constituents().size()) {
//       // _trigger_jet = PseudoJet();
//       _trigger_jet_is_matched = false;
//       cout << "Not matched " << _event_number << endl;
//       return;
//    }


//    // Next, compare if the 4-vector matches upto 10e-4 precision or not.
//    double tolerance = pow(10., -3.);
//    if ( ( abs(_trigger_jet.px() - _closest_fastjet_jet_to_trigger_jet.px()) < tolerance ) && ( abs(_trigger_jet.py() - _closest_fastjet_jet_to_trigger_jet.py()) < tolerance ) && ( abs(_trigger_jet.pz() - _closest_fastjet_jet_to_trigger_jet.pz()) < tolerance ) && ( abs(_trigger_jet.E() - _closest_fastjet_jet_to_trigger_jet.E()) < tolerance ) ) {      
//       _trigger_jet_is_matched = true;
//       return;
//    }



//    _trigger_jet_is_matched = false;
//    // _trigger_jet = PseudoJet();

//    cout << "Not matched " << _event_number << endl;

//    return;
// }



bool MOD::Event::is_trigger_jet_matched() {


   // cout << "Energy: " << _trigger_jet.E() << endl;

   if (_trigger_jet.E() == 0.0) {
      return false;
   }


   // Compare the number of constituents first.
   if ((unsigned) _trigger_jet.user_info<MOD::InfoCalibratedJet>().number_of_constituents() != (unsigned) _closest_fastjet_jet_to_trigger_jet.constituents().size()) {
      return false;
   }


   // Next, compare if the 4-vector matches upto 10e-4 precision or not.
   double tolerance = pow(10., -3.);
   if ( ( abs(_trigger_jet.px() - _closest_fastjet_jet_to_trigger_jet.px()) < tolerance ) && ( abs(_trigger_jet.py() - _closest_fastjet_jet_to_trigger_jet.py()) < tolerance ) && ( abs(_trigger_jet.pz() - _closest_fastjet_jet_to_trigger_jet.pz()) < tolerance ) && ( abs(_trigger_jet.E() - _closest_fastjet_jet_to_trigger_jet.E()) < tolerance ) ) {      
      return true;
   }


   return false;
   
}





void MOD::Event::set_hardest_jet() {
   _hardest_jet = sorted_by_pt(_cms_jets)[0];
}



vector<PseudoJet> MOD::Event::apply_jet_energy_corrections(vector<PseudoJet> jets) {

   vector<PseudoJet> jec_corrected_jets;

   for (unsigned i = 0; i < jets.size(); i++) {
      jec_corrected_jets.push_back( jets[i] * jets[i].user_info<InfoCalibratedJet>().JEC() );
   }

   if (jets.size() > 0) {
      _hardest_jet_jec = sorted_by_pt(jec_corrected_jets)[0].user_info<InfoCalibratedJet>().JEC();
   }

   return jec_corrected_jets;
}






namespace MOD {
   ostream& operator<< (ostream& os, const Event& event) {
      os << event.make_string();
      return os;
   }
}
