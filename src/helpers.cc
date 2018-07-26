#include "helpers.h"
#include <algorithm>

using namespace std;
using namespace fastjet;




std::vector<fastjet::PseudoJet> MOD::filter_charged(std::vector<fastjet::PseudoJet> particles) {
  vector<PseudoJet> filtered;

  vector<int> pdg_ids { 1, -1, 2, -2, -11, -12, -13, -14, -16, -211, -2112, -2212, -321, 11, 12, 13, 130, 14, 16, 211, 2112, 22, 2212, 321, 310, 3222, 3122, 3312, 3322, -310, -3222, -3122, -3312, -3322, 3112, -3112, 3334, -3334};

  vector<int> charged_pdg_ids {1, 2, 11, 13, 15, 211, 321, 2212, 3112, 3222, 3312, 3334};

  for (unsigned i = 0; i < particles.size(); i++) {
    // if ( (abs(particles[i].user_info<MOD::InfoPFC>().pdgId()) == 211) || (abs(particles[i].user_info<MOD::InfoPFC>().pdgId()) == 11) || (abs(particles[i].user_info<MOD::InfoPFC>().pdgId()) == 13) ) {

    int pdgId = abs( particles[i].user_info<MOD::InfoPFC>().pdgId() );

    std::vector<int>::iterator it = find (charged_pdg_ids.begin(), charged_pdg_ids.end(), pdgId);

	if (it != charged_pdg_ids.end()){
		filtered.push_back(particles[i]);

  }

  }

  return filtered;
}
