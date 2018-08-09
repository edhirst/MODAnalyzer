#include "InfoPFC.h"

using namespace std;
using namespace fastjet;


MOD::InfoPFC::InfoPFC(int pdgId, string tag, int vertex) : _pdgId(pdgId), _tag(tag), _vertex(vertex){}

const int MOD::InfoPFC::pdgId() const {
  return _pdgId;
}

const string MOD::InfoPFC::tag() const {
  return _tag;
}

const int MOD::InfoPFC::vertex() const {
    return _vertex;
}


const std::string MOD::InfoPFC::header() const {
	stringstream ss;
	ss << "#    " << _tag << "              px              py              pz          energy           pdgId" << endl;
	return ss.str();
}
