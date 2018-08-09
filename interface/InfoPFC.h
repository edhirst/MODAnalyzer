#ifndef InfoPFC_H
#define InfoPFC_H

#include <iostream>

#include "fastjet/ClusterSequence.hh"

namespace MOD {

	class InfoPFC : public fastjet::PseudoJet::UserInfoBase {

		public:
			InfoPFC(int pdgId, std::string tag, int vertex);
			
			const int pdgId() const;
			const std::string tag() const;
			const std::string header() const;
            const int vertex() const;
			
		private:
			int _pdgId;
			std::string _tag;
            int _vertex;

	};
}


#endif /* InfoPFC_H */
