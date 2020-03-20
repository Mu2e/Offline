//
// struct to record RecoQual information
// Andy Edmonds, September 2019
//
#ifndef RecoQualInfo_HH
#define RecoQualInfo_HH
#include "Rtypes.h"
#include <string>
namespace mu2e
{
  struct RecoQualInfo {
    static const int MAX_QUALS = 50;
    const std::string leafnames(std::vector<std::string> labels) {
      std::string leaves = "nquals/I:";
      for (std::vector<std::string>::const_iterator i_label = labels.begin(); i_label != labels.end(); ++i_label) {
	leaves += *i_label + "/F:" + *i_label + "Eff/F";
	if (i_label != labels.end()-1) {
	  leaves += ":";
	}
      }
      n_quals = labels.size()*2;
      return leaves;
    }

    void setQuals(const std::vector<Float_t>& qualsAndEffs) { 
      for (unsigned int i_qual = 0; i_qual < qualsAndEffs.size(); ++i_qual) {
	_qualsAndEffs[i_qual] = qualsAndEffs.at(i_qual);
      }
    }

    void reset() {
      for (auto& i_qualAndEff : _qualsAndEffs) {
	i_qualAndEff = -1.0;
      }
    }

    Int_t n_quals;
    Float_t _qualsAndEffs[MAX_QUALS];
  };
}
#endif
