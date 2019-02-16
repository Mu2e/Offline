//
// struct to record event weight information
// Andy Edmonds, September 2019
//
#ifndef EventWeightInfo_HH
#define EventWeightInfo_HH
#include "Rtypes.h"
#include <string>
namespace mu2e
{
  struct EventWeightInfo {
    static const int MAX_WEIGHTS = 50;
    const std::string leafnames(std::vector<std::string> labels) {
      std::string leaves = "nwts/I:";
      for (std::vector<std::string>::const_iterator i_label = labels.begin(); i_label != labels.end(); ++i_label) {
	leaves += *i_label + "/F";
	if (i_label != labels.end()-1) {
	  leaves += ":";
	}
      }
      n_weights = labels.size();
      return leaves;
    }

    void setWeights(const std::vector<Float_t>& weights) { 
      for (unsigned int i_weight = 0; i_weight < weights.size(); ++i_weight) {
	_weights[i_weight] = weights.at(i_weight);
      }
    }

    void reset() {
      for (auto& i_weight : _weights) {
	i_weight = -1.0;
      }
    }

    Int_t n_weights;
    Float_t _weights[MAX_WEIGHTS];
  };
}
#endif
