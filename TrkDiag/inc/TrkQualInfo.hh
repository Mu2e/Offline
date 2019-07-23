//
// struct to place variables that are used for TrkQual calculation
// 
#ifndef TrkQualInfo_HH
#define TrkQualInfo_HH
#include "RecoDataProducts/inc/TrkQual.hh"
#include "Rtypes.h"
namespace mu2e
{
// general information about a track
  struct TrkQualInfo {
    Float_t _trkqualvars[TrkQualDetail::n_vars];
    Float_t _mvaout;
    Int_t _mvastat;

    TrkQualInfo() { reset(); }

    void reset() {
      int n_trkqual_vars = TrkQual::n_vars;
      for (int i_trkqual_var = 0; i_trkqual_var < n_trkqual_vars; ++i_trkqual_var) {
	_trkqualvars[i_trkqual_var] = -1000.0;
      }
      _mvaout = -1000.0;
      _mvastat = -1;
    }

    static std::string const leafnames() { 
      std::string leaves = "";
      int n_trkqual_vars = TrkQual::n_vars;
      for (int i_trkqual_var = 0; i_trkqual_var < n_trkqual_vars; ++i_trkqual_var) {
	TrkQual::MVA_varindex i_index =TrkQual::MVA_varindex(i_trkqual_var);
	std::string varname = TrkQual::varName(i_index);
	leaves += varname + "/F:";
      }
      leaves += "mvaout/F:mvastat/I";
      return leaves;
    }
  };
}
#endif

