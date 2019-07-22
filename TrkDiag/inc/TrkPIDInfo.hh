//
// struct to place variables that are used for TrkPID calculation
// 
#ifndef TrkPIDInfo_HH
#define TrkPIDInfo_HH
#include "RecoDataProducts/inc/TrkCaloHitPID.hh"
#include "Rtypes.h"
namespace mu2e
{
// general information about a track
  struct TrkPIDInfo {
    Float_t _tchpvars[TrkCaloHitPID::n_vars];
    Float_t _mvaout;
    Int_t _mvastat;
// extrapolation radii to the disks; front and back
    Float_t _diskfrad[2], _diskbrad[2];

    // eventually add dE/dx PID FIXME!

    TrkPIDInfo() { reset(); }

    void reset() {
      int n_tchp_vars = TrkCaloHitPID::n_vars;
      for (int ivar = 0; ivar < n_tchp_vars; ++ivar) {
	_tchpvars[ivar] = -1000.0;
      }
      _mvaout = -1000.0;
      _mvastat = -1;
      for(size_t idisk=0;idisk<2;idisk++){
	_diskfrad[idisk]=-1000.0;
	_diskbrad[idisk]=-1000.0;
      }
    }

    static std::string const leafnames() { 
      std::string leaves = "";
      for (int ivar = 0; ivar < TrkCaloHitPID::n_vars; ++ivar) {
	TrkCaloHitPID::MVA_varindex index =TrkCaloHitPID::MVA_varindex(ivar);
	std::string varname = TrkCaloHitPID::varName(index);
	leaves += varname + "/F:";
      }
      leaves += "mvaout/F:mvastat/I:disk0frad/F:disk1frad/F:disk0brad/F:disk1brad/F";
      return leaves;
    }
  };
}
#endif

