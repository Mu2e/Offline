// Detail class used to implement the MVAStruct for TrkCaloHitPID
#include "RecoDataProducts/inc/TrkCaloHitPID.hh"
#include <stdexcept>
#include <iostream>
#include <stdio.h>

namespace mu2e {
  std::string const& TrkCaloHitPIDDetail::typeName() {
    static std::string type("TrkCaloHitPID");
    return type;
  }
  std::map<std::string,TrkCaloHitPIDDetail::MVA_varindex> const& TrkCaloHitPIDDetail::varNames() {
    static std::map<std::string,TrkCaloHitPIDDetail::MVA_varindex> varnames;
    if(varnames.size()==0){
      varnames[std::string("DeltaE")]	  = TrkCaloHitPIDDetail::DeltaE;
      varnames[std::string("ClusterLength")]  = TrkCaloHitPIDDetail::ClusterLen;
      varnames[std::string("RPOCA")]	  = TrkCaloHitPIDDetail::RPOCA;
      varnames[std::string("TrackDirection")]		  = TrkCaloHitPIDDetail::TrkDir;
      varnames[std::string("DeltaT")]		  = TrkCaloHitPIDDetail::DeltaT;
    }
    return varnames;
  }

  std::string TrkCaloHitPIDDetail::varName(MVA_varindex vindex) {
    static std::string nullret("Unknown");
    std::map<std::string,TrkCaloHitPIDDetail::MVA_varindex> const& varnames = varNames();
    for(const auto& ivar : varnames) {
      if(ivar.second == vindex){
	return ivar.first;
      }
    }
    return nullret;
  }
}
