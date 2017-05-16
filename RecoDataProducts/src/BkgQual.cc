// Detail class used to implement the MVAStruct for BktClusterQual, used to
// estimate the 'quality' of a background cluster (low-E electron)
// D. Brown, LBNL 1/30/2017
#include "RecoDataProducts/inc/BktClusterQual.hh"
#include <stdexcept>
#include <iostream>
#include <stdio.h>

namespace mu2e {
  std::string const& BktClusterQualDetail::typeName() {
    static std::string type("BktClusterQual");
    return type;
  }
  std::map<std::string,BktClusterQualDetail::MVA_varindex> const& BktClusterQualDetail::varNames() {
    static std::map<std::string,BktClusterQualDetail::MVA_varindex> varnames;
    if(varnames.size()==0){
      varnames[std::string("RhoPosition")]	  = BktClusterQualDetail::prho;
      varnames[std::string("RhoSpread")]  = BktClusterQualDetail::srho;
      varnames[std::string("MinZ")]	  = BktClusterQualDetail::zmin;
      varnames[std::string("MaxZ")]		  = BktClusterQualDetail::zmax;
      varnames[std::string("ZGap")]		  = BktClusterQualDetail::zgap;
      varnames[std::string("NStations")]		  = BktClusterQualDetail::ns;
      varnames[std::string("MissingStations")]	  = BktClusterQualDetail::nsmiss; 
      varnames[std::string("PhiSpread")]  = BktClusterQualDetail::sphi; 
      varnames[std::string("NHits")]= BktClusterQualDetail::nhits; 
    }
    return varnames;
  }

  std::string BktClusterQualDetail::varName(MVA_varindex vindex) {
    static std::string nullret("Unknown");
    std::map<std::string,BktClusterQualDetail::MVA_varindex> const& varnames = varNames();
    for(const auto& ivar : varnames) {
      if(ivar.second == vindex){
	return ivar.first;
      }
    }
    return nullret;
  }
}
