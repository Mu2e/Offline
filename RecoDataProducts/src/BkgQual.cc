// Detail class used to implement the MVAStruct for BkgQual, used to
// estimate the 'quality' of a background cluster (low-E electron)
// D. Brown, LBNL 1/30/2017
#include "RecoDataProducts/inc/BkgQual.hh"
#include <stdexcept>
#include <iostream>
#include <stdio.h>

namespace mu2e {
  std::string const& BkgQualDetail::typeName() {
    static std::string type("BkgQual");
    return type;
  }
  std::map<std::string,BkgQualDetail::MVA_varindex> const& BkgQualDetail::varNames() {
    static std::map<std::string,BkgQualDetail::MVA_varindex> varnames;
    if(varnames.size()==0){
      varnames[std::string("HitRho")]		= BkgQualDetail::hrho;
      varnames[std::string("HitRhoSpread")]	= BkgQualDetail::shrho;
      varnames[std::string("ClusterRho")]	= BkgQualDetail::crho;
      varnames[std::string("TimeSpread")]	= BkgQualDetail::sdt;
      varnames[std::string("ZMin")]		= BkgQualDetail::zmin;
      varnames[std::string("ZMax")]		= BkgQualDetail::zmax;
      varnames[std::string("ZGap")]		= BkgQualDetail::zgap;
      varnames[std::string("NPlanes")]		= BkgQualDetail::np;
      varnames[std::string("NExpectedPlanes")]	= BkgQualDetail::npexp; 
      varnames[std::string("PlaneFraction")]	= BkgQualDetail::npfrac; 
      varnames[std::string("NPlaneHits")]	= BkgQualDetail::nphits; 
      varnames[std::string("NHits")]		= BkgQualDetail::nhits; 
      varnames[std::string("StereoFraction")]	= BkgQualDetail::sfrac; 
    }
    return varnames;
  }

  std::string BkgQualDetail::varName(MVA_varindex vindex) {
    static std::string nullret("Unknown");
    std::map<std::string,BkgQualDetail::MVA_varindex> const& varnames = varNames();
    for(const auto& ivar : varnames) {
      if(ivar.second == vindex){
	return ivar.first;
      }
    }
    return nullret;
  }
}
