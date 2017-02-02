// Detail class used to implement the MVAStruct for TrkQual, the track quality
// measure
// D. Brown, LBNL 1/30/2017
#include "RecoDataProducts/inc/TrkQual.hh"
#include <stdexcept>
#include <iostream>
#include <stdio.h>

namespace mu2e {
  std::string const& TrkQualDetail::typeName() {
    static std::string type("TrkQual");
    return type;
  }
  std::map<std::string,TrkQualDetail::MVA_varindex> const& TrkQualDetail::varNames() {
    static std::map<std::string,TrkQualDetail::MVA_varindex> varnames;
    if(varnames.size()==0){
      varnames[std::string("NActiveHits")]	  = TrkQualDetail::nactive;
      varnames[std::string("ActiveHitFraction")]  = TrkQualDetail::factive;
      varnames[std::string("Log10FitCon")]	  = TrkQualDetail::log10fitcon;
      varnames[std::string("MomError")]		  = TrkQualDetail::momerr;
      varnames[std::string("T0Error")]		  = TrkQualDetail::t0err;
      varnames[std::string("d0")]		  = TrkQualDetail::d0;
      varnames[std::string("MaxRadius")]	  = TrkQualDetail::rmax; 
      varnames[std::string("DoubleHitFraction")]  = TrkQualDetail::fdouble; 
      varnames[std::string("NullAmbigHitFraction")]= TrkQualDetail::fnullambig; 
      varnames[std::string("StrawHitFraction")]   = TrkQualDetail::fstraws; 
    }
    return varnames;
  }

}
