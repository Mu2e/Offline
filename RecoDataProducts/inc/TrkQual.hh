// Detail class used to implement the MVAStruct for TrkQual, the track quality
// measure
// D. Brown, LBNL 1/30/2017
//
#ifndef RecoDataProducts_TrkQual_hh
#define RecoDataProducts_TrkQual_hh
// 
//
// Original author David Brown
//
// Mu2e includes
#include "GeneralUtilities/inc/MVAStruct.hh"
#include <string>
#include <map>
#include <vector>
namespace mu2e {
  struct TrkQualDetail {
// enumerate the input varibles used in TrkQual.  The order should match that used in
// the MVA configuration (XML).  The names may optionally be required to match exactly what's in the
// MVA configuration
    enum MVA_varindex { nactive=0, factive, log10fitcon, momerr, t0err, d0, rmax,
      fdouble, fnullambig, fstraws, n_vars};
    typedef std::map<std::string,MVA_varindex> map_type;
    static std::string const& typeName();
    static std::map<std::string,MVA_varindex> const& varNames();
    static std::string varName(MVA_varindex vindex);
  };
  // define the types
  typedef MVAStruct<TrkQualDetail> TrkQual;
  typedef std::vector<TrkQual> TrkQualCollection;

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

  std::string TrkQualDetail::varName(MVA_varindex vindex) {
    static std::string nullret("Unknown");
    std::map<std::string,TrkQualDetail::MVA_varindex> const& varnames = varNames();
    for(const auto& ivar : varnames) {
      if(ivar.second == vindex){
	return ivar.first;
      }
    }
    return nullret;
  }
}
#endif
