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
#include "Offline/GeneralUtilities/inc/MVAStruct.hh"
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
    constexpr static const char* cxname = "TrkQualCatalog";
  };
  // define the types
  typedef MVAStruct<TrkQualDetail> TrkQual;
  typedef std::vector<TrkQual> TrkQualCollection;

}
#endif
