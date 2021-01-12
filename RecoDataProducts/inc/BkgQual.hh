// 
// Detail class used to implement the MVAStruct used to qualify low-energy electron spirals in the tracker
// D. Brown, LBNL 3/28/2017
//
#ifndef RecoDataProducts_BkgQual_hh
#define RecoDataProducts_BkgQual_hh
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
  struct BkgQualDetail {
// enumerate the input varibles used in BkgQual.  The order should match that used in
// the MVA configuration (XML).  The names may optionally be required to match exactly what's in the
// MVA configuration
    enum MVA_varindex {hrho=0, shrho, crho, sdt, zmin, zmax, zgap, np, npexp, npfrac, nphits, nhits, sfrac, n_vars};
    typedef std::map<std::string,MVA_varindex> map_type;
    static std::string const& typeName();
    static std::map<std::string,MVA_varindex> const& varNames();
    static std::string varName(MVA_varindex vindex);
  };
  // define the types
  typedef MVAStruct<BkgQualDetail> BkgQual;
  typedef std::vector<BkgQual> BkgQualCollection;
}
#endif
