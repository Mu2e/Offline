// 
// Detail class used to implement the MVAStruct used to qualify low-energy electron spirals in the tracker
// D. Brown, LBNL 3/28/2017
//
#ifndef RecoDataProducts_BktClusterQual_hh
#define RecoDataProducts_BktClusterQual_hh
//
// $Id: BktClusterQual.hh,v 1.4 2013/04/04 01:08:20 brownd Exp $
// $Author: brownd $
// $Date: 2013/04/04 01:08:20 $
//
// Original author David Brown
//
// Mu2e includes
#include "GeneralUtilities/inc/MVAStruct.hh"
#include <string>
#include <map>
#include <vector>
namespace mu2e {
  struct BktClusterQualDetail {
// enumerate the input varibles used in BktClusterQual.  The order should match that used in
// the MVA configuration (XML).  The names may optionally be required to match exactly what's in the
// MVA configuration
    enum MVA_varindex {prho=0, srho, zmin, zmax, zgap, ns, nsmiss, sphi, nhits};
    typedef std::map<std::string,MVA_varindex> map_type;
    static std::string const& typeName();
    static std::map<std::string,MVA_varindex> const& varNames();
    static std::string varName(MVA_varindex vindex);
  };
  // define the types
  typedef MVAStruct<BktClusterQualDetail> BktClusterQual;
  typedef std::vector<BktClusterQual> BktClusterQualCollection;
}
#endif
