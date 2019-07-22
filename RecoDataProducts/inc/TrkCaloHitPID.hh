// Detail class used to implement the MVAStruct for TrkCaloHitPID, a PID
// measure based on a TrkCaloHit
// D. Brown, LBNL 5/28/2019
//
#ifndef RecoDataProducts_TrkCaloHitPID_hh
#define RecoDataProducts_TrkCaloHitPID_hh
//
// Mu2e includes
#include "GeneralUtilities/inc/MVAStruct.hh"
#include <string>
#include <map>
#include <vector>
namespace mu2e {
  struct TrkCaloHitPIDDetail {
// enumerate the input varibles used in TrkCaloHitPID.  The order should match that used in
// the MVA configuration (XML).  The names may optionally be required to match exactly what's in the
// MVA configuration
    enum MVA_varindex { DeltaE= 0,ClusterLen, RPOCA, TrkDir, DeltaT , n_vars};
    typedef std::map<std::string,MVA_varindex> map_type;
    static std::string const& typeName();
    static std::map<std::string,MVA_varindex> const& varNames();
    static std::string varName(MVA_varindex vindex);
  };
  // define the types
  typedef MVAStruct<TrkCaloHitPIDDetail> TrkCaloHitPID;
  typedef std::vector<TrkCaloHitPID> TrkCaloHitPIDCollection;
}
#endif
