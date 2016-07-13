//
// structs used to record per-track calorimeter information
// 
#ifndef TrkCaloInfo_HH
#define TrkCaloInfo_HH
#include "DataProducts/inc/threevec.hh"
#include "Rtypes.h"
#include <string>
namespace mu2e
{
  using std::string;
  // struct for Calorimeter information associated with a track
  struct TrkCaloInfo {
  // construct from a track-calo match
    Float_t _dt; // extrapolated track - cluster time
    Float_t _du; // track-cluster match separation in 'u' direction
    Float_t _dv; // track-cluster match separation in 'v' direction
    Float_t _ds; // path through crystal (?) 
    Float_t _ep; // track-cluster match E/p
    Float_t _uvChisq; // chisquared from UV match
    Float_t _tChisq; // chisquared from time match
    Float_t _dtllr; // log likelihood ratio for this dt (electron/muon??)
    Float_t _epllr; // log likelihood ratio for this e/p
    Float_t _eclust; // cluster energy
    Float_t _tclust; // cluster time
    threevec _cpos; // calorimeter cluster position
    threevec _tpos; // extrapolated track position near calorimeter cluster
    threevec _tdir; // extrapolated track position near calorimeter cluster
    Float_t _ttrk; // track time at intersection point
    static string const& leafnames() { 
      static const string leaves =
      string("dt/F:du/F:dv/F:ds/F:ep/F:uvchisq/F:tchisq/F:dtllr/F:epllr/F:") + // matching info
      string("eclust/F:tclust/F:") + // cluster information
      threevec::leafnames("cpos").c_str() + string(":") +
      threevec::leafnames("tpos").c_str() + string(":") +
      threevec::leafnames("tdir").c_str() + string(":") +
      string("ttrk/F");
      return leaves;
    }

    void reset() {
      _dt = _du = _dv = _ds = _ep = _uvChisq = _tChisq = _dtllr = _epllr = _eclust = _tclust = _ttrk -1.0;
      _cpos = _tpos = _tdir = threevec();
    }
  };
}
#endif
