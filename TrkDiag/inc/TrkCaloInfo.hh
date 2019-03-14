//
// structs used to record calorimeter information matched to a track
// All energies are in units of MeV,  momenta are in units of MeV/c,
// time in nsec WRT when the proton bunch pulse peak hits the production target,
// positions are in mm WRT the center of the tracker.
 
// Dave Brown (LBNL)
// 
#ifndef TrkCaloInfo_HH
#define TrkCaloInfo_HH
#include "DataProducts/inc/XYZVec.hh"
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
    Int_t _section; // cluster section
    XYZVec _cpos; // calorimeter cluster position
    XYZVec _tpos; // extrapolated track position near calorimeter cluster
    XYZVec _tdir; // extrapolated track position near calorimeter cluster
    Float_t _ttrk; // track time at intersection point
    static string leafnames() { 
      static string leaves;
      leaves = 
      string("dt/F:du/F:dv/F:ds/F:ep/F:uvchisq/F:tchisq/F:dtllr/F:epllr/F:") + // matching info
      string("eclust/F:tclust/F:section/I:") + // cluster information
      Geom::XYZnames("cpos").c_str() + string(":") +
      Geom::XYZnames("tpos").c_str() + string(":") +
      Geom::XYZnames("tdir").c_str() + string(":") +
      string("ttrk/F");
      return leaves;
    }

    void reset() {
      _dt = _du = _dv = _ds = _ep = _uvChisq = _tChisq = _dtllr = _epllr = _eclust = _tclust = _ttrk = 0.0;
      _cpos = _tpos = _tdir = XYZVec();
    }
  };
}
#endif
