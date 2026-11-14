#ifndef ValTrkCrvMatch_HH_
#define ValTrkCrvMatch_HH_

#include "Offline/RecoDataProducts/inc/TrkCrvMatch.hh"
#include "Offline/Validation/inc/ValTrkCrvSectorAxis.hh"
#include "KinKal/General/Vectors.hh"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileDirectory.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include <cstdint>
#include <map>
#include <string>

namespace mu2e {

// track-CRV association: candidate counts, residuals and chi2 per CRV sector, dv against slope (alignment monitor),
// predicted vs matched sector, layer probe
class ValTrkCrvMatch {
 public:
  ValTrkCrvMatch(std::string name) : _name(name) {}
  int declare(const art::TFileDirectory& tfs);
  int fill(const TrkCrvMatchCollection& coll, art::Event const& event);
  std::string& name() { return _name; }

 private:
  std::string _name;
  struct Plane {                      // a sector of this job's geometry: centre, across-bar direction, normal, dv vs slope
    KinKal::VEC3 c, v, w;
    TProfile* p = nullptr;
  };
  ValTrkCrvSectorAxis _axis;
  std::map<int, Plane> _plane;        // by axis bin
  std::map<int16_t, int> _shieldBin;  // product shield index (sector_) -> axis bin, from the crossings
  static constexpr double _planeTol = 1.0;  // max distance of a crossing from the job's sector mid-plane [mm]

  TH1D* _hVer;
  TH1D* _hNTrk;
  TH1D* _hNCand;
  TH1D* _hNGeom;
  TH1D* _hNXing;
  TH1D* _hNXingProj;
  TH1D* _hDtRef;
  TH1D* _hDtRefBest;
  TH1D* _hDtRefToF;
  TH1D* _hChi2;
  TH1D* _hChi2Proj;
  TH1D* _hNShare;
  TH1D* _hFlags;
  TH1D* _hXingShield;
  TH1D* _hXingGeomShield;
  TH1D* _hXingProjShield;
  TH2D* _hDtXing;
  TH2D* _hDu;
  TH2D* _hDv;
  TH2D* _hPredMatched;
  TH1D* _hSlopeNoGeom;
  TH1D* _hProbePred;
  TH1D* _hProbeHit;
  TH1D* _hProbeInCluster;
  TH1D* _hProbeDBar;
  TH1D* _hProbeDt;
  TH1D* _hProbePE;
  TH1D* _hProbeDu;
};
}  // namespace mu2e

#endif
