#ifndef ValTrkCrvMatch_HH_
#define ValTrkCrvMatch_HH_

#include "Offline/RecoDataProducts/inc/TrkCrvMatch.hh"
#include "KinKal/General/Vectors.hh"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileDirectory.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include <string>
#include <vector>

namespace mu2e {

// track-CRV association: candidate counts, residuals and chi2 per CRV shield, dv against slope (alignment monitor),
// predicted vs matched shield, layer probe
class ValTrkCrvMatch {
 public:
  ValTrkCrvMatch(std::string name) : _name(name) {}
  int declare(const art::TFileDirectory& tfs);
  int fill(const TrkCrvMatchCollection& coll, art::Event const& event);
  std::string& name() { return _name; }

 private:
  std::string _name;
  std::vector<KinKal::VEC3> _v, _w;   // per shield: across-bar direction and normal, detector frame

  TH1D* _hVer;
  TH1D* _hNTrk;
  TH1D* _hNCand;
  TH1D* _hNGeom;
  TH1D* _hNXing;
  TH1D* _hDtRef;
  TH1D* _hDtRefBest;
  TH1D* _hDtRefToF;
  TH1D* _hChi2;
  TH1D* _hNShare;
  TH1D* _hFlags;
  TH1D* _hXingShield;
  TH1D* _hXingGeomShield;
  TH2D* _hDtXing;
  TH2D* _hDu;
  TH2D* _hDv;
  TH2D* _hPredMatched;
  std::vector<TProfile*> _pDvSlope;
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
