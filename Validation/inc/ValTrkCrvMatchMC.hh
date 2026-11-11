#ifndef ValTrkCrvMatchMC_HH_
#define ValTrkCrvMatchMC_HH_

#include "Offline/MCDataProducts/inc/TrkCrvMatchMC.hh"
#include "Offline/Validation/inc/ValTrkCrvSectorAxis.hh"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileDirectory.h"
#include "TH1D.h"
#include "TH2D.h"
#include <string>

namespace mu2e {

// track-CRV association truth: relation of the best candidates, own-cluster association and crossing truth per
// CRV sector, reco crossings against the true crossings
class ValTrkCrvMatchMC {
 public:
  ValTrkCrvMatchMC(std::string name) : _name(name) {}
  int declare(const art::TFileDirectory& tfs);
  int fill(const TrkCrvMatchMCCollection& coll, art::Event const& event);
  std::string& name() { return _name; }

 private:
  std::string _name;
  ValTrkCrvSectorAxis _axis;

  TH1D* _hVer;
  TH1D* _hBestGeomRel;
  TH1D* _hBestTimeRel;
  TH1D* _hGeomCandRel;
  TH1D* _hTimeCandRel;
  TH1D* _hTrueClu;
  TH1D* _hTrueCluGeom;
  TH1D* _hTrueCluBest;
  TH1D* _hRecoXing;
  TH1D* _hRecoXingTrue;
  TH1D* _hTrueXing;
  TH1D* _hTrueXingFound;
  TH2D* _hXingDist;
  TH2D* _hXingDt;
};
}  // namespace mu2e

#endif
