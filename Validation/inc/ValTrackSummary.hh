
#ifndef ValTrackSummary_HH_
#define ValTrackSummary_HH_

#include "art/Framework/Principal/Event.h"
#include "RecoDataProducts/inc/TrackSummary.hh"
#include "TH1D.h"
#include <string>

namespace mu2e {

  class ValTrackSummary {

  public:
    ValTrackSummary(std::string name):_name(name){}
    int declare( art::TFileDirectory tfs);
    int fill(const TrackSummaryCollection & coll, art::Event const& event);
    double mcTrkP(art::Event const& event);
    std::string& name() { return _name; }

  private:
    std::string _name;

    TH1D* _hVer;
    TH1D* _hNTr;
    TH1D* _hNState;
    TH1D* _hp;
    TH1D* _hpce;
    TH1D* _hD0;
    TH1D* _hPhi0;
    TH1D* _hOmega;
    TH1D* _hZ0;
    TH1D* _hTan;
    TH1D* _hT0;
    TH1D* _hNActive;
    TH1D* _hNDof;
    TH1D* _hChi2N;
    TH1D* _hCL;
    TH1D* _hStatus;
    TH1D* _hCuts;
    TH1D* _hPRes;
  };
}


#endif
