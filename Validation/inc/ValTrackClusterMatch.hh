
#ifndef ValTrackClusterMatch_HH_
#define ValTrackClusterMatch_HH_

#include "art/Framework/Principal/Event.h"
#include "RecoDataProducts/inc/TrackClusterMatch.hh"
#include "TH1D.h"
#include <string>

namespace mu2e {

  class ValTrackClusterMatch {

  public:
    ValTrackClusterMatch(std::string name):_name(name){}
    int declare( art::TFileDirectory tfs);
    int fill(const TrackClusterMatchCollection & coll, art::Event const& event);
    double mcTrkP(art::Event const& event);
    std::string& name() { return _name; }

  private:
    std::string _name;

    TH1D* _hVer;
    TH1D* _hNMatch;
    TH1D* _hdu;
    TH1D* _hdv;
    TH1D* _hdt;
    TH1D* _hep;
    TH1D* _hchi2;
    TH1D* _hchi2t;
    TH1D* _hchi2t2;
  };
}


#endif
