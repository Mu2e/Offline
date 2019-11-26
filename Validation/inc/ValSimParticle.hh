
#ifndef ValSimParticle_HH_
#define ValSimParticle_HH_

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "Validation/inc/ValId.hh"
#include "TH1D.h"
#include <string>

namespace mu2e {

  class ValSimParticle {

  public:
    ValSimParticle(std::string name):_name(name){}
    int declare( art::TFileDirectory tfs);
    int fill(const SimParticleCollection & coll, art::Event const& event);
    std::string& name() { return _name; }

  private:

    std::string _name;
    
    TH1D* _hVer;
    TH1D* _hN;
    TH1D* _hN2;
    ValId _id;
    TH1D* _hp;
    TH1D* _hendKE;
    TH1D* _hpe;
    TH1D* _hpe2;
    TH1D* _hpg;
    TH1D* _hpg2;
    TH1D* _hpm;
    TH1D* _hp0;
    TH1D* _hpi;
    TH1D* _hpk0;
    TH1D* _hpk;
    TH1D* _hpn0;
    TH1D* _hpn02;
    TH1D* _hpn;
    TH1D* _hpn2;
    TH1D* _hsx;
    TH1D* _hsy;
    TH1D* _hsz;
    TH1D* _hsxDS;
    TH1D* _hsyDS;
    TH1D* _hszDS;
    TH1D* _hex;
    TH1D* _hey;
    TH1D* _hez;
    TH1D* _hexDS;
    TH1D* _heyDS;
    TH1D* _hezDS;
    TH1D* _hscode;
    TH1D* _hecode;
    TH1D* _hNS;
    TH1D* _hNS2;
    TH1D* _hl1;
    TH1D* _hl2;
    TH1D* _hl3;
    TH1D* _hND;
    TH1D* _hND2;
    ValId _idh;
    TH1D* _hscodeh;
    TH1D* _hecodeh;
    ValId _idhendKE;
    TH1D* _hscodehendKE;
    TH1D* _hecodehendKE;
    ValId _idh9endKE;
    TH1D* _hscodeh9endKE;
    TH1D* _hecodeh9endKE;
    TH1D* _htx;
    TH1D* _hty;
    TH1D* _htz;
    TH1D* _tgtmux;
    TH1D* _tgtmuy;
    TH1D* _tgtmuz;


  };
}


#endif
