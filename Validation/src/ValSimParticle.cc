
#include "Validation/inc/ValSimParticle.hh"


int mu2e::ValSimParticle::declare(art::TFileDirectory tfs) {
  _hVer = tfs.make<TH1D>( "Ver", "Version Number", 101, -0.5, 100.0);
  _hN = tfs.make<TH1D>( "Nsim", "N particle", 100, -0.05, 1000.0);
  _hN2 = tfs.make<TH1D>( "Nsim2", "log10(N particle)", 100, 0.0, 6.00);
  _id.declare(tfs,"id","id fold");
  _hp = tfs.make<TH1D>( "p", "P", 100, 0.0, 200.0);
  _pos.declare(tfs);
  _hscode = tfs.make<TH1D>( "scode", "start code", 151, -0.5, 150.0);
  _hecode = tfs.make<TH1D>( "ecode", "end code", 151, -0.5, 150.0);
  _idh.declare(tfs,"idh","id fold, p>10");
  _hscodeh = tfs.make<TH1D>( "scodeh", "start code, p>10", 151, -0.5, 150.0);
  _hecodeh = tfs.make<TH1D>( "ecodeh", "end code, p>10", 151, -0.5, 150.0);
  return 0;
}

int mu2e::ValSimParticle::fill(const mu2e::SimParticleCollection & coll,
				art::Event const& event) {

  // increment this by 1 any time the defnitions of the histograms or the 
  // histogram contents change, and will not match previous versions
  _hVer->Fill(1.0);

  _hN->Fill(coll.size()); 
  double x = (coll.size()<=0 ? 0 : log10(coll.size()) );
  _hN2->Fill(x);
  for(auto sp : coll) {
    const mu2e::SimParticle& part = sp.second;
    double pstart = part.startMomentum().vect().mag();
    //int idc = _id.fill(part.pdgId());
    _id.fill(part.pdgId());
    _hp->Fill(part.startMomentum().vect().mag()); 
    _pos.fill(part.endPosition());
    _hscode->Fill(part.realCreationCode().id());
    _hecode->Fill(part.stoppingCode().id());

    if(pstart>10.0) {
      _idh.fill(part.pdgId()); 
      _hscodeh->Fill(part.realCreationCode().id());
      _hecodeh->Fill(part.stoppingCode().id());
    }
  }
  return 0;
}
