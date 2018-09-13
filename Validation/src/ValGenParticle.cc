
#include "Validation/inc/ValGenParticle.hh"


int mu2e::ValGenParticle::declare(art::TFileDirectory tfs) {
  _hVer = tfs.make<TH1D>( "Ver", "Version Number", 101, -0.5, 100.0);
  _hN = tfs.make<TH1D>( "Ngen", "N Particle", 100, -0.05, 1000.0);
  _id.declare(tfs);
  _hp = tfs.make<TH1D>( "p", "P", 100, 0.0, 200.0);
  _hlogp = tfs.make<TH1D>( "logp", "log10 P", 100, 0.0, 13.0);
  _hx = tfs.make<TH1D>( "X", "X", 100, -14000.0, 6000.0);
  _hxt = tfs.make<TH1D>( "Xtgt", "X stopping target", 100, -3904-120.0, -3904+120.0);
  _hy = tfs.make<TH1D>( "Y", "Y", 100, -2000.0, 2000.0);
  _hyt = tfs.make<TH1D>( "Ytgt", "Y stopping target", 100, -120.0, 120.0);
  _hz = tfs.make<TH1D>( "Z", "Z", 100, -20000.0, 25000.0);
  _hzt = tfs.make<TH1D>( "Ztgt", "Z stopping target", 100, 5400, 6350);
  return 0;
}

int mu2e::ValGenParticle::fill(const mu2e::GenParticleCollection & coll,
				art::Event const& event) {

  // increment this by 1 any time the defnitions of the histograms or the 
  // histogram contents change, and will not match previous versions
  _hVer->Fill(1.0);

  _hN->Fill(coll.size()); 
  for(auto sp : coll) {
    _id.fill(sp.pdgId()); 
    double p = sp.momentum().vect().mag();
    _hp->Fill(p); 
    _hlogp->Fill( (p>0.0? log10(p) : -1.0) );
    _hx->Fill(sp.position().x());
    _hxt->Fill(sp.position().x());
    _hy->Fill(sp.position().y());
    _hyt->Fill(sp.position().y());
    _hz->Fill(sp.position().z());
    _hzt->Fill(sp.position().z());
    
  }
  return 0;
}
