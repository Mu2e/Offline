
#include "Validation/inc/ValCaloRecoDigi.hh"


int mu2e::ValCaloRecoDigi::declare(art::TFileDirectory tfs) {
  _hVer = tfs.make<TH1D>( "Ver", "Version Number", 101, -0.5, 100.0);
  _hN = tfs.make<TH1D>( "NDigis", "N RecoDigis", 101, -0.5, 100.5);
  _hN2= tfs.make<TH1D>( "NDigis2", "N RecoDigis", 100, -0.5,4999.5);
  _hI = tfs.make<TH1D>( "ROID", "ROID",100, 0.0, 3000.0);
  _ht = tfs.make<TH1D>( "t", "time", 100, 0.0, 2000.0);
  _hE = tfs.make<TH1D>( "E", "Energy",100, 0.0, 150.0);
  _hc = tfs.make<TH1D>( "chi2", "chi2",100, 0.0, 10.0);
  _hp = tfs.make<TH1D>( "pileup", "pileup flag",2, -0.5, 1.5);

  return 0;
}

int mu2e::ValCaloRecoDigi::fill(const mu2e::CaloRecoDigiCollection & coll,
				art::Event const& event) {

  // increment this by 1 any time the defnitions of the histograms or the 
  // histogram contents change, and will not match previous versions
  _hVer->Fill(0.0);

  _hN->Fill(coll.size()); 
  _hN2->Fill(coll.size()); 
  for(auto rd : coll) {
    _hI->Fill(rd.ROid());
    _ht->Fill(rd.time());
    _hE->Fill(rd.energyDep());
    _hc->Fill(rd.chi2());
    _hp->Fill(rd.pileUp());
  }
  return 0;
}
