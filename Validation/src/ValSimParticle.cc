
#include "Validation/inc/ValSimParticle.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"
#include "Mu2eUtilities/inc/PhysicalVolumeMultiHelper.hh"
#include <vector>

int mu2e::ValSimParticle::declare(art::TFileDirectory tfs) {
  _hVer = tfs.make<TH1D>( "Ver", "Version Number", 101, -0.5, 100.0);
  _hN = tfs.make<TH1D>( "Nsim", "N particle", 100, -0.05, 1000.0);
  _hN2 = tfs.make<TH1D>( "Nsim2", "log10(N particle)", 100, 0.0, 6.00);
  _id.declare(tfs,"id","id fold");
  _hp = tfs.make<TH1D>( "p", "P", 100, 0.0, 200.0);
  _hpe = tfs.make<TH1D>( "pe", "P ele", 100, 0.0, 200.0);
  _hpm = tfs.make<TH1D>( "pm", "P muon", 100, 0.0, 600.0);
  _hp0 = tfs.make<TH1D>( "p0", "P pi0", 100, 0.0, 600.0);
  _hpi = tfs.make<TH1D>( "pi", "P pi+/-", 100, 0.0, 600.0);
  _hpn = tfs.make<TH1D>( "pn", "P nuclei", 100, 0.0, 200.0);
  _spos.declare(tfs,"start","start");
  _epos.declare(tfs,"end","end");
  _hscode = tfs.make<TH1D>( "scode", "start code", 151, -0.5, 150.0);
  _hecode = tfs.make<TH1D>( "ecode", "end code", 151, -0.5, 150.0);
  _idh.declare(tfs,"idh","id fold, p>10");
  _hscodeh = tfs.make<TH1D>( "scodeh", "start code, p>10", 151, -0.5, 150.0);
  _hecodeh = tfs.make<TH1D>( "ecodeh", "end code, p>10", 151, -0.5, 150.0);
  _stopPos.declare(tfs,"stop","stopped particle");
  _tgtmux = tfs.make<TH1D>( "tgtmux", "target stopped mu x", 50, -4000.0, -3800.0);
  _tgtmuy = tfs.make<TH1D>( "tgtmuy", "target stopped mu y", 50, -100.0, 100.0);
  _tgtmuz = tfs.make<TH1D>( "tgtmuz", "target stopped mu z", 100, -6350.0, -5400.0);

  return 0;
}

int mu2e::ValSimParticle::fill(const mu2e::SimParticleCollection & coll,
				art::Event const& event) {

  // increment this by 1 any time the defnitions of the histograms or the 
  // histogram contents change, and will not match previous versions
  _hVer->Fill(3.0);

  _hN->Fill(coll.size()); 
  double x = (coll.size()<=0 ? 0 : log10(coll.size()) );
  _hN2->Fill(x);

  // find geometry volumes, usually saved from production
  art::SubRun const& subrun = event.getSubRun();
  std::vector<art::Handle< PhysicalVolumeInfoMultiCollection >> vah;
  subrun.getManyByType(vah);
  PhysicalVolumeMultiHelper *vhelper = nullptr;
  if(vah.size()>0) vhelper = new PhysicalVolumeMultiHelper(*(vah[0]));


  for(auto sp : coll) {
    const mu2e::SimParticle& part = sp.second;
    double pstart = part.startMomentum().vect().mag();
    //int idc = _id.fill(part.pdgId());
    int idc =_id.fill(part.pdgId());
    double p = part.startMomentum().vect().mag();
    _hp->Fill(p);
    if(abs(idc)==11) _hpe->Fill(p);
    if(abs(idc)==13) _hpm->Fill(p);
    if(abs(idc)==30) _hp0->Fill(p);
    if(abs(idc)==31) _hpi->Fill(p);
    if(abs(idc)==51) _hpn->Fill(p);
    _spos.fill(part.startPosition());
    _epos.fill(part.endPosition());
    _hscode->Fill(part.originParticle().creationCode().id());
    _hecode->Fill(part.stoppingCode().id());
    
    if(pstart>10.0) {
      _idh.fill(part.pdgId()); 
      _hscodeh->Fill(part.originParticle().creationCode().id());
      _hecodeh->Fill(part.stoppingCode().id());
    }

    // stopped muons
    if ( part.endMomentum().v().mag2() <= 
	      std::numeric_limits<double>::epsilon() &&
	      part.pdgId() == 13 && vhelper != nullptr) {  // mu-
      _stopPos.fill(part.endPosition());
      const std::string material = vhelper->endVolume(part).materialName();
      if (material.find("StoppingTarget_Al")!=std::string::npos) {
	_tgtmux->Fill( part.endPosition().x() );
	_tgtmuy->Fill( part.endPosition().y() );
	_tgtmuz->Fill( part.endPosition().z() );
      }
    }

  }
  return 0;
}
