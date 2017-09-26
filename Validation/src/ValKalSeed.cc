#include <cmath>
#include "Validation/inc/ValKalSeed.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TrackerGeom/inc/Tracker.hh"

int mu2e::ValKalSeed::declare(art::TFileDirectory tfs) {
  _hVer = tfs.make<TH1D>( "Ver", "Version Number", 101, -0.5, 100.0);
  _hN = tfs.make<TH1D>( "NSeed", "N KalSeed", 11, -0.5, 10.0);
  _hNStraw = tfs.make<TH1D>( "NHit", "N Straw", 101, -0.5, 100.0);
  _hNSeg = tfs.make<TH1D>( "NSeg", "N KalSegment", 21, -0.5, 20.0);
  _hStatus = tfs.make<TH1D>( "Status", "Status", 32, -0.5, 31.0);
  _hflt0 = tfs.make<TH1D>( "flt0", "flt0", 100, -1200.0, 1200.0);
  _ht0 = tfs.make<TH1D>( "t0", "t0", 100, 400.0, 1800.0);
  _hchi2 = tfs.make<TH1D>( "Chi2N", "Chi2/DOF", 100, 0.0, 100.0);
  _hhasCal = tfs.make<TH1D>( "hasCal", "CalCluster attached", 2, -0.5, 1.5);
  _hfitCon = tfs.make<TH1D>( "FitConn", "Fit CL", 100, 0.0, 1.0);
  _hp = tfs.make<TH1D>( "p", "p", 100, 0., 110.);
  _hpce = tfs.make<TH1D>( "pce", "p CE", 100, 95.0, 110.);
  _hpe = tfs.make<TH1D>( "pe", "p error", 100, 0.0, 1.0);
  _hD0 = tfs.make<TH1D>( "d0", "d0", 100, -100., 100.);
  _hPhi0 = tfs.make<TH1D>( "phi0", "phi0", 100, -M_PI, M_PI);
  _hOmega = tfs.make<TH1D>( "omega", "omega", 100, -0.02, 0.02);
  _hZ0 = tfs.make<TH1D>( "Z0", "Z0", 100, -1000.0, 1000.0);
  _hTan = tfs.make<TH1D>( "tanDip", "tanDip", 100, -1.8, 1.8);

  return 0;
}

int mu2e::ValKalSeed::fill(const mu2e::KalSeedCollection & coll,
				art::Event const& event) {

  // increment this by 1 any time the defnitions of the histograms or the 
  // histogram contents change, and will not match previous versions
  _hVer->Fill(0.0);
  _hN->Fill(coll.size()); 
  for(auto const& ks : coll) {
    _hNStraw->Fill(ks.straws().size());
    _hNSeg->Fill(ks.segments().size());
    const TrkFitFlag& f = ks.status();

    //    for(mu2e::TrkFitFlagDetail::bit_type i=0; i<f.size(); i++) 
    //  if(f.hasAnyProperty(i)) _hStatus->Fill(i); 

    int i=0;
    for(auto sn: f.bitNames()) { 
      if(f.hasAnyProperty(TrkFitFlag(sn.first))) _hStatus->Fill(i); 
      i++;
    }

    _hflt0->Fill(ks.flt0());
    _ht0->Fill(ks.t0().t0());
    _hchi2->Fill(ks.chisquared());
    int q = (ks.caloCluster().isNull()?0:1);
    _hhasCal->Fill(q);
    _hfitCon->Fill(ks.fitConsistency());

    if( ks.segments().size()>0 ) {
      std::size_t i = ks.segments().size()/2;
      auto const& ss = ks.segments()[i]; //KalSegment
      auto const& h = ss.helix(); // HelixVal
      _hp->Fill(ss.mom());
      _hpce->Fill(ss.mom());
      _hpe->Fill(ss.momerr());
      _hD0->Fill(h.d0());
      _hPhi0->Fill(h.phi0());
      _hOmega->Fill(h.omega());
      _hZ0->Fill(h.z0());
      _hTan->Fill(h.tanDip());
    }
    
  }
  return 0;
}
