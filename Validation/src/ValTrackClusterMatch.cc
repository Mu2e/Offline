
#include "art_root_io/TFileDirectory.h"
#include "Validation/inc/ValTrackClusterMatch.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "TMath.h"

int mu2e::ValTrackClusterMatch::declare(art::TFileDirectory tfs) {

  _hVer = tfs.make<TH1D>( "Ver", "Version Number", 101, -0.5, 100.0);

  _hNMatch = tfs.make<TH1D>( "NMatch", "N Matches", 11, -0.5, 10.5);
  _hdu = tfs.make<TH1D>( "du", "du", 100, -1000., 1000.);
  _hdv = tfs.make<TH1D>( "dv", "dv", 100, -1000., 1000.);
  _hdt = tfs.make<TH1D>( "dt", "dt", 100, -15., 15.);
  _hep = tfs.make<TH1D>( "ep", "E/p", 100, 0., 1.6);
  _hchi2 = tfs.make<TH1D>( "Chi2", "Chi2", 100, 0.0, 10.0);
  _hchi2t = tfs.make<TH1D>( "Chi2t", "Chi2 time", 100, 0.0, 10.0);
  _hchi2t2 = tfs.make<TH1D>( "Chi2t2", "Chi2 time", 100, 0.0, 100.0);

  return 0;
}

int mu2e::ValTrackClusterMatch::fill(
                       const mu2e::TrackClusterMatchCollection & coll,
	               art::Event const& event) {

  // increment this by 1 any time the defnitions of the histograms or the 
  // histogram contents change, and will not match previous versions
  _hVer->Fill(1.0);

  _hNMatch->Fill(coll.size());
  for(auto match : coll) {

    _hdu->Fill(match.du());
    _hdv->Fill(match.dv());
    _hdv->Fill(match.dv());
    _hdt->Fill(match.dt());
    _hep->Fill(match.ep());
    _hchi2->Fill(match.chi2());
    _hchi2t->Fill(match.chi2_time());
    _hchi2t2->Fill(match.chi2_time());

  }

  return 0;
}

