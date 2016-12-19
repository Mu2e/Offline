
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "Validation/inc/ValPosition.hh"

int mu2e::ValPosition::declare(art::TFileDirectory tfs) {
  _hx = tfs.make<TH1D>( "x", "X", 100, -6000., 6000.);
  _hy = tfs.make<TH1D>( "y", "Y", 100, -2000., 2000.);
  _hz = tfs.make<TH1D>( "z", "Z", 100, -20000., 20000.);
  _hxf = tfs.make<TH1D>( "xfold", "X fold", 100, 0.0, 1000.0);
  _hyf = tfs.make<TH1D>( "yfold", "Y fold", 100, 0.0, 1000.0);
  _hzf = tfs.make<TH1D>( "zfold", "Z fold", 100, 0.0, 1000.0);
  return 0;
}

int mu2e::ValPosition::fill(CLHEP::Hep3Vector const& position) {
  _hx->Fill(position.x()); 
  _hy->Fill(position.y()); 
  _hz->Fill(position.z()); 
  _hxf->Fill( fold(position.x()) ); 
  _hyf->Fill( fold(position.y()) ); 
  _hzf->Fill( fold(position.z()) ); 
  return 0;
}
