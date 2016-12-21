
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "Validation/inc/ValPosition.hh"

int mu2e::ValPosition::declare(art::TFileDirectory tfs, 
			       std::string name, std::string title) {
  _hx = tfs.make<TH1D>( (name+"x").c_str(), 
			(title+"X").c_str(), 
			100, -6000., 6000.);
  _hy = tfs.make<TH1D>( (name+"y").c_str(), 
			(title+"Y").c_str(), 
			100, -2000., 2000.);
  _hz = tfs.make<TH1D>( (name+"z").c_str(), 
			(title+"Z").c_str(), 
			100, -20000., 20000.);
  _hxf = tfs.make<TH1D>( (name+"xfold").c_str(), 
			 (title+"X fold").c_str(), 
			 100, 0.0, 1000.0);
  _hyf = tfs.make<TH1D>( (name+"yfold").c_str(), 
			 (title+"Y fold").c_str(), 
			 100, 0.0, 1000.0);
  _hzf = tfs.make<TH1D>( (name+"zfold").c_str(), 
			 (title+"Z fold").c_str(), 
			 100, 0.0, 1000.0);
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
