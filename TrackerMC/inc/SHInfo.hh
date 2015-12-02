// diagnostic structs
#include "Rtypes.h"
#include "DataProducts/inc/StrawId.hh"
namespace mu2e {

  struct SHID {
    Int_t _device, _panel, _layer, _straw;
    SHID() : _device(-1), _panel(-1), _layer(-1), _straw(-1) {}
    SHID(StrawId const& strawid) : 
	  _device(strawid.getDevice()),
	  _panel(strawid.getPanel()),
	  _layer(strawid.getLayer()),
	  _straw(strawid.getStraw()) {}
  };

  struct SHMCInfo {
    Float_t _energy,_trigenergy, _threshenergy;
    Int_t _pdg, _proc, _gen, _nmcpart;
    Float_t _mom, _dperp, _len;
    Int_t _ambig;
    Bool_t _xtalk;
  };
}

