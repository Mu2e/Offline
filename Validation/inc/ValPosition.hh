#ifndef ValPosition_HH_
#define ValPosition_HH_

//
// A helper class to create, hold, and fill position
// histograms to avoid copying this code several places
//

#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "TH1D.h"
#include <string>

namespace mu2e {

  class ValPosition {

  public:
    int declare( art::TFileDirectory tfs, 
		 std::string name="", std::string title="");
    int fill(CLHEP::Hep3Vector const& position);
    double fold(double x) { return x - 1000.0*floor(x/1000.0); }

  private:
    
    TH1D* _hx;
    TH1D* _hy;
    TH1D* _hz;
    TH1D* _hxf;
    TH1D* _hyf;
    TH1D* _hzf;
  };
}


#endif
