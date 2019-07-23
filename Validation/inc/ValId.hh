#ifndef ValId_HH_
#define ValId_HH_

//
// A helper class to create, hold, and fill particle ID
// histograms to avoid copying this code several places
//
#include "art_root_io/TFileDirectory.h"
#include "TH1D.h"
#include <string>

namespace mu2e {

  class ValId {

  public:
    int declare( art::TFileDirectory tfs, 
		 std::string name="id", std::string title="id fold");
    int fill(int id);
    int compress(int id);
  private:
    
    TH1D* _hid;
  };
}


#endif
