//
// Original author Stefano Roberto Soleti
//

// C++ includes
#include <ostream>

// Framework includes.
#include "cetlib/exception.h"

// Mu2e includes
#include "RecoDataProducts/inc/CaloDigi.hh"

using namespace std;

namespace mu2e {

  // Print the information found in this hit. TODO
  void CaloDigi::print( ) const {

    printf(" CaloDigi output: roId = %i t0 = %5.2f nSamples = %i waveform = ", _roId, _t0, _nSamples);
    
    for (int i=0; i<int(_waveform.size()); ++i){
      printf("%i ", _waveform.at(i));
    }
    printf("\n");
  }

} // namespace mu2e
