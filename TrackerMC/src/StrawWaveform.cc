//
// StrawWaveform represents post-amplification voltage as a function of time at one end of a
// a straw, over the time period of 1 microbunch.  It includes all physical and electronics
// effects prior to digitization.
//
// $Id: StrawWaveform.cc,v 1.1 2013/12/07 19:51:42 brownd Exp $
// $Author: brownd $
// $Date: 2013/12/07 19:51:42 $
//
// Original author David Brown, LBNL
//
#include "TrackerMC/inc/StrawWaveform.hh"
#include <math.h>

using namespace std;
namespace mu2e {

  StrawWaveform::StrawWaveform(StrawHitletSequence const& hseq, StrawElectronics const& strawele) :
   _hseq(hseq), _strawele(strawele) 
  {}

  bool StrawWaveform::crossesThreshold(double threshold,crossdir cdir,WFX& wfx) const {
   // find the maximum fr
    // If the current hitlet is already about threshold, move to the next
    // scan forward from the currently referenced hitlet to the next crossing point
    // needs implementation, FIXME!!
    return false;
  }

  double StrawWaveform::sampleWaveform(double time) const {
// loop over all hitlets and add their response at this time
// needs implementation, FIXME!!!
    return 0.0;
  }
  
  void StrawWaveform::sampleWaveform(std::vector<double> const& times,std::vector<double>& volts) const {
    volts.clear();
    volts.reserve(times.size());
    for(auto itime=times.begin();itime!=times.end();++itime){
      volts.push_back(sampleWaveform(*itime));
    }
  }
  
}

