//
// StrawWaveform represents post-amplification voltage as a function of time at one end of a
// a straw, over the time period of 1 microbunch.  It includes all physical and electronics
// effects prior to digitization.
//
// $Id: StrawWaveform.cc,v 1.4 2013/12/10 21:43:45 brownd Exp $
// $Author: brownd $
// $Date: 2013/12/10 21:43:45 $
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

  bool StrawWaveform::crossesThreshold(double threshold,WFX& wfx) const {
    bool retval(false);
    static double epsilon(1.0e-3); // small time to avoid problems with 0 rise times
    static double timestep(0.010); // interpolation minimum to use linear threshold crossing calculation
  // loop forward from this hitlet
    while(wfx._ihitlet != _hseq.hitletList().end()){
  // sample the wavefom just before the referenced hitlet
      double pretime =wfx._ihitlet->time()-epsilon;
      double presample = sampleWaveform(pretime);
// if pre-hitlet response is below threshold, check for crossing
      if(presample < threshold){
// check response at maximum
	double posttime = wfx._ihitlet->time() + _strawele.maxResponseTime();
	double postsample = sampleWaveform(posttime);
	if(postsample > threshold){
// this hitlet pushes the waveform over threshold
	  retval = true;
// compute the actual crossing time, using linear interpolation
	  double dt = posttime-pretime;
	  do {
	    dt *= 0.5;
	    double t = pretime + dt;
	    double s = sampleWaveform(t);
	    if(s>threshold){
	      posttime = t;
	      postsample = s;
	    } else {
	      pretime = t;
	      presample = s;
	    }
// linear interpolation
	    double slope = dt/(postsample-presample);
	    wfx._time = pretime + slope*(threshold-presample);
	  } while(dt > timestep);
// smear the time for electronics noise: FIXME!!
	  break;
	}
      }
  // advance to next hitlet
      ++(wfx._ihitlet);
    }
    return retval;
  }

  double StrawWaveform::sampleWaveform(double time) const {
// loop over all hitlets and add their response at this time
    HitletList const& hlist = _hseq.hitletList();
    double retval(0.0);
    auto ihitlet = hlist.begin();
    while(ihitlet != hlist.end() && ihitlet->time() < time){
    // cmopute the straw electronics response to this charge.  This is pre-saturation 
      retval += _strawele.hitletResponse(time,*ihitlet);
    // move to next hitlet
      ++ihitlet;
    }
    // apply saturation effects
    retval = _strawele.saturatedResponse(retval);
    return retval;
  }
  
  void StrawWaveform::sampleWaveform(std::vector<double> const& times,std::vector<double>& volts) const {
    volts.clear();
    volts.reserve(times.size());
    for(auto itime=times.begin();itime!=times.end();++itime){
      volts.push_back(sampleWaveform(*itime));
    }
  }
  
}

