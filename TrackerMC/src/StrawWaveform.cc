//
// StrawWaveform represents post-amplification voltage as a function of time at one end of a
// a straw, over the time period of 1 microbunch.  It includes all physical and electronics
// effects prior to digitization.
//
// $Id: StrawWaveform.cc,v 1.8 2014/01/18 17:33:34 brownd Exp $
// $Author: brownd $
// $Date: 2014/01/18 17:33:34 $
//
// Original author David Brown, LBNL
//
#include "TrackerMC/inc/StrawWaveform.hh"
#include <math.h>

using namespace std;
namespace mu2e {

  StrawWaveform::StrawWaveform(StrawHitletSequence const& hseq, ConditionsHandle<StrawElectronics> const& strawele) :
   _hseq(hseq), _strawele(strawele) 
  {}

  bool StrawWaveform::crossesThreshold(double threshold,WFX& wfx) const {
    bool retval(false);
    static double timestep(0.010); // interpolation minimum to use linear threshold crossing calculation
  // loop forward from this hitlet
    while(wfx._ihitlet != _hseq.hitletList().end() ){
  // require the hitlet be beyond the input time
      if(wfx._ihitlet->time() > wfx._time){
	// sample the wavefom just before the referenced hitlet
	double pretime =wfx._ihitlet->time();
	double presample = sampleWaveform(pretime);
	// if pre-hitlet response is below threshold, check for crossing
	if(presample < threshold){
	  // check response at maximum
	  double posttime = pretime + _strawele->maxResponseTime();
	  double postsample = sampleWaveform(posttime);
	  if(postsample > threshold){
	    // this hitlet pushes the waveform over threshold
	    retval = true;
	    // compute the actual crossing time, using linear interpolation
	    double dt = posttime-pretime;
	    while(dt > timestep) {
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
	    }
	    // linear interpolation
	    double slope = dt/(postsample-presample);
	    wfx._time = pretime + slope*(threshold-presample);
	    break;
	  }
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
    // cmopute the linear straw electronics response to this charge.  This is pre-saturation 
      retval += _strawele->linearResponse(time-ihitlet->time(),ihitlet->charge());
    // move to next hitlet
      ++ihitlet;
    }
    // apply saturation effects
    retval = _strawele->saturatedResponse(retval);
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

