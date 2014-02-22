//
// StrawWaveform represents post-amplification voltage as a function of time at one end of a
// a straw, over the time period of 1 microbunch.  It includes all physical and electronics
// effects prior to digitization.
//
// $Id: StrawWaveform.cc,v 1.9 2014/02/22 02:17:09 brownd Exp $
// $Author: brownd $
// $Date: 2014/02/22 02:17:09 $
//
// Original author David Brown, LBNL
//
#include "TrackerMC/inc/StrawWaveform.hh"
#include <math.h>
#include <boost/math/special_functions/binomial.hpp>

using namespace std;
namespace mu2e {

  StrawWaveform::StrawWaveform(StrawHitletSequence const& hseq, ConditionsHandle<StrawElectronics> const& strawele) :
   _hseq(hseq), _strawele(strawele) 
  {}

  bool StrawWaveform::crossesThreshold(double threshold,WFX& wfx) const {
    bool retval(false);
    // advance to the hitlet just before the time specified.  Buffer should be a parameter FIXME!!!
    while(wfx._ihitlet != _hseq.hitletList().end() && 
      wfx._ihitlet->time() + 1.0 < wfx._time){
	++(wfx._ihitlet);
    }
    if(wfx._ihitlet != _hseq.hitletList().end()){
      // sample initial voltage for this hitlet
      wfx._vstart = sampleWaveform(wfx._ihitlet->time());
      // if we start above threhold, scan forward till we're below
      if(wfx._vstart > threshold)
	returnCrossing(threshold, wfx);
      // scan ahead quickly from there to where there's enough integral charge to possibly cross threshold.
      if(roughCrossing(threshold,wfx)) {
	// start the fine scan from this hitlet.  
	while(wfx._ihitlet != _hseq.hitletList().end() ){
	  // First, get the starting voltage
	  wfx._vstart = sampleWaveform(wfx._ihitlet->time());
	  // check if this hitlet could cross threshold
	  if(wfx._vstart + _strawele->maxLinearResponse(wfx._ihitlet->charge()) > threshold){
	    // check the actual response 
	    double maxtime = wfx._ihitlet->time() + _strawele->maxResponseTime();
	    double maxresp = sampleWaveform(maxtime);
	    if(maxresp > threshold){
	      retval = true;
	      // interpolate to find the precise crossing
	      fineCrossing(threshold,maxresp,wfx);
	      break;
	    }
	  }
	  // advance to next hitlet
	  ++(wfx._ihitlet);
	}
      }
    }
    return retval;
  }

  void StrawWaveform::returnCrossing(double threshold, WFX& wfx) const {
    while(wfx._ihitlet != _hseq.hitletList().end() && wfx._vstart > threshold) {
      // move forward in time at least as twice the time to the maxium for this hitlet
      double time = wfx._ihitlet->time() + 2*_strawele->maxResponseTime(); 
      while(wfx._ihitlet != _hseq.hitletList().end() && 
	  wfx._ihitlet->time() < time){
	++(wfx._ihitlet);
      }
      if(wfx._ihitlet != _hseq.hitletList().end())
	wfx._vstart = sampleWaveform(wfx._ihitlet->time());
	wfx._time =wfx._ihitlet->time();
    }
  }

  bool StrawWaveform::roughCrossing(double threshold, WFX& wfx) const {
    // add voltage till we go over threshold by simple linear sum.  That's a minimum requirement
    // for actually crossing threshold
    double resp = wfx._vstart;
    while(wfx._ihitlet != _hseq.hitletList().end()){
      resp += _strawele->maxLinearResponse(wfx._ihitlet->charge());
      if(resp > threshold)break;
      ++(wfx._ihitlet);
    }
    // update time
    if(wfx._ihitlet != _hseq.hitletList().end() )
      wfx._time = wfx._ihitlet->time();

    return wfx._ihitlet != _hseq.hitletList().end() && resp > threshold;
  }

  bool StrawWaveform::fineCrossing(double threshold,double maxresp, WFX& wfx) const {
    static double timestep(0.010); // interpolation minimum to use linear threshold crossing calculation
    double pretime = wfx._ihitlet->time();
    double posttime = pretime + _strawele->maxResponseTime();
    double presample = wfx._vstart;
    double postsample = maxresp;
    static const unsigned maxstep(10); // 10 steps max
    unsigned nstep(0);
    double deltat = posttime-pretime;
    double dt = deltat;
    double slope = deltat/(postsample-presample);
    double time = pretime + slope*(threshold-presample);
    // linear interpolation
    while(dt > timestep && nstep < maxstep) {
      double sample = sampleWaveform(time);
      if(sample > threshold){
	posttime = time;
	postsample = sample;
      } else {
	pretime = time;
	presample = sample;
      }
      deltat = posttime-pretime;
      slope = deltat/(postsample-presample);
      double oldtime = time;
      time = pretime + slope*(threshold-presample);
      dt = fabs(time-oldtime);
      ++nstep;
    }     
    // set crossing time
    wfx._time = time;
    // update the referenced hitlet: this can be different than the one we started with!
    while(wfx._ihitlet != _hseq.hitletList().end() && 
      wfx._ihitlet->time() < wfx._time){
      ++(wfx._ihitlet);
    }
    // back off one
    if(wfx._ihitlet != _hseq.hitletList().begin())--(wfx._ihitlet);
    // return on convergence
    return dt < timestep;
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

