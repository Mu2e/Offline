//
// StrawWaveform represents post-amplification voltage as a function of time at one end of a
// a straw, over the time period of 1 microbunch.  It includes all physical and electronics
// effects prior to digitization.
//
// $Id: StrawWaveform.cc,v 1.14 2014/03/04 00:29:17 brownd Exp $
// $Author: brownd $
// $Date: 2014/03/04 00:29:17 $
//
// Original author David Brown, LBNL
//
#include "TrackerMC/inc/StrawWaveform.hh"
#include <math.h>
#include <boost/math/special_functions/binomial.hpp>

using namespace std;
namespace mu2e {
  namespace TrackerMC {
    StrawWaveform::StrawWaveform(StrawClusterSequence const& hseq, ConditionsHandle<StrawElectronics> const& strawele, XTalk const& xtalk) :
      _cseq(hseq), _strawele(strawele), _xtalk(xtalk)
    {}

    StrawWaveform::StrawWaveform(StrawWaveform const& other) : _cseq(other._cseq),
    _strawele(other._strawele), _xtalk(other._xtalk)
    {}

    bool StrawWaveform::crossesThreshold(double threshold,WFX& wfx) const {
      bool retval(false);
      // make sure we start past the input time
      while(wfx._iclust != _cseq.clustList().end() && wfx._iclust->time()< wfx._time ){
	++(wfx._iclust);
      }
      // loop till we're at the end or we go over threshold
      if(wfx._iclust != _cseq.clustList().end()){
	// sample initial voltage for this clust
	wfx._vstart = sampleWaveform(StrawElectronics::thresh,wfx._iclust->time());
	// if we start above threhold, scan forward till we're below
	if(wfx._vstart > threshold)
	  returnCrossing(threshold, wfx);
	// scan ahead quickly from there to where there's enough integral charge to possibly cross threshold.
	if(roughCrossing(threshold,wfx)) {
	  // start the fine scan from this clust.  
	  while(wfx._iclust != _cseq.clustList().end() ){
	    // First, get the starting voltage
	    wfx._vstart = sampleWaveform(StrawElectronics::thresh,wfx._iclust->time());
	    // check if this clust could cross threshold
	    if(wfx._vstart + maxLinearResponse(wfx._iclust) > threshold){
	      // check the actual response 
	      double maxtime = wfx._iclust->time() + _strawele->maxResponseTime(StrawElectronics::thresh);
	      double maxresp = sampleWaveform(StrawElectronics::thresh,maxtime);
	      if(maxresp > threshold){
		// interpolate to find the precise crossing
		fineCrossing(threshold,maxresp,wfx);

		retval = true;
		break;
	      }
	    }
	    // advance to next clust
	    ++(wfx._iclust);
	  }
	}
      }
      return retval;
    }

    void StrawWaveform::returnCrossing(double threshold, WFX& wfx) const {
      while(wfx._iclust != _cseq.clustList().end() && wfx._vstart > threshold) {
	// move forward in time at least as twice the time to the maxium for this clust
	double time = wfx._iclust->time() + 2*_strawele->maxResponseTime(StrawElectronics::thresh); 
	while(wfx._iclust != _cseq.clustList().end() && 
	    wfx._iclust->time() < time){
	  ++(wfx._iclust);
	}
	if(wfx._iclust != _cseq.clustList().end()){
	  wfx._vstart = sampleWaveform(StrawElectronics::thresh,wfx._iclust->time());
	  wfx._time =wfx._iclust->time();
	}
      }
    }

    bool StrawWaveform::roughCrossing(double threshold, WFX& wfx) const {
      // add voltage till we go over threshold by simple linear sum.  That's a minimum requirement
      // for actually crossing threshold
      double resp = wfx._vstart;
      while(wfx._iclust != _cseq.clustList().end()){
	resp += maxLinearResponse(wfx._iclust);
	if(resp > threshold)break;
	++(wfx._iclust);
      }
      // update time
      if(wfx._iclust != _cseq.clustList().end() )
	wfx._time = wfx._iclust->time();

      return wfx._iclust != _cseq.clustList().end() && resp > threshold;
    }

    bool StrawWaveform::fineCrossing(double threshold,double maxresp, WFX& wfx) const {
      static double timestep(0.020); // interpolation minimum to use linear threshold crossing calculation
      double pretime = wfx._iclust->time();
      double posttime = pretime + _strawele->maxResponseTime(StrawElectronics::thresh);
      double presample = wfx._vstart;
      double postsample = maxresp;
      static const unsigned maxstep(10); // 10 steps max
      unsigned nstep(0);
      double deltat = posttime-pretime;
      double dt = deltat;
      double slope = deltat/(postsample-presample);
      double time = pretime + slope*(threshold-presample);
      // linear interpolation
      while(fabs(dt) > timestep && nstep < maxstep) {
	double sample = sampleWaveform(StrawElectronics::thresh,time);
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
	dt = time-oldtime;
	++nstep;
      }     
      // set crossing time
      wfx._time = time;
      // record the threshold 
      wfx._vcross = threshold;
      // update the referenced clust: this can be different than the one we started with!
      while(wfx._iclust != _cseq.clustList().end() && 
	  wfx._iclust->time() < wfx._time){
	++(wfx._iclust);
      }
      // back off one
      if(wfx._iclust != _cseq.clustList().begin())--(wfx._iclust);

      // apply dispersion effects.  This changes the slope of the voltage response FIXME!!!

      // return on convergence
      return dt < timestep;
    }

    double StrawWaveform::maxLinearResponse(ClusterList::const_iterator const& iclust) const {
      // ignore saturation effects
      double linresp = _strawele->maxLinearResponse(StrawElectronics::thresh,iclust->charge());
      linresp *= (_xtalk._preamp + _xtalk._postamp);
      return linresp;
    }

    double StrawWaveform::sampleWaveform(StrawElectronics::path ipath,double time) const {
      // loop over all clusts and add their response at this time
      ClusterList const& hlist = _cseq.clustList();
      double linresp(0.0);
      auto iclust = hlist.begin();
      while(iclust != hlist.end() && iclust->time() < time){
	// compute the linear straw electronics response to this charge.  This is pre-saturation 
	linresp += _strawele->linearResponse(ipath,time-iclust->time(),iclust->charge());
	// move to next clust
	++iclust;
      }
      // apply saturation effects and x-talk
      double satresp = _strawele->saturatedResponse(linresp);
      satresp *= _xtalk._postamp;
      if(_xtalk._preamp>0.0)
	satresp += _strawele->saturatedResponse(_xtalk._preamp*linresp);
      return satresp;
    }

    void StrawWaveform::sampleWaveform(StrawElectronics::path ipath,std::vector<double> const& times,std::vector<double>& volts) const {
      volts.clear();
      volts.reserve(times.size());
      for(auto itime=times.begin();itime!=times.end();++itime){
	volts.push_back(sampleWaveform(ipath,*itime));
      }
    }

    StrawEnd StrawWaveform::strawEnd() const {
      if(!_cseq.clustList().empty())
	return _cseq.clustList().begin()->strawEnd();
      else
	return StrawEnd(StrawEnd::unknown);

    }
  }
}

