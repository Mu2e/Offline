//
// StrawWaveform represents post-amplification voltage as a function of time at one end of a
// a straw, over the time period of 1 microbunch.  It includes all physical and electronics
// effects prior to digitization.
//
// Original author David Brown, LBNL
//
#include "TrackerMC/inc/StrawWaveform.hh"
#include <math.h>
#include <boost/math/special_functions/binomial.hpp>

using namespace std;
namespace mu2e {
  using namespace TrkTypes;
  namespace TrackerMC {
    StrawWaveform::StrawWaveform(StrawClusterSequence const& hseq, XTalk const& xtalk) :
      _cseq(hseq), _xtalk(xtalk), _sid(hseq.strawId())
    {}

    StrawWaveform::StrawWaveform(StrawWaveform const& other) : _cseq(other._cseq),
    _xtalk(other._xtalk), _sid(other._sid)
    {}

    bool StrawWaveform::crossesThreshold(StrawElectronics const& strawele,double threshold,WFX& wfx) const {
      bool retval(false);
      // make sure we start past the input time
      while(wfx._iclust != _cseq.clustList().end() && wfx._iclust->time()-strawele.clusterLookbackTime()< wfx._time ){
	++(wfx._iclust);
      }
      // loop till we're at the end or we go over threshold
      if(wfx._iclust != _cseq.clustList().end()){
	// sample initial voltage for this clust
	wfx._vstart = sampleWaveform(strawele,StrawElectronics::thresh,wfx._iclust->time()-strawele.clusterLookbackTime());
	// if we start above threhold, scan forward till we're below
	if(wfx._vstart > threshold)
	  returnCrossing(strawele,threshold, wfx);
	// scan ahead quickly from there to where there's enough integral charge to possibly cross threshold.
	if(roughCrossing(strawele,threshold,wfx)) {
	  // start the fine scan from this clust.
	  while(wfx._iclust != _cseq.clustList().end() ){
	    // First, get the starting voltage
	    wfx._vstart = sampleWaveform(strawele,StrawElectronics::thresh,wfx._iclust->time()-strawele.clusterLookbackTime());
	    //// check if this clust could cross threshold
	    //if(wfx._vstart + maxLinearResponse(wfx._iclust) > threshold){
	      // check the actual response
	      double maxtime = wfx._iclust->time()+strawele.maxResponseTime(StrawElectronics::thresh,wfx._iclust->wireDistance());
	      double maxresp = sampleWaveform(strawele,StrawElectronics::thresh,maxtime);
	      if(maxresp > threshold){
		// interpolate to find the precise crossing
		fineCrossing(strawele,threshold,maxresp,wfx);

		retval = true;
		break;
	      }
            //}
	    // advance to next clust
	    ++(wfx._iclust);
	  }
	}
      }
      return retval;
    }

    void StrawWaveform::returnCrossing(StrawElectronics const& strawele, double threshold, WFX& wfx) const {
      while(wfx._iclust != _cseq.clustList().end() && wfx._vstart > threshold) {
	// move forward in time at least as twice the time to the maxium for this clust
	double time = wfx._iclust->time()+strawele.clusterLookbackTime() + 2*strawele.maxResponseTime(StrawElectronics::thresh,wfx._iclust->wireDistance());
	while(wfx._iclust != _cseq.clustList().end() &&
	    wfx._iclust->time()-strawele.clusterLookbackTime() < time){
	  ++(wfx._iclust);
	}
	if(wfx._iclust != _cseq.clustList().end()){
	  wfx._vstart = sampleWaveform(strawele,StrawElectronics::thresh,wfx._iclust->time()-strawele.clusterLookbackTime());
	  wfx._time =wfx._iclust->time()-strawele.clusterLookbackTime();
	}
      }
    }

    bool StrawWaveform::roughCrossing(StrawElectronics const& strawele, double threshold, WFX& wfx) const {
      // add voltage till we go over threshold by simple linear sum.  That's a minimum requirement
      // for actually crossing threshold
      double resp = wfx._vstart;
      while(wfx._iclust != _cseq.clustList().end()){
	resp += maxLinearResponse(strawele,wfx._iclust);
	if(resp > threshold)break;
	++(wfx._iclust);
      }
      // update time
      if(wfx._iclust != _cseq.clustList().end() )
	wfx._time = wfx._iclust->time()-strawele.clusterLookbackTime();

      return wfx._iclust != _cseq.clustList().end() && resp > threshold;
    }

    bool StrawWaveform::fineCrossing(StrawElectronics const& strawele, double threshold,double maxresp, WFX& wfx) const {
      static double timestep(0.020); // interpolation minimum to use linear threshold crossing calculation
      double pretime = wfx._iclust->time()-strawele.clusterLookbackTime();
      double posttime = pretime + strawele.clusterLookbackTime() + strawele.maxResponseTime(StrawElectronics::thresh,wfx._iclust->wireDistance());
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
	double sample = sampleWaveform(strawele,StrawElectronics::thresh,time);
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
	  wfx._iclust->time()-strawele.clusterLookbackTime() < wfx._time){
	++(wfx._iclust);
      }
      // back off one
      if(wfx._iclust != _cseq.clustList().begin())--(wfx._iclust);

      // apply dispersion effects.  This changes the slope of the voltage response FIXME!!!

      // return on convergence
      return dt < timestep;
    }

    double StrawWaveform::maxLinearResponse(StrawElectronics const& strawele,ClusterList::const_iterator const& iclust) const {
      // ignore saturation effects
      double linresp = strawele.maxLinearResponse(_sid,StrawElectronics::thresh,iclust->wireDistance(),iclust->charge());
      linresp *= (_xtalk._preamp + _xtalk._postamp);
      return linresp;
    }

    double StrawWaveform::sampleWaveform(StrawElectronics const& strawele,StrawElectronics::Path ipath,double time) const {
      // loop over all clusts and add their response at this time
      ClusterList const& hlist = _cseq.clustList();
      double linresp(0.0);
      auto iclust = hlist.begin();
      while(iclust != hlist.end() && iclust->time()-strawele.clusterLookbackTime() < time){
	// compute the linear straw electronics response to this charge.  This is pre-saturation
	linresp += strawele.linearResponse(_sid,ipath,time-iclust->time(),iclust->charge(),iclust->wireDistance());
	// move to next clust
	++iclust;
      }
      double totresp = linresp * _xtalk._postamp;
      if(_xtalk._preamp>0.0)
	totresp += _xtalk._preamp*linresp;
      return totresp;
    }

    void StrawWaveform::sampleADCWaveform(StrawElectronics const& strawele,ADCTimes const& times,ADCVoltages& volts) const {
      volts.clear();
      volts.reserve(times.size());
      if (_xtalk._dest != _xtalk._source){
        //FIXME doesn't deal with cross talk hits yet
        for (size_t j=0;j<times.size();j++){
          volts.push_back(0);
        }
        return;
      }

      // check if going to be saturated
      double max_possible_voltage = 0;
      for (auto iclust = _cseq.clustList().begin();iclust != _cseq.clustList().end();++iclust){
        max_possible_voltage += maxLinearResponse(strawele,iclust);
      }
      if (max_possible_voltage > strawele.saturationVoltage()){
        // create waveform of threshold circuit output
        // step along waveform and apply saturation
        // for each time, get contribution from each step in waveform using impulse response

        // skip to the first cluster that matters for the first adc time
        auto iclust = _cseq.clustList().begin();
        while (iclust != _cseq.clustList().end()){
          double time = iclust->time()-strawele.clusterLookbackTime();
          if (time + strawele.truncationTime(StrawElectronics::thresh) > times[0])
            break;
          else
            ++iclust;
        }

        // step through time
        for (size_t j=0;j<times.size();j++){
          volts.push_back(0);
        }

        int num_steps = (int)ceil((times[times.size()-1]-iclust->time()-strawele.clusterLookbackTime())/strawele.saturationTimeStep());

        for (int i=0;i<num_steps;i++){
          double time = iclust->time()-strawele.clusterLookbackTime() + i*strawele.saturationTimeStep();
          // sum up the preamp response at this step
          double response = 0;
          auto jclust = iclust;
          while(jclust != _cseq.clustList().end() && jclust->time()-strawele.clusterLookbackTime() < time){
            response += strawele.linearResponse(_sid,StrawElectronics::thresh,time-jclust->time(),jclust->charge(),jclust->wireDistance(),true);
            ++jclust;
          }
          // now saturate it
          double sat_response = strawele.saturatedResponse(response);
          // then calculate the impulse response at each of the adctimes and add it to that
          for (size_t j=0;j<times.size();j++){
            // this function includes multiplication by number of steps in saturationTimeStep
            volts[j] += strawele.adcImpulseResponse(_sid,times[j]-time,sat_response);
          }
        }
      }else{
        for(auto itime=times.begin();itime!=times.end();++itime){
          volts.push_back(sampleWaveform(strawele,StrawElectronics::adc,*itime));
        }
      }
    }

  unsigned short StrawWaveform::digitizeTOT(StrawElectronics const& strawele, double threshold, double time) const {
      for (size_t i=1;i<strawele.maxTOT();i++){
        if (sampleWaveform(strawele,StrawElectronics::thresh,time + i*strawele.totLSB()) < threshold)
          return static_cast<unsigned short>(i);
      }
      return static_cast<unsigned short>(strawele.maxTOT());
    }
  }


}
