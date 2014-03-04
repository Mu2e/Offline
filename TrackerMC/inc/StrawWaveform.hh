#ifndef TrackerMC_StrawWaveform_hh
#define TrackerMC_StrawWaveform_hh
//
// StrawWaveform integrates post-amplification voltage as a function of time at one end of a
// a straw, over the time period of 1 microbunch.  It includes all physical and electronics
// effects prior to digitization.
//
// $Id: StrawWaveform.hh,v 1.11 2014/03/04 00:29:17 brownd Exp $
// $Author: brownd $
// $Date: 2014/03/04 00:29:17 $
//
// Original author David Brown, LBNL
//

// C++ includes
#include <iostream>
#include <array>
#include <vector>
#include <utility>

// Mu2e includes
#include "DataProducts/inc/StrawIndex.hh"
#include "RecoDataProducts/inc/StrawEnd.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "TrackerConditions/inc/StrawElectronics.hh"
#include "TrackerMC/inc/StrawHitletSequence.hh"

namespace mu2e {
// struct to represent cross-talk.
  struct XTalk {
    XTalk(StrawIndex const& sid) : _source(sid), _dest(sid), _preamp(0.0), _postamp(1.0)  {} // self x-talk constructor
    XTalk(StrawIndex const& source, StrawIndex const& dest,double preamp, double postamp) :
      _source(source), _dest(dest), _preamp(preamp), _postamp(postamp) {} // true x-talk constructor
    bool self() const { return _dest==_source; }
    StrawIndex _source; // index of the straw which is the source of the x-talk
    StrawIndex _dest; // index of the straw in which x-talk is observed
    double _preamp; // scaling before amplification
    double _postamp; // scaling after amplificiation
  };

  struct WFX;
  class StrawWaveform{
    public:
      // construct from a hitlet sequence and response object.  Scale affects the voltage
      StrawWaveform(StrawHitletSequence const& hseqq, ConditionsHandle<StrawElectronics> const& sresponse,XTalk const& xtalk); 
// disallow copy and assignment
      StrawWaveform() = delete; // don't allow default constructor, references can't be assigned empty
      StrawWaveform(StrawWaveform const& other);
      StrawWaveform & operator=(StrawWaveform const& other) = delete;
// find the next point the waveform crosses threhold.  Waveform crossing
// is both input (determines starting point) and output
      bool crossesThreshold(double threshold,WFX& wfx) const;
// sample the waveform at a given time.  Return value is in units of volts 
      double sampleWaveform(StrawElectronics::path ipath,double time) const;
// sample the waveform at a series of points
      void sampleWaveform(StrawElectronics::path ipath, std::vector<double> const& times,std::vector<double>& volts) const;
  //accessors
      StrawHitletSequence const& hitlets() const { return _hseq; }
      ConditionsHandle<StrawElectronics> const& strawElectronics() const { return _strawele; }
      XTalk const& xtalk() const { return _xtalk; }
      StrawEnd strawEnd() const;
    private:
// hitlet sequence used in this waveform
      StrawHitletSequence const& _hseq;
      ConditionsHandle<StrawElectronics> const& _strawele; // straw response object
      XTalk _xtalk; // X-talk applied to all voltages
// helper functions
      void returnCrossing(double threshold, WFX& wfx) const;
      bool roughCrossing(double threshold, WFX& wfx) const;
      bool fineCrossing(double threshold, double vmax, WFX& wfx) const;
      double maxLinearResponse(HitletList::const_iterator const& ihitlet) const;
};

  struct WFX { // waveform crossing
    double _time; // time waveform crossed threhold.  Note this includes noise effects
    double _vstart; // starting voltage, at time 0 of the referenced hitlet
    double _vcross; // crossing voltage
    HitletList::const_iterator _ihitlet; // iterator to hitlet associated with this crossing
    WFX() = delete; // disallow
    WFX(StrawWaveform const& wf, double time) : _time(time), _vstart(0.0), _vcross(0.0),
    _ihitlet(wf.hitlets().hitletList().begin()) {}
    // sorting function
    bool operator < (WFX const& other) { return _time < other._time; }
  };
}
#endif

