#ifndef TrackerMC_StrawWaveform_hh
#define TrackerMC_StrawWaveform_hh
//
// StrawWaveform integrates post-amplification voltage as a function of time at one end of a
// a straw, over the time period of 1 microbunch.  It includes all physical and electronics
// effects prior to digitization.
//
// $Id: StrawWaveform.hh,v 1.5 2014/01/18 17:33:34 brownd Exp $
// $Author: brownd $
// $Date: 2014/01/18 17:33:34 $
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

  struct WFX;
  class StrawWaveform{
    public:
      // construct from a hitlet sequence and response object
      StrawWaveform(StrawHitletSequence const& hseqq, ConditionsHandle<StrawElectronics> const& sresponse); 
// disallow copy and assignment
      StrawWaveform() = delete; // don't allow default constructor, references can't be assigned empty
      StrawWaveform(StrawWaveform const& other) = delete;
      StrawWaveform & operator=(StrawWaveform const& other) = delete;
// find the next point the waveform crosses threhold.  Waveform crossing
// is both input (determines starting point) and output
      bool crossesThreshold(double threshold,WFX& wfx) const;
// sample the waveform at a given time.  Return value is in units of volts 
      double sampleWaveform(double time) const;
// sample the waveform at a series of points
      void sampleWaveform(std::vector<double> const& times,std::vector<double>& volts) const;
  //accessors
      StrawHitletSequence const& hitlets() const { return _hseq; }
      ConditionsHandle<StrawElectronics> const& strawElectronics() const { return _strawele; }
    private:
// hitlet sequence used in this waveform
      StrawHitletSequence const& _hseq;
      ConditionsHandle<StrawElectronics> const& _strawele; // straw response object
  };

  struct WFX { // waveform crossing
    double _time; // time waveform crossed threhold.  Note this includes noise effects
    HitletList::const_iterator _ihitlet; // iterator to hitlet associated with this crossing
    WFX() = delete; // disallow
    WFX(StrawWaveform const& wf,double time=0.0) : _time(time), _ihitlet(wf.hitlets().hitletList().begin())  {}
    // sorting function
    bool operator < (WFX const& other) { return _time < other._time; }
  };
}
#endif

