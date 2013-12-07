#ifndef TrackerMC_StrawWaveform_hh
#define TrackerMC_StrawWaveform_hh
//
// StrawWaveform integrates post-amplification voltage as a function of time at one end of a
// a straw, over the time period of 1 microbunch.  It includes all physical and electronics
// effects prior to digitization.
//
// $Id: StrawWaveform.hh,v 1.1 2013/12/07 19:51:42 brownd Exp $
// $Author: brownd $
// $Date: 2013/12/07 19:51:42 $
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
#include "TrackerMC/inc/StrawElectronics.hh"
#include "TrackerMC/inc/StrawHitletSequence.hh"

namespace mu2e {

  struct WFX { // waveform crossing
    double _time; // time waveform crossed threhold.  Note this includes noise effects
    HitletList::const_iterator _ihitlet; // iterator to hitlet associated with this time
    WFX(double t, HitletList::const_iterator ihitlet) : _time(t), _ihitlet(ihitlet) {}
    // sorting function
    bool operator < (WFX const& other) { return _time < other._time; }
};

  class StrawWaveform{
    public:
      enum crossdir {increasing=0,decreasing};
      StrawWaveform() = delete; // don't allow default constructor, references can't be assigned empty
      // construct from a hitlet sequence and response object
      StrawWaveform(StrawHitletSequence const& hseqq, StrawElectronics const& sresponse); 
// disallow copy and assignment, this object is too big.  In future, we could implement move
      StrawWaveform(StrawWaveform const& other) = delete;
      StrawWaveform & operator=(StrawWaveform const& other) = delete;
// find the next point the waveform crosses threhold, in either direction.  Waveform
// is both input (determines starting point) and output
      bool crossesThreshold(double threshold,crossdir cdir,WFX& wfx) const;
// sample the waveform at a given time.  Return value is in units of volts 
      double sampleWaveform(double time) const;
// sample the waveform at a series of points
      void sampleWaveform(std::vector<double> const& times,std::vector<double>& volts) const;
      StrawHitletSequence const& hitlets() const { return _hseq; }
    private:
// hitlet sequence used in this waveform
      StrawHitletSequence const& _hseq;
      StrawElectronics const& _strawele; // straw response object
  };

}
#endif

