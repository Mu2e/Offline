#ifndef TrackerMC_StrawWaveform_hh
#define TrackerMC_StrawWaveform_hh
//
// StrawWaveform integrates post-amplification voltage as a function of time at one end of a
// a straw, over the time period of 1 microbunch.  It includes all physical and electronics
// effects prior to digitization.
//
// Original author David Brown, LBNL
//

// C++ includes
#include <iostream>
#include <array>
#include <vector>
#include <utility>

// Mu2e includes
#include "Offline/DataProducts/inc/StrawEnd.hh"
#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/TrackerConditions/inc/StrawElectronics.hh"
#include "Offline/DataProducts/inc/TrkTypes.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"
#include "Offline/TrackerMC/inc/StrawClusterSequence.hh"

namespace mu2e {
  namespace TrackerMC {
    // struct to represent cross-talk.
    struct XTalk {
      XTalk(StrawId const& sid) : _source(sid), _dest(sid), _preamp(0.0), _postamp(1.0)  {} // self x-talk constructor
      XTalk(StrawId const& source, StrawId const& dest,double preamp, double postamp) :
        _source(source), _dest(dest), _preamp(preamp), _postamp(postamp) {} // true x-talk constructor
      bool self() const { return _dest==_source; }
      StrawId _source; // index of the straw which is the source of the x-talk
      StrawId _dest; // index of the straw in which x-talk is observed
      double _preamp; // scaling before amplification
      double _postamp; // scaling after amplificiation
    };

    struct WFX;
    class StrawWaveform{
      public:
        // construct from a clust sequence and response object.  Scale affects the voltage
        StrawWaveform(Straw const& straw, StrawClusterSequence const& hseqq, XTalk const& xtalk);
        // disallow copy and assignment
        StrawWaveform() = delete; // don't allow default constructor, references can't be assigned empty
        StrawWaveform(StrawWaveform const& other);
        StrawWaveform & operator=(StrawWaveform const& other) = delete;
        // find the next point the waveform crosses threhold.  Waveform crossing
        // is both input (determines starting point) and output
        bool crossesThreshold(StrawElectronics const& strawele, double threshold,WFX& wfx) const;
        // sample the waveform at a given time, no saturation included.  Return value is in units of volts
        double sampleWaveform(StrawElectronics const& strawele,StrawElectronics::Path ipath,double time) const;
        // sample the waveform at a series of points allowing saturation to occur after preamp stage
        // FIXME no cross talk yet
        void sampleADCWaveform(StrawElectronics const& strawele,TrkTypes::ADCTimes const& times,TrkTypes::ADCVoltages& volts) const;
        unsigned short digitizeTOT(StrawElectronics const& strawele, double threshold, double time) const;
        //accessors
        StrawClusterSequence const& clusts() const { return _cseq; }
        XTalk const& xtalk() const { return _xtalk; }
        StrawEnd const& strawEnd() const { return _cseq.strawEnd(); }
        Straw const& straw() const { return _straw;}
      private:
        // clust sequence used in this waveform
        StrawClusterSequence const& _cseq;
        XTalk _xtalk; // X-talk applied to all voltages
        Straw const& _straw;
        // helper functions
        void returnCrossing(StrawElectronics const& strawele, double threshold, WFX& wfx) const;
        bool roughCrossing(StrawElectronics const& strawele, double threshold, WFX& wfx) const;
        bool fineCrossing(StrawElectronics const& strawele, double threshold, double vmax, WFX& wfx) const;
        double maxLinearResponse(StrawElectronics const& strawele,StrawClusterList::const_iterator const& iclust) const;
    };

    struct WFX { // waveform crossing
      double _time; // time waveform crossed threhold.  Note this includes noise effects
      double _vstart; // starting voltage, at time 0 of the referenced clust
      double _vcross; // crossing voltage
      StrawClusterList::const_iterator _iclust; // iterator to clust associated with this crossing
      WFX() = delete; // disallow
      WFX(StrawWaveform const& wf, double time) : _time(time), _vstart(0.0), _vcross(0.0),
      _iclust(wf.clusts().clustList().begin()) {}
      // sorting function
      bool operator < (WFX const& other) { return _time < other._time; }
    };
  }
}
#endif
