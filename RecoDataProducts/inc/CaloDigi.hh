#ifndef RecoDataProducts_CaloDigi_hh
#define RecoDataProducts_CaloDigi_hh

#include <vector>
#include <cstddef>

namespace mu2e
{
  class CaloDigi
  {
      public:
          // Format version of the digi content, stored per digi. Bump it whenever the meaning of a member
          // changes, and record the change here; consumers branch on format() or reject formats they cannot read.
          //   0  written before this member existed (ROOT schema evolution leaves the default):
          //      t0 units are not recorded (ns before Offline#2022, ticks after)
          //   1  t0 in digitizer clock ticks (Offline#2022)
          static constexpr int currentFormat = 1;

          CaloDigi() :
            SiPMID_(-1), t0_(0.), waveform_(0), peakpos_(0)
          {}

          CaloDigi(int SiPMID, int t0, const std::vector<int>& waveform, size_t peakpos):
             SiPMID_(SiPMID),t0_(t0), waveform_(waveform), peakpos_(peakpos), format_(currentFormat)
          {}

          // for code that rebuilds a digi from existing ones: keep the format of the input, not currentFormat
          CaloDigi(int SiPMID, int t0, const std::vector<int>& waveform, size_t peakpos, int format):
             SiPMID_(SiPMID),t0_(t0), waveform_(waveform), peakpos_(peakpos), format_(format)
          {}

          CaloDigi(int SiPMID, int t0, const std::vector<int>& waveform):
            SiPMID_(SiPMID),t0_(t0),waveform_(waveform), peakpos_(0), format_(currentFormat)
          {}

          int                     SiPMID()   const {return SiPMID_;}
          // start of the waveform in digitizer clock ticks (CaloConst::_digitizationPeriod, 5 ns), in the
          // digitizer (DR marker) frame; the raw hit-packet Time. CaloRecoDigiMaker converts to ns.
          int                     t0()       const {return t0_;}
          int                     peakpos()  const {return peakpos_;}
          const std::vector<int>& waveform() const {return waveform_;}
          int                     format()   const {return format_;}


        private:
          int               SiPMID_;
          int               t0_;
          std::vector<int>  waveform_;
          int               peakpos_;
          int               format_ = 0;
  };


  using CaloDigiCollection = std::vector<mu2e::CaloDigi>;
}

#endif
