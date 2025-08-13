#ifndef RecoDataProducts_CaloDigi_hh
#define RecoDataProducts_CaloDigi_hh

#include <vector>
#include <cstddef>

namespace mu2e
{
  class CaloDigi
  {
      public:
          CaloDigi() :
            SiPMID_(-1), t0_(0.), waveform_(0), peakpos_(0)
          {}

          CaloDigi(int SiPMID, int t0, const std::vector<int>& waveform, size_t peakpos):
             SiPMID_(SiPMID),t0_(t0), waveform_(waveform), peakpos_(peakpos)
          {}

          CaloDigi(int SiPMID, int t0, const std::vector<int>& waveform):
            SiPMID_(SiPMID),t0_(t0),waveform_(waveform), peakpos_(0)
          {}

          int                     SiPMID()   const {return SiPMID_;}
          int                     t0()       const {return t0_;}
          int                     peakpos()  const {return peakpos_;}
          const std::vector<int>& waveform() const {return waveform_;}


        private:
          int               SiPMID_;
          int               t0_;
          std::vector<int>  waveform_;
          int               peakpos_;
  };


  using CaloDigiCollection = std::vector<mu2e::CaloDigi>;
}

#endif
