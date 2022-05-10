#ifndef CaloMC_CaloWFExtractor_hh
#define CaloMC_CaloWFExtractor_hh
//
// Utility to simulate waveform hit extraction in FPGA
//
#include <vector>

namespace mu2e {

     class CaloWFExtractor
     {
           public:
               CaloWFExtractor(unsigned bufferDigi, unsigned  nBinsPeak, int minPeakADC, unsigned startOffset) :
                  bufferDigi_(bufferDigi),nBinsPeak_(nBinsPeak),minPeakADC_(minPeakADC), startOffset_(startOffset)
               {};

               void extract(const std::vector<int>& wf, std::vector<size_t>& starts, std::vector<size_t>& stops) const;

           private:
               unsigned  bufferDigi_;
               unsigned  nBinsPeak_;
               int           minPeakADC_;
               unsigned  startOffset_;
     };
}

#endif
