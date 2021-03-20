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
               CaloWFExtractor(unsigned bufferDigi, int minPeakADC,unsigned  nBinsPeak) : 
                  bufferDigi_(bufferDigi),minPeakADC_(minPeakADC),nBinsPeak_(nBinsPeak)
               {};

               void extract(const std::vector<int>& wf, std::vector<unsigned>& starts, std::vector<unsigned>& stops) const;

           private:               
               unsigned  bufferDigi_;
               int  	 minPeakADC_;
               unsigned  nBinsPeak_;
     };
} 

#endif
