#ifndef RecoDataProducts_CaloDigi_hh
#define RecoDataProducts_CaloDigi_hh

#include <vector>

namespace mu2e
{
  class CaloDigi 
  {
      public:
	  CaloDigi() : 
            roId_(-1), t0_(0.), waveform_(0), peakpos_(0) 
          {}

	  CaloDigi(int ROId, int t0, const std::vector<int>& waveform, size_t peakpos):
	     roId_(ROId),t0_(t0), waveform_(waveform), peakpos_(peakpos)
	  {}

          CaloDigi(int ROId, int t0, const std::vector<int>& waveform):
	    roId_(ROId),t0_(t0),waveform_(waveform), peakpos_(0)
	  {}

	        int               roId()     const { return roId_;}    
	        int               t0()       const { return t0_;}
	        int               peakpos()  const { return peakpos_;	}
	  const std::vector<int>& waveform() const { return waveform_; }
	 
         
	private:
	  int               roId_;      
	  int               t0_;        
	  std::vector<int>  waveform_;  
	  int               peakpos_;	
  };


  typedef std::vector<mu2e::CaloDigi> CaloDigiCollection;
}

#endif
