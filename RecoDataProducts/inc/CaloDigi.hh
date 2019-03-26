#ifndef RecoDataProducts_CaloDigi_hh
#define RecoDataProducts_CaloDigi_hh
// Original author B. Echenard

#include <vector>

namespace mu2e
{

  class CaloDigi {

      public:

	  CaloDigi(): _roId(-1),_t0(0.),_waveform(0) {}

	  CaloDigi(int ROId, int t0, std::vector<int>& vaveform):
	    _roId(ROId),
	    _t0(t0),
	    _waveform(vaveform)
	  {}

	  int                     roId()      const { return _roId;}    
	  int                     t0()        const { return _t0;}
	  const std::vector<int>& waveform()  const { return _waveform; }



      private:

	  int               _roId;      
	  int               _t0;        //time of the first digitezd bin of the signal
	  std::vector<int>  _waveform;  //array of the samples associated with the digitezed signal
  };

 
   typedef std::vector<mu2e::CaloDigi> CaloDigiCollection;

}

#endif
