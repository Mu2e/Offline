// Original author B. Echenard, editted by S Middleton
//Date: Nov 2019
//Purpose: New Digi product
#ifndef RecoDataProducts_CaloDigi_hh
#define RecoDataProducts_CaloDigi_hh

#include <vector>

namespace mu2e
{

  class CaloDigi {

      public:
  
	  CaloDigi(): _roId(-1), _t0(0.), _waveform(0), _peakpos(0.){}
      
	 CaloDigi(int ROId, int t0, std::vector<int>& waveform, size_t peakpos):
	    _roId(ROId),
	    _t0(t0),
	    _waveform(waveform),
	    _peakpos(peakpos){}

          //For schema evolution:
          CaloDigi(int ROId, int t0, std::vector<int>& waveform):
	    _roId(ROId),
	    _t0(t0),
	    _waveform(waveform)
	  {}

	  int                     roId()      const { return _roId;}    
	  int                     t0()        const { return _t0;}
	  const std::vector<int>& waveform()  const { return _waveform; }
	  size_t 	  peakpos()   const { return _peakpos;	}
	 

	private:

	  int               _roId;      
	  int               _t0;        //time of the first digitezd bin of the signal
	  std::vector<int>  _waveform;  //array of the samples associated with the digitezed signal
	  size_t     _peakpos;	//peak position	for fast estimate of total charge and hit time
  };

   typedef std::vector<mu2e::CaloDigi> CaloDigiCollection;

}

#endif
