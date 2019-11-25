//Author: S Middleton
//Date: Nov 2019
//Purpose: New Digi product
#ifndef RecoDataProducts_NewCaloDigi_hh
#define RecoDataProducts_NewCaloDigi_hh

#include <vector>

namespace mu2e
{

  class NewCaloDigi {

      public:

	  NewCaloDigi(): _roId(-1), _t0(0.), _waveform(0), _peakpos(0.), _errorFlag("False"), _onspillFlag("On") {}

	  NewCaloDigi(int ROId, int t0, std::vector<int>& waveform, peakpos, errorFlag, onspillFlag ):
	    _roId(ROId),
	    _t0(t0),
	    _waveform(waveform),
	    _peakpos(peakpos),
	    _errorFlag(errorFlag),
	    _onspillFlag(onspillFlag)
	  {}

	  int                     roId()      const { return _roId;}    
	  int                     t0()        const { return _t0;}
	  const std::vector<int>& waveform()  const { return _waveform; }
	  double 		  peakpos()   const { return _peakpos;	}
	  std::string 		  errorFlag() const { return _errorFlag;			
	  std::string	          onspillFlag() const { return _onspillFlag; }
	  std::string 		  eventMode()   const return _eventMode; }
      private:

	  int               _roId;      
	  int               _t0;        //time of the first digitezd bin of the signal
	  std::vector<int>  _waveform;  //array of the samples associated with the digitezed signal
	  double	    _peakpos;	//peak position	
	  std::string	    _errorFlag; //flag for errors
	  std::string       _onspillFlag; //flag to say whether this is on or off spill data 
	  std::string       _eventMode; //gives info on event mode
  };

 
   typedef std::vector<mu2e::NewCaloDigi> NewCaloDigiCollection;

}

#endif
