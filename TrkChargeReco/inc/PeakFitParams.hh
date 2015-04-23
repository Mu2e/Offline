#ifndef PeakFitParams_hh
#define PeakFitParams_hh

// Parameter structs used by functions in FitModel.hh
#include "Rtypes.h"

namespace mu2e{

  namespace TrkChargeReco{

    struct PeakFitParams {
      enum paramIndex {earlyCharge=0,pedestal,time,scale,width,lateShift,lateScale,nParams};
      // explicit data members
      Float_t _earlyCharge; // decaying charge from earlier hit, units of pC
      Float_t _pedestal; // units of ADC counts (??)
      Float_t _time; // primary peak threshold crossing time, units of nSec WRT beam crossing (??)
      Float_t _scale;  // primary peak charge, normalized charge units (???)
      Float_t _width; // Additional width of primary peak, in units of nSec
      Float_t _lateShift; // 2nd (late) peak time WRT primary peak time, same units as _time
      Float_t _lateScale; // 2nd (late) peak scale, same units as scale
      UInt_t _free; // bitmap of free/fixed parameters

      // default constructor
      PeakFitParams(): _earlyCharge(0.0), _pedestal(0.0), _time(0.0), _scale(0.0), _width(0.0), _lateShift(0.0), _lateScale(0.0), _free(0) {}
      // explicit constructor
      PeakFitParams(Float_t e, Float_t p, Float_t t, Float_t s,Float_t w, Float_t lsh, Float_t lsc, UInt_t f) : _earlyCharge(e), _pedestal(p), _time(t),_scale(s), _width(w), _lateShift(lsh), _lateScale(lsc) , _free(f) {}
      // transform from Double_t array for root functions
      PeakFitParams(Double_t array[]) : _earlyCharge(array[earlyCharge]),
      _pedestal(array[pedestal]),
      _time(array[time]),
      _scale(array[scale]),
      _width(array[width]),
      _lateShift(array[lateShift]),
      _lateScale(array[lateScale]) {}
// convert to array for root fit
      void fillArray(Double_t array[] ) const {
	array[earlyCharge]=_earlyCharge;
	array[pedestal]=_pedestal;
	array[time]=_time;
	array[scale]=_scale;
	array[width]=_width;
	array[lateShift]=_lateShift;
	array[lateScale]=_lateScale;
      }
// access fixed/free
    bool isFree(paramIndex param) const { return (_free & (param<<1)) != 0; }
    bool isFixed(paramIndex param) const { return !isFree(param);}
    void fixParam(paramIndex param) { _free &= (~(param<<1)); } 
    void freeParam(paramIndex param) { _free |= (param<<1); } 

    };
  }  // TrkChargeReco namespace
}// mu2e namespace
#endif
