#ifndef PeakFitParams_hh
#define PeakFitParams_hh

// Parameter structs used by functions in FitModel.hh
#include "Rtypes.h"

namespace mu2e{

  namespace TrkChargeReco{

    struct PeakFitParams {
      enum paramIndex {earlyCharge=0,pedestal,time,scale,width,lateShift,lateScale,nParams};
      // explicit data members
      Double_t _earlyCharge; // decaying charge from earlier hit, units of pC
      Double_t _pedestal; // units of ADC counts (??)
      Double_t _time; // primary peak threshold crossing time, units of nSec WRT beam crossing (??)
      Double_t _scale;  // primary peak charge, normalized charge units (???)
      Double_t _width; // Additional width of primary peak, in units of nSec
      Double_t _lateShift; // 2nd (late) peak time WRT primary peak time, same units as _time
      Double_t _lateScale; // 2nd (late) peak scale, same units as scale
      // transform from/to Double_t array for root functions
      PeakFitParams(Double_t array[]) : _earlyCharge(array[earlyCharge]),
      _pedestal(array[pedestal]),
      _time(array[time]),
      _scale(array[scale]),
      _width(array[width]),
      _lateShift(array[lateShift]),
      _lateScale(array[lateScale]) {}

      void fillArray(Double_t array[] ) const {
	array[earlyCharge]=_earlyCharge;
	array[pedestal]=_pedestal;
	array[time]=_time;
	array[scale]=_scale;
	array[width]=_width;
	array[lateShift]=_lateShift;
	array[lateScale]=_lateScale;
      }


    };
  }  // TrkChargeReco namespace
}// mu2e namespace
#endif
