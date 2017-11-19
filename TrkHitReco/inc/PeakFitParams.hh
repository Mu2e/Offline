#ifndef TrkHitReco_PeakFitParams_hh
#define TrkHitReco_PeakFitParams_hh

// Parameter structs used by functions in FitModel.hh
#include "Rtypes.h"
#include <string>

namespace mu2e{

  namespace TrkHitReco{

    struct PeakFitParams
    {      
      // CINT can't handle C++14 enum types FIXME!!!
      //      enum paramIndex : size_t {earlyCharge=1,pedestal,time,charge,width,lateShift,lateCharge,nParams};
      enum  paramIndex {earlyCharge=0,pedestal,time,charge,width,lateShift,lateCharge,nParams};
      
      Float_t _earlyCharge; // decaying charge from earlier hit, units of pC
      Float_t _pedestal; // units of ADC counts (??)
      Float_t _time; // primary peak threshold crossing time, units of nSec WRT beam crossing (??)
      Float_t _charge;  // primary peak charge in pC
      Float_t _width; // Additional width of primary peak, in units of nSec
      Float_t _lateShift; // 2nd (late) peak time WRT primary peak time, same units as _time
      Float_t _lateCharge; // 2nd (late) peak charge, same units as charge
      UInt_t  _free; // bitmap of free/fixed parameters
      Float_t _chi2;  // chisquared of the resultant fit
      UInt_t  _ndf; // number of degrees of freedom of fit
      Int_t   _status; // fit status

      // default constructor
      PeakFitParams(): _earlyCharge(0.0), _pedestal(0.0), _time(0.0), _charge(0.0), _width(0.0), _lateShift(0.0), _lateCharge(0.0), 
                       _free(0),_chi2(-1.0), _ndf(0), _status(-1) 
      {}

      PeakFitParams(Float_t e, Float_t p, Float_t t, Float_t s,Float_t w, Float_t lsh, Float_t lsc, UInt_t f, Float_t chi2, 
                    UInt_t ndf, Int_t status) :
	_earlyCharge(e), _pedestal(p), _time(t),_charge(s), _width(w),
	_lateShift(lsh), _lateCharge(lsc) , _free(f),
	_chi2(chi2), _ndf(ndf), _status(status) 
      {}

      // transform from Double_t array for root functions      
      PeakFitParams(Double_t array[],Double_t chi2=-1., UInt_t ndf=0, Int_t status=-1) : 
        _earlyCharge(array[earlyCharge]),
        _pedestal(array[pedestal]),
        _time(array[time]),
        _charge(array[charge]),
        _width(array[width]),
        _lateShift(array[lateShift]),
        _lateCharge(array[lateCharge]),_free(0),
        _chi2(chi2), _ndf(ndf), _status(status)
      {}
      
      // convert to array for root fit
      void fillArray(Double_t array[] ) const {
	array[earlyCharge]=_earlyCharge;
	array[pedestal]=_pedestal;
	array[time]=_time;
	array[charge]=_charge;
	array[width]=_width;
	array[lateShift]=_lateShift;
	array[lateCharge]=_lateCharge;
      }
      
      bool isFree(paramIndex param)  const { return (_free & (1<<param)) != 0; }
      bool isFixed(paramIndex param) const { return !isFree(param);}
      void fixParam(paramIndex param)      { _free &= (~(1<<param)); } 
      void freeParam(paramIndex param)     { _free |= (1<<param); } 
      static std::vector<std::string> _pnames;
      static std::string const& parameterName(paramIndex pindex) { return _pnames[pindex]; }
    };

    // add a substruct that includes parameter limit information
    struct PeakFitParamsLimits : public PeakFitParams {
      Double_t _parmin[nParams], _parmax[nParams];
    };

  } 
}
#endif
