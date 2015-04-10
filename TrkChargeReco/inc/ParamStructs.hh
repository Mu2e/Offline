#ifndef ParamStructs_hh
#define ParamStructs_hh

// Parameter structs used by functions in FitModel.hh

namespace mu2e{

namespace TrkChargeReco{

struct EarlyPeakParamStruct{
   Double_t _Q;

   EarlyPeakParamStruct(){};

   EarlyPeakParamStruct(Double_t Q) : _Q(Q){};
};

struct SinglePeakParamStruct
{
  Double_t _shiftedTime;
  Double_t _scalingFactor;
  Double_t _sigma;

  SinglePeakParamStruct(){};

  SinglePeakParamStruct(Double_t shiftedTime, Double_t scalingFactor, Double_t sigma = 0.0) : _shiftedTime(shiftedTime), _scalingFactor(scalingFactor), _sigma(sigma){};
};

struct SinglePeakFloatingPedestalParamStruct : SinglePeakParamStruct
{
  Double_t _verticalShift;

  SinglePeakFloatingPedestalParamStruct(){};

  SinglePeakFloatingPedestalParamStruct(Double_t shiftedTime, Double_t scalingFactor, Double_t verticalShift, Double_t sigma) :
  			SinglePeakParamStruct(shiftedTime, scalingFactor, sigma), _verticalShift(verticalShift){};
};

struct EXPeakParamStruct : SinglePeakParamStruct, EarlyPeakParamStruct
{
  EXPeakParamStruct(){};

	EXPeakParamStruct(Double_t shiftedTime, Double_t scalingFactor, Double_t Q, Double_t sigma) : 
			SinglePeakParamStruct(shiftedTime, scalingFactor, sigma), EarlyPeakParamStruct(Q){};
};

struct LXPeakParamStruct
{
  Double_t _shiftedTimeFirstPeak;
  Double_t _scalingFactorFirstPeak;
  Double_t _shiftedTimeSecondPeak;
  Double_t _scalingFactorSecondPeak;

  LXPeakParamStruct(){};

  LXPeakParamStruct(Double_t shiftedTimeFirstPeak, Double_t scalingFactorFirstPeak, 
  						Double_t shiftedTimeSecondPeak, Double_t scalingFactorSecondPeak) : 
  	_shiftedTimeFirstPeak(shiftedTimeFirstPeak), _scalingFactorFirstPeak(scalingFactorFirstPeak), 
  	_shiftedTimeSecondPeak(shiftedTimeSecondPeak), _scalingFactorSecondPeak(scalingFactorSecondPeak){};
};


struct LXPeakFloatingPedestalParamStruct : LXPeakParamStruct
{
  Double_t _verticalShift;

  LXPeakFloatingPedestalParamStruct(){};
  LXPeakFloatingPedestalParamStruct(Double_t shiftedTimeFirstPeak, Double_t scalingFactorFirstPeak, Double_t verticalShift,
  						Double_t shiftedTimeSecondPeak, Double_t scalingFactorSecondPeak) :
  	LXPeakParamStruct(shiftedTimeFirstPeak, scalingFactorFirstPeak, shiftedTimeSecondPeak, scalingFactorSecondPeak), _verticalShift(verticalShift){};
};

struct DoublePeakWithEarlyPeakParamStruct : LXPeakParamStruct, EarlyPeakParamStruct
{
  DoublePeakWithEarlyPeakParamStruct(){};

	DoublePeakWithEarlyPeakParamStruct(Double_t shiftedTimeFirstPeak, Double_t scalingFactorFirstPeak, Double_t Q,
  						Double_t shiftedTimeSecondPeak, Double_t scalingFactorSecondPeak ) :
	LXPeakParamStruct(shiftedTimeFirstPeak, scalingFactorFirstPeak, shiftedTimeSecondPeak, scalingFactorSecondPeak), 
	EarlyPeakParamStruct(Q){}
};
}
}
#endif
