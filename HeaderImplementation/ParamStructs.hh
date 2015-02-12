#ifndef ParamStructs_hh
#define ParamStructs_hh

struct DynamicPedestalParamStruct{
   Double_t _Q;

   DynamicPedestalParamStruct(){};

   DynamicPedestalParamStruct(Double_t Q) : _Q(Q){};
};

struct SinglePeakParamStruct
{
  Double_t _shiftedTime;
  Double_t _scalingFactor;
  Double_t _sigma;

  SinglePeakParamStruct(){};

  SinglePeakParamStruct(Double_t shiftedTime, Double_t scalingFactor, Double_t sigma = 0.0) : _shiftedTime(shiftedTime), _scalingFactor(scalingFactor), _sigma(sigma){};
};

struct SinglePeakWithConstantPedestalParamStruct : SinglePeakParamStruct
{
  Double_t _verticalShift;

  SinglePeakWithConstantPedestalParamStruct(){};

  SinglePeakWithConstantPedestalParamStruct(Double_t shiftedTime, Double_t scalingFactor, Double_t verticalShift, Double_t sigma) :
  			SinglePeakParamStruct(shiftedTime, scalingFactor, sigma), _verticalShift(verticalShift){};
};

struct SinglePeakWithDynamicPedestalParamStruct : SinglePeakParamStruct, DynamicPedestalParamStruct
{
  SinglePeakWithDynamicPedestalParamStruct(){};

	SinglePeakWithDynamicPedestalParamStruct(Double_t shiftedTime, Double_t scalingFactor, Double_t Q, Double_t sigma) : 
			SinglePeakParamStruct(shiftedTime, scalingFactor, sigma), DynamicPedestalParamStruct(Q){};
};

struct DoublePeakParamStruct
{
  Double_t _shiftedTimeFirstPeak;
  Double_t _scalingFactorFirstPeak;
  Double_t _shiftedTimeSecondPeak;
  Double_t _scalingFactorSecondPeak;

  DoublePeakParamStruct(){};

  DoublePeakParamStruct(Double_t shiftedTimeFirstPeak, Double_t scalingFactorFirstPeak, 
  						Double_t shiftedTimeSecondPeak, Double_t scalingFactorSecondPeak) : 
  	_shiftedTimeFirstPeak(shiftedTimeFirstPeak), _scalingFactorFirstPeak(scalingFactorFirstPeak), 
  	_shiftedTimeSecondPeak(shiftedTimeSecondPeak), _scalingFactorSecondPeak(scalingFactorSecondPeak){};
};


struct DoublePeakWithConstantPedestalParamStruct : DoublePeakParamStruct
{
  Double_t _verticalShift;

  DoublePeakWithConstantPedestalParamStruct(){};
  DoublePeakWithConstantPedestalParamStruct(Double_t shiftedTimeFirstPeak, Double_t scalingFactorFirstPeak, Double_t verticalShift,
  						Double_t shiftedTimeSecondPeak, Double_t scalingFactorSecondPeak) :
  	DoublePeakParamStruct(shiftedTimeFirstPeak, scalingFactorFirstPeak, shiftedTimeSecondPeak, scalingFactorSecondPeak), _verticalShift(verticalShift){};
};

struct DoublePeakWithDynamicPedestalParamStruct : DoublePeakParamStruct, DynamicPedestalParamStruct
{
  DoublePeakWithDynamicPedestalParamStruct(){};

	DoublePeakWithDynamicPedestalParamStruct(Double_t shiftedTimeFirstPeak, Double_t scalingFactorFirstPeak, Double_t Q,
  						Double_t shiftedTimeSecondPeak, Double_t scalingFactorSecondPeak ) :
	DoublePeakParamStruct(shiftedTimeFirstPeak, scalingFactorFirstPeak, shiftedTimeSecondPeak, scalingFactorSecondPeak), 
	DynamicPedestalParamStruct(Q){}
};

#endif