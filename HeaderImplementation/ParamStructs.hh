#ifndef ParamStructs_hh
#define ParamStructs_hh

struct dynamicPedestalParamStruct{
   Double_t _Q;

   dynamicPedestalParamStruct(){};

   dynamicPedestalParamStruct(Double_t Q) : _Q(Q){};
};

struct singlePeakParamStruct
{
  Double_t _shiftedTime;
  Double_t _scalingFactor;
  Double_t _sigma;

  singlePeakParamStruct(){};

  singlePeakParamStruct(Double_t shiftedTime, Double_t scalingFactor, Double_t sigma = 0.0) : _shiftedTime(shiftedTime), _scalingFactor(scalingFactor), _sigma(sigma){};
};

struct singlePeakWithConstantPedestalParamStruct : singlePeakParamStruct
{
  Double_t _verticalShift;

  singlePeakWithConstantPedestalParamStruct(){};

  singlePeakWithConstantPedestalParamStruct(Double_t shiftedTime, Double_t scalingFactor, Double_t verticalShift, Double_t sigma) :
  			singlePeakParamStruct(shiftedTime, scalingFactor, sigma), _verticalShift(verticalShift){};
};

struct singlePeakWithDynamicPedestalParamStruct : singlePeakParamStruct, dynamicPedestalParamStruct
{
  singlePeakWithDynamicPedestalParamStruct(){};

	singlePeakWithDynamicPedestalParamStruct(Double_t shiftedTime, Double_t scalingFactor, Double_t Q, Double_t sigma) : 
			singlePeakParamStruct(shiftedTime, scalingFactor, sigma), dynamicPedestalParamStruct(Q){};
};

struct doublePeakParamStruct
{
  Double_t _shiftedTimeFirstPeak;
  Double_t _scalingFactorFirstPeak;
  Double_t _shiftedTimeSecondPeak;
  Double_t _scalingFactorSecondPeak;

  doublePeakParamStruct(){};

  doublePeakParamStruct(Double_t shiftedTimeFirstPeak, Double_t scalingFactorFirstPeak, 
  						Double_t shiftedTimeSecondPeak, Double_t scalingFactorSecondPeak) : 
  	_shiftedTimeFirstPeak(shiftedTimeFirstPeak), _scalingFactorFirstPeak(scalingFactorFirstPeak), 
  	_shiftedTimeSecondPeak(shiftedTimeSecondPeak), _scalingFactorSecondPeak(scalingFactorSecondPeak){};
};


struct doublePeakWithConstantPedestalParamStruct : doublePeakParamStruct
{
  Double_t _verticalShift;

  doublePeakWithConstantPedestalParamStruct(){};
  doublePeakWithConstantPedestalParamStruct(Double_t shiftedTimeFirstPeak, Double_t scalingFactorFirstPeak, Double_t verticalShift,
  						Double_t shiftedTimeSecondPeak, Double_t scalingFactorSecondPeak) :
  	doublePeakParamStruct(shiftedTimeFirstPeak, scalingFactorFirstPeak, shiftedTimeSecondPeak, scalingFactorSecondPeak), _verticalShift(verticalShift){};
};

struct doublePeakWithDynamicPedestalParamStruct : doublePeakParamStruct, dynamicPedestalParamStruct
{
  doublePeakWithDynamicPedestalParamStruct(){};

	doublePeakWithDynamicPedestalParamStruct(Double_t shiftedTimeFirstPeak, Double_t scalingFactorFirstPeak, Double_t Q,
  						Double_t shiftedTimeSecondPeak, Double_t scalingFactorSecondPeak ) :
	doublePeakParamStruct(shiftedTimeFirstPeak, scalingFactorFirstPeak, shiftedTimeSecondPeak, scalingFactorSecondPeak), 
	dynamicPedestalParamStruct(Q){}
};

#endif