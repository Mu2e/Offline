struct resultantPeakData
{
	Double_t scalingFactor;
	Double_t peakTime;

	resultantPeakData() : scalingFactor(0.0), peakTime(0.0){};
	resultantPeakData(Double_t scalingFactor, Double_t peakTime) : scalingFactor(scalingFactor), peakTime(peakTime){};
};