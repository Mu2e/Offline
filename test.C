#include <vector>
#include <iostream>
struct resultantPeakData
{
	double scalingFactor;
	double peakTime;

	resultantPeakData() : scalingFactor(0.0), peakTime(0.0){};
	resultantPeakData(double scalingFactor, double peakTime) : scalingFactor(scalingFactor), peakTime(peakTime){};
};

typedef std::vector<resultantPeakData> resultantHitData;

void testit() {

//resultantPeakData a(10.0, 20.0);
resultantPeakData b(30.0,40.0);

resultantHitData c;
c.push_back(resultantPeakData a(10.0, 20.0));
}
