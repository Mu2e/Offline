#ifndef ConfidenceLevelComputer_h
#define ConfidenceLevelComputer_h 1
//-------------------------------------------------------
// ConfidenceLevelComputer: A class to compute confidence
// Levels.  Author:  Joe Boudreau
//-------------------------------------------------------
class ConfidenceLevelComputer  {

public:

  // Calculate a confidence level:
  static double getConfidenceLevel(double chiSquared, 
				   unsigned int degreesOfFreedom);


private:
  
  static inline float gammp(float,float);
  static inline float gammq(float,float);
  static inline float gamser(float, float, float);
  static inline float gammcf(float, float, float);
  static const unsigned int ITMAX;
  static const double       EPS;
  static const double       FPMIN;

};

#include "CDFFitting/inc/ConfidenceLevelComputer.icc"

#endif


