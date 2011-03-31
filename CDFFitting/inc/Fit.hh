#ifndef FIT_HH_
#define FIT_HH_
#include "CDFFitting/inc/Measurement.hh"
#include "CDFFitting/inc/ConfidenceLevelComputer.hh"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Matrix.h"

//	The Fit is the output of a fitter.  It derrives from
//	measurement, and so is also usable in a subsequent fit.
//
//	The methods getDerivativeMatrix, get ErrorMatrix,
//	getDimensionality and getDisplacement from are
//	calculated automatically and are necessary for this
//	class to be incorporated in a subsequent fit step.
//
//	The other methods are self-explanatory.

template <class Measureable>
class Fit : public Measurement<Measureable> {

public:
  // Constructor
  Fit();
  
  // Copy constructor
  Fit(const Fit<Measureable> &right);

  // Constructor
      Fit(double chiSquared,
	  unsigned  int degreesOfFreedom,
	  const HepVector     &estimate,
	  const HepSymMatrix &covariance);

  // Destructor
  ~Fit();
  
  // Assignment
  const Fit<Measureable> & operator=(const Fit<Measureable> &right);
  
  //	Retrieves the parameters of the fit, i.e., the estimate.
  const HepVector & getParameters() const;
  
  //	Return a non-const reference to the fit parameters.
  //    WARNING:  modifying parameters can lead to an inconsistent fit
  HepVector & parameters();
  
  //	Returns the number of degrees of freedom
  unsigned int getNumDegreesOfFreedom() const;
  
  //    Return a non-const reference to the number of degrees of freedom
  //    WARNING:  modifying degrees of freedom can lead to an inconsistent fit
  unsigned int & numDegreesOfFreedom();

  //	Returns the Chi-Squared of  the fit.
  double getChiSquared() const;
  
  //    Return a non-const reference to the chisquared
  //    WARNING:  modifying chisquared can lead to an inconsistent fit
  double & chiSquared();

  //	Returns the confidence level of the fit.
  double getProbability() const;
  
  //	For use in subsequent fitting.
  HepVector getDisplacementFrom(const Measureable &measureable) const;
  
  //	Gets the error matrix of the fit.
  HepSymMatrix getErrorMatrix(const Measureable &measureable) const;
  
  //    Return the error matrix 
  const HepSymMatrix & getErrorMatrix() const;

  //    Return a non-const reference to the error matrix.
  //    WARNING:  modifying error matrix can lead to an inconsistent fit.
  HepSymMatrix & errorMatrix();

  //	For use in subsequent fitting.  Returns the identity
  //	matrix.
  HepMatrix getDerivativeMatrix(const Measureable &measureable) const;
  
  //	Returns the dimensionality of the parameter space.
  //	Needed for subsequent fitting.  May also be generally
  //	useful.
  unsigned int getDimensionality() const;
  
  Measureable getMeasureable() const;

private:

  unsigned  int                        _degreesOfFreedom;
  double                               _chiSquared;
  HepVector                            _estimate;
  HepSymMatrix                        _covariance;
  HepSymMatrix                        _derivativeMatrix;

};

#include "CDFFitting/inc/Fit.icc"


#endif //FIT_HH_
