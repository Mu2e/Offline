#ifndef CDFFitting_Measurement_hh
#define CDFFitting_Measurement_hh
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CDFFitting/inc/FitAction.hh"
template <class Measureable> class Fit;

//	Abstract Base Class representing the measurement of a
//	measureable object.  Subclasses must  provide:
//	residuals w.r.t. the measureable object (not necessarily
//	in Cartesian 3-space), the error matrix of the
//	measurement, and the derrivates w.r.t. the parameters of
//	the object under measurement.
//
//	In addition, the class itself can compute the chisquared
//	contribution either w.r.t. the measurable object itself,
//	or w.r.t some fit of the measureable object.  This can
//	be useful in deciding, on the basis of chi-squared or
//	confidence level, whether or not to include the
//	measurement in a fit.

template <class Measureable>
class Measurement : public FitAction<Measureable> {

public:

  // Constructor
  Measurement();

  // Destructor
  virtual ~Measurement();


  //	An abstract function that must be provided by the
  //	subclass.  The "displacement from" is the distance from
  //	the object in some vector space.  Examples:
  //
  //	A Hit is being used as a measurement on a track.  The
  //	displacement is the vector between the track and the
  //	point of closest approach of the track, in 3D cartesian
  //	space.
  //
  //	A Track is being used as a measurement of a vertex.  The
  //	displacement in this case is in the five-dimensional
  //	space of the track parameters.
  //
  //	A Track fit is being used as a measurement of a track.
  //	Again, the displacement is in the five-dimensional space
  //	of the track parameters.
  virtual HepVector getDisplacementFrom(const Measureable &measureable) const = 0;

  //	Gets the error matrix of the displacement; again, not
  //	necessarily in Cartesian 3-space (see getDisplacement).
  virtual HepSymMatrix getErrorMatrix(const Measureable &measureable) const = 0;

  //	Abstract function which gets the matrix of derivatives
  //	between the space in which the measurement takes place
  //	and the parameters space of the measureable object.
  virtual HepMatrix getDerivativeMatrix(const Measureable &measureable) const = 0;

  //	Abstract function.  Subclasses must provide thie
  //	dimensionality of the measurement.
  virtual unsigned int getDimensionality() const = 0;

  //	Calculates the Chi-squared difference between this
  //	measurement and a measurable object.
  double getDeltaChiSquared(const Measureable &measureable) const;

  //	Calculates the Chi-squared difference between this
  //	measurment and a fit to a measureable object.
  double getDeltaChiSquared(const Fit<Measureable> &fit) const;

  // Half of a double dispatch pair:
  virtual void applyYourselfTo(Fitter<Measureable> *theFitter) const;

};
#include "CDFFitting/inc/Measurement.icc"

#endif /* CDFFitting_Measurement_hh */
