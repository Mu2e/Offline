#ifndef CDFFitting_LinearFitter_hh
#define CDFFitting_LinearFitter_hh
#include <vector>
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CDFFitting/inc/Fitter.hh"
#include "CDFFitting/inc/Fit.hh"
#include "CDFFitting/inc/Residual.hh"
#include "CDFFitting/inc/Constraint.hh"

template <class Measureable>
class LinearFitter : public Fitter<Measureable> {

public:

  // Constructor
  LinearFitter();

  // Copy Constructor
  LinearFitter(const LinearFitter<Measureable> &right);

  // Destructor
  virtual ~LinearFitter();

  // Assignment operator
  const LinearFitter<Measureable> & operator=(const LinearFitter<Measureable> &right);

  // Access number of applied fit actions
  virtual unsigned int getNumAppliedFitActions() const;

  // Fit
  virtual void fit();

  int status() const {return _status;}

  // Return the fit result
  virtual const Fit<Measureable> * getFit() const;

  // Return a residual.. Client is responsible for deleting it
  virtual const Residual * newResidual(const Measurement<Measureable> *measurement) ;

  // Reset the fit.
  virtual void resetFit();

private:

  int                  _degreesOfFreedom;
  unsigned int         _numAppliedFitActions;
  HepMatrix            _derivatives;
  HepSymMatrix        _weight;
  HepVector            _residuals;
  HepMatrix            _constraintMatrices;
  HepVector            _constraintVectors;
  Fit<Measureable> *   _fit;
  vector<unsigned int> _measurementPtr;
  int _status;

protected:
  // Apply the various actions; only Measurement & constraint actually do anything
  // for a linear fitter.
  virtual void apply(const Measurement<Measureable> *theMeasurement);
  virtual void apply(const Constraint<Measureable> *theConstraint);
  virtual void apply(const Scatter<Measureable> *theScatter);
  virtual void apply(const Transport<Measureable> *theTransport);
  virtual void apply(const StartingPoint<Measureable> *theStartingPoint);
  virtual void apply(const FitAction<Measureable> *theFitterAction);

};

#include "CDFFitting/inc/LinearFitter.icc"

#endif /* CDFFitting_LinearFitter_hh */



