#ifndef CDFFitting_ProgressiveFitter_hh
#define CDFFitting_ProgressiveFitter_hh
#include "CDFFitting/inc/Fitter.hh"
#include "CDFFitting/inc/Fit.hh"
#include "CDFFitting/inc/Residual.hh"
#include "CDFFitting/inc/ProgressiveFitter.hh"
#include "CDFFitting/inc/Scatter.hh"
#include "CDFFitting/inc/Transport.hh"
#include "CDFFitting/inc/Constraint.hh"
#include "CDFFitting/inc/StartingPoint.hh"
#include "CLHEP/Matrix/DiagMatrix.h"

//      This class represents a progressive fitter.  A
//      progressive fitter incorporates information
//      sequentially, but in the presence of "system noise"
//      cannot provide accurate estimates of the state vector at
//      any position except the latest.  For cases where this is
//      important it is better to use the Kalman Filter.
//
//      The advantage of the Progressive Fitter is a possiblly
//      faster execution time.
template <class Measureable>
class ProgressiveFitter : public Fitter<Measureable> {


public:

  // Constructor
  ProgressiveFitter();

  // Copy Constructor
  ProgressiveFitter(const ProgressiveFitter<Measureable> &right);

  // Destructor
  virtual ~ProgressiveFitter();

  // Assignment
  const ProgressiveFitter<Measureable> & operator=(const ProgressiveFitter<Measureable> &right);

  // Returns the number of Fit Actions applied so far by the fitter.
  virtual unsigned int getNumAppliedFitActions() const;

  // Apply all of the fit actions sequentially, first through last.
  virtual void fit();

  // Incorporate the next unincorporated fit action; return
  // true if an action remained to be incorporated, false
  // otherwise.
  virtual bool incorporateNextAction();

  // Return the fit.  Return nullptr if no fit has been
  // performed.
  virtual const Fit<Measureable> * getFit() const;

  // return the Residual of the ith measurement.  In the
  // presence of system noise, this will be an inaccurate
  // estimate.
  virtual const Residual * newResidual(const Measurement<Measureable> *measurement) ;

  // Reset the fit
  virtual void resetFit();


protected:

  HepVector _estimate;
  HepSymMatrix _covariance;
  double _chiSquared;
  int _degreesOfFreedom;
  int _numAppliedFitActions;
  Fit<Measureable> *   _fit;


  // Some workspace:  These are not part of the state and would
  // in fact be *local* variables (in the Progressive Fitter apply
  // measurement method) except:  1) let's allocate them once and
  // for all, and 2) certain subclasses (kalman) may need to peak
  // at these variables.   In order to stow them away for itself.
  HepVector      _residual;
  HepSymMatrix   _residualCov;
  HepSymMatrix   _errorMatrix;
  HepMatrix      _derivativeMatrix;
  HepVector      _displacementVector;


protected:

  virtual void apply(const Measurement<Measureable> *theFitterAction);
  virtual void apply(const Constraint<Measureable> *theFitterAction);
  virtual void apply(const Scatter<Measureable> *theFitterAction);
  virtual void apply(const Transport<Measureable> *theFitterAction);
  virtual void apply(const StartingPoint<Measureable> *theFitterAction);
  virtual void apply(const FitAction<Measureable> *theFitterAction);

};

#include "CDFFitting/inc/ProgressiveFitter.icc"

#endif /* CDFFitting_ProgressiveFitter_hh */


