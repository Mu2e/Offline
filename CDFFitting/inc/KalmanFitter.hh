#ifndef CDFFitting_KalmanFitter_hh
#define CDFFitting_KalmanFitter_hh
#include <vector>
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CDFFitting/inc/Fit.hh"
#include "CDFFitting/inc/ProgressiveFitter.hh"
#include "CDFFitting/inc/Residual.hh"
#include "CDFFitting/inc/KalmanFitter.hh"

//      This class represents a Kalman Filter fitter.  The
//      Kalman filter gets the residuals right, and gets the
//      right estimate of parameters at any point along the fit,
//      including at the beginning.
//
//      The Kalman buffers all information it needs to calculate
//      the so-called "smoothed estimate" at prior instants, and
//      actually calculates the estimate upon demand.
template <class Measureable>
class KalmanFitter : public ProgressiveFitter<Measureable> {

public:

  // constructor
  KalmanFitter();

  // copy constructor
  KalmanFitter(const KalmanFitter<Measureable> &right);

  // destructor
  virtual ~KalmanFitter();

  // assignment
  const KalmanFitter<Measureable> & operator=(const KalmanFitter<Measureable> &right);

  // get the fit as of the most "recent" innovation
  virtual const Fit<Measureable> * getFit() const;

  // get the fit at the "time" of a particular fit action.
  virtual const Fit<Measureable> * getFit(const FitAction<Measureable> *action) ;

  // get the residual, in this case, the "smoothed residual" (correct residual).
  virtual const Residual * newResidual(const Measurement<Measureable> *measurement);

  // reset the fit
  virtual void resetFit();

private:

  vector <HepVector>                                    _FilteredEstimate;
  vector <HepSymMatrix>                                _FilteredCovariance;
  vector <HepVector>                                    _FilteredResidual;
  vector <HepSymMatrix>                                _FilteredResidualCovariance;
  vector <HepMatrix>                                    _SmootherGain;
  vector <HepMatrix>                                    _DerivativeMatrix;
  vector <const Measurement<Measureable> *>             _MeasurementArray;
  vector <HepVector>                                    _SmoothedEstimate;
  vector <HepSymMatrix>                                _SmoothedCovariance;
  vector <HepVector>                                    _SmoothedResidual;
  vector <HepSymMatrix>                                _SmoothedResidualCovariance;
  vector <Fit <Measureable> >                            _FitArray;
  vector <Residual>                                     _ResidualArray;
  unsigned int                                          _validDepth;

protected:

  virtual void                                          smooth(unsigned int depth=0);

protected:
  // apply a measurement action
  virtual void apply(const Constraint<Measureable> *theFitterAction);
  virtual void apply(const Measurement<Measureable> *theFitterAction);
  virtual void apply(const Scatter<Measureable> *theFitterAction);
  virtual void apply(const Transport<Measureable> *theFitterAction);
  virtual void apply(const StartingPoint<Measureable> *theFitterAction);
  virtual void apply(const FitAction<Measureable> *theFitterAction);

};

#include "CDFFitting/inc/KalmanFitter.icc"

#endif /* CDFFitting_KalmanFitter_hh */



