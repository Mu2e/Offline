#ifndef CDFFitting_FitAction_hh
#define CDFFitting_FitAction_hh
#include "CDFFitting/inc/Fitter.hh"
#include "CDFFitting/inc/FitAction.hh"

//	This is an abstract base class for Fit Actions.
//	Subclasses of FitAction that are automatically provided
//	are Measurement (probably the most important),
//	Constraint, Transport, and Scatter.  These  fit actions
//	are referred back to the fitter for (double dispatch)
//	for appropriate action.
//
//	Fit actions can be extended with "User Fit Actions", that
//	are invoked at the appropriate point during the fitting
//	procedure.

template <class Measureable>
class FitAction  {
  
public:

  // Constructor
  FitAction();
  
  // Destructor
  virtual ~FitAction();

  //	This method is part of a double-dispatch pair.  The
  //	fitter requests that a fit action be incorporated.
  //	particular kind of fitter is then requested to
  //	incorporate the particular kind of fit.
  virtual void applyYourselfTo(Fitter<Measureable> *theFitter) const=0;
  
};

#include "CDFFitting/inc/FitAction.icc"

#endif /* CDFFitting_FitAction_hh */
