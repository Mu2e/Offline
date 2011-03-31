#ifndef FITTER_HH_
#define FITTER_HH_

#ifdef DEFECT_OLD_STL_HEADERS
#include <vector>
#else
#include <vector>
#endif

#ifndef DEFECT_NO_NAMESPACES
#ifndef DEFECT_NO_STDLIB_NAMESPACES
using std::vector;
#endif
#endif

template <class Measureable> class FitAction;
template <class Measureable> class Measurement;
template <class Measureable> class Constraint;
template <class Measureable> class Transport;
template <class Measureable> class Scatter;
template <class Measureable> class StartingPoint;
template <class Measureable> class Fit;
class Residual;

//	This is a very general interface for a fitter.  It
//	basically provides the overhead for adding various fit
//	actions (measuerment, constraint, transport, scatter, as
//	well as user fit interventions) to the fit, for doing
//	the fit, and for retrieving the result.

template <class Measureable>
class Fitter 
{

friend class FitAction<Measureable>;  
friend class Measurement<Measureable>;
friend class StartingPoint<Measureable>;
friend class Transport<Measureable>;
friend class Scatter<Measureable>;
friend class Constraint<Measureable>;

public:

  Fitter();
  
  Fitter(const Fitter<Measureable> &right);
  
  virtual ~Fitter();
  
  const Fitter<Measureable> & operator=(const Fitter<Measureable> &right);
  
  
  //	Adds a fit action to the fitter.  The order of the fit
  //	action is in general important, but whether it really
  //	makes any difference or not depends on the fitter
  //	subclass.
  void addFitAction(const FitAction<Measureable> *action);
  
  //	Gets the number of fit actions currenty attributed to
  //	the fitter.  They may or not be incorporated into the
  //	fit.
  unsigned int getNumFitActions() const;
  
  //	Abstract function.  Subclasses are expected to implement
  //	this function to provide the number of fit actions they
  //	have actually incorporated into the fit.
  virtual unsigned int getNumAppliedFitActions() const = 0;
  
  //	Returns a pointer to the ith fit action.
  const FitAction<Measureable> * getFitAction(unsigned int index) const;
  
  //	Sets the reference "measurable"; all calculations of
  //	residuals, derivates & cetera are w. r. t. this
  //	reference measurable.
  void setReference(Measureable *reference);
  
  //	Returns the reference measureable
  const Measureable * getReference() const;
  
  //	Abstract function for performing a fit.  The sequence of
  //	operations that will be carried out when this function
  //	is invoked is determined by the subclass.
  virtual void fit() = 0;
  
  
  //	An abstract function returning the "fit".  Provided by
  //	the subclasses.
  virtual const Fit<Measureable> * getFit() const = 0;
  
  //	An abstract function that returns the residual for a
  //	particular measurement.  The implementation is provided
  //	in the subclasses.
  virtual const Residual * newResidual(const Measurement<Measureable> *measurement) = 0;
  
  // Resets the fit
  virtual void resetFit()=0;

  // Failed flags:
  void setFailed(bool flag) {_failed=flag;}
  bool failed() const {return _failed;}

private:  

  vector<const FitAction<Measureable> *> _fitActionSet;
  Measureable *_reference;
  bool         _failed;

protected:

  // A Non-const version of getReference for subclasses
  Measureable * getReference();


protected:

  //	The standard set of fit actions is invoked using
  //	double-dispatch.  Subclasses of fitter handle the
  //	incorporation of specific fit actions.  In case the
  //	subclass doesn't handle the fit action (and if the fit
  //	action doesn't apply itself directly), this method
  //	simply prints out a hopefully not too bothersome warning
  //	message.
  virtual void apply(const Measurement<Measureable> *theMeasurement);  
  virtual void apply(const Constraint<Measureable> *theConstraint);  
  virtual void apply(const Scatter<Measureable> *theScatter);  
  virtual void apply(const Transport<Measureable> *theTransport);  
  virtual void apply(const StartingPoint<Measureable> *theStartingPoint);
  virtual void apply(const FitAction<Measureable> *theFitterAction);


};
#include "CDFFitting/inc/Fitter.icc"

#endif //FITTER_HH_

