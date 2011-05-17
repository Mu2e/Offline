#ifndef CDFFitting_StartingPoint_hh
#define CDFFitting_StartingPoint_hh
#include "CDFFitting/inc/FitAction.hh"
#include "CDFFitting/inc/StartingPoint.hh"


class HepVector;
class HepSymMatrix;
//	This class is an Abstract Base Class representing a
//	starting point, typically of use only for one of the
//	recursive fitters.  This class is rather like a
//	measurement but adds no degrees of freedom to the fit.

template <class Measureable>
class StartingPoint : public FitAction<Measureable> {

public:
  
  // constructor
  StartingPoint();

  // destructor
  virtual ~StartingPoint();
  
  // Gets the error matrix of the reference value at the start of the fit;  not
  // necessarily in Cartesian 3-space 
  virtual HepSymMatrix getErrorMatrix(const Measureable &measureable) const = 0;
  
  // Get the dimensionality of this constraint
  virtual unsigned int getDimensionality() const = 0;

  // Half of a double dispatch pair:
  virtual void applyYourselfTo(Fitter<Measureable> *theFitter) const;
};
#include "CDFFitting/inc/StartingPoint.icc"
#endif /* CDFFitting_StartingPoint_hh */


