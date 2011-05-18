#ifndef CDFFitting_Constraint_hh
#define CDFFitting_Constraint_hh
#include "CDFFitting/inc/FitAction.hh"


//	This class is an Abstract Base Class representing
//	constraints.  Constraints are handled through the use of
//	Lagrange Multipliers.  In this formalism the constraint
//	is expressed in terms of a matrix and a vector which are
//	to be provided by the subclass.

template <class Measureable>
class Constraint : public FitAction<Measureable> {

public:

  // Constructor
  Constraint();

  // Destructor
  virtual ~Constraint();

  // Get the matrix
  virtual HepMatrix getConstraintMatrix(const Measureable &  measurable) const = 0;

  // Get the vector
  virtual HepVector getConstraintVector(const Measureable &  measurable) const = 0;

  // Get the dimensionality of this constraint
  virtual unsigned int getDimensionality() const = 0;

  // Half of a double dispatch pair:
  virtual void applyYourselfTo(Fitter<Measureable> *theFitter) const;
};

#include "CDFFitting/inc/Constraint.icc"




#endif /* CDFFitting_Constraint_hh */


