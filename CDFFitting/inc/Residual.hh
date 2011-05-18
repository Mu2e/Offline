
#ifndef CDFFitting_Residual_hh
#define CDFFitting_Residual_hh

#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"
//	This class encapsulates the residual of a measurement,
//	and its covariance matrix, and also handles requests for
//	the dimensionality of these two objects, and for the
//	contribution to Chi-Squared that the residual makes.

class Residual {

public:

  inline Residual(const Residual &right);

  inline Residual(const HepVector &magnitude, const HepSymMatrix &error);

  inline ~Residual();

  inline const Residual & operator=(const Residual &right);


  inline const HepVector & getDisplacement() const;

  inline const HepSymMatrix & getErrorMatrix() const;

  inline unsigned int getDimensionality() const;

  inline double getDeltaChiSquared() const;

  inline Residual();


private:

  HepVector     _magnitude;
  HepSymMatrix  _error;

};

#include "CDFFitting/inc/Residual.icc"

#endif /* CDFFitting_Residual_hh */


