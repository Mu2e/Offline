#ifndef CDFFitting_Transport_hh
#define CDFFitting_Transport_hh
#include "CDFFitting/inc/FitAction.hh"
#include "CLHEP/Matrix/Matrix.h"

//	This class is an Abstract Base Class that specifies the
//	interface for "transports", e.g, in the case of track
//	fitting,
//	energy loss, magnetic field inhomogeneities.  In general
//	a transport produces a deterministic shift of the state
//	vector.
template <class Measureable>
class Transport : public FitAction<Measureable> {

public:

  // Constructor
  Transport();

  // Destructor
  virtual ~Transport();

  // Transportation Matrix
  virtual HepMatrix getTransportMatrix(const Measureable &  measureable) const = 0;

  // Dimensionality
  virtual unsigned int getDimensionality() const = 0;

  // Half of a double dispatch pair:
  virtual void applyYourselfTo(Fitter<Measureable> *theFitter) const;
};
#include "CDFFitting/inc/Transport.icc"


#endif /* CDFFitting_Transport_hh */


