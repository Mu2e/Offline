//
//  Class recording the intersection (geometric and physical) of a KinKal fit with a
//  geometric surface.
//
#ifndef RecoDataProducts_KalIntersection_HH
#define RecoDataProducts_KalIntersection_HH
#include "KinKal/Geometry/InterData.hh"
#include "CLHEP/Matrix/Vector.h"
#include "Offline/DataProducts/inc/GenVector.hh"
namespace mu2e {
  struct KalIntersection {
    // record which system the intersection was in
    enum system{tracker=0,  calo, crv, target, ds, crvtest};
    KinKal::InterData data_; // main payload from KinKal intersection
    system sys_; // which system the intersection surface is in
    int elem_; // which element in the system
  };
}

#endif
