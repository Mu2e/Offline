//--------------------------------------------------------------------------
// Name:
//   DetStrawHitType:  Trivial class to define types for TrkStrawHits
//        Copyright (C) 2011    Lawrence Berkeley Laboratory
// Author List:
//      Dave Brown 9 Sept. 2011
//------------------------------------------------------------------------
#ifndef DetStrawHitType_hh
#define DetStrawHitType_hh

#include "BTrk/DetectorModel/DetType.hh"
#include "BTrk/MatEnv/MatDBInfo.hh"

namespace mu2e {
class DetStrawHitType : public DetType {
  public:
    DetStrawHitType(MatDBInfo const* matdbinfo,const char* strawwallmat) {
      _mat = matdbinfo->findDetMaterial(strawwallmat); }
    virtual ~DetStrawHitType() {}
// DetType interface
    virtual bool physicalMaterial(const TypeCoord*) const { return true; }
    virtual const DetMaterial& material(const TypeCoord*) const { return *_mat; }
  private:
    const DetMaterial* _mat;
  };
}  
#endif
