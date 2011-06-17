//--------------------------------------------------------------------------
// Name:
//   DetStrawHitElem: Dummy class to represent a straw element.
//        Copyright (C) 2008    Lawrence Berkeley Laboratory
// Author List:
//      Dave Brown 27 May 2011
//------------------------------------------------------------------------
#ifndef DetStrawHitElem_HH
#define DetStrawHitElem_HH

#include "DetectorModel/DetElem.hh"
#include "DetectorModel/DetType.hh"
#include "MatEnv/MatDBInfo.hh"

namespace mu2e {
  class TrkStrawHit;
  enum elemType {wall, gas };
// define the type for this class first  
  class DetStrawHitType : public DetType {
  public:
    DetStrawHitType(elemType etype, MatDBInfo const* matdbinfo) {
      if(etype == wall)
        _mat = matdbinfo->findDetMaterial("straw-wall");
      else if(etype == gas)
        _mat = matdbinfo->findDetMaterial("straw-gas"); }
    virtual ~DetStrawHitType() {}
// DetType interface
    virtual bool physicalMaterial(const TypeCoord*) const { return true; }
    virtual const DetMaterial& material(const TypeCoord*) const { return *_mat; }
  private:
    const DetMaterial* _mat;
  };
  
// element class  
  class DetStrawHitElem : public DetElem {
  public:
// construct from a TrkStrawHit
    DetStrawHitElem(elemType etype, MatDBInfo const* matdbinfo) :
    DetElem(new DetStrawHitType(etype,matdbinfo),"DetStrawHitElem",_ielem++) {}
    virtual ~DetStrawHitElem() { delete const_cast<DetType*>(detectorType()); }
// DetElem interface
    virtual int intersect(const Trajectory*,DetIntersection&) const{ return 0; }
  protected:
    virtual HepPoint coordToPoint( const TypeCoord* aCoord ) const { return HepPoint(0.0,0.0,0.0); }
  private:
    static unsigned _ielem;
  };
}

#endif
