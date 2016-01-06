//--------------------------------------------------------------------------
// Name:
//   DetStrawType:  Class to define DetectorModel type for a straw
//        Copyright (C) 2011    Lawrence Berkeley Laboratory
// Author List:
//      Dave Brown 9 Sept. 2011
//------------------------------------------------------------------------
#ifndef DetStrawType_hh
#define DetStrawType_hh

#include "BTrk/DetectorModel/DetType.hh"
#include "BTrk/DetectorModel/DetMaterial.hh"

namespace mu2e {
class DetStrawType : public DetType {
  public:
    DetStrawType(const DetMaterial* gasmat, const DetMaterial* wallmat, const DetMaterial* wiremat,
      double offset, double tol, double rfrac) :
      _gasmat(gasmat), _wallmat(wallmat), _wiremat(wiremat), _offset(offset), _tol(tol), _rfrac(rfrac) {}
    virtual ~DetStrawType() {}
// Trivial implementation of DetType interface; these are not used
    virtual bool physicalMaterial(const TypeCoord*) const { return true; }
    virtual const DetMaterial& material(const TypeCoord*) const { return *_gasmat; }
// specific functions
    const DetMaterial* gasMaterial() const { return _gasmat; }
    const DetMaterial* wallMaterial() const { return _wallmat; }
    const DetMaterial* wireMaterial() const { return _wiremat; }
    double offset() const { return _offset; }
    double tolerance() const { return _tol; }
    double maxRadiusFraction() const { return _rfrac; }
  private:
    const DetMaterial* _gasmat;
    const DetMaterial* _wallmat;
    const DetMaterial* _wiremat;
    double _offset;
    double _tol;
    double _rfrac;
  };
}  
#endif
