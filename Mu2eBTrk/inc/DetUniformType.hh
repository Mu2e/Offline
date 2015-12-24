//--------------------------------------------------------------------------
// Name:
//   DetUniformType:  Trivial class to define types for uniform elements
//        Copyright (C) 2011    Lawrence Berkeley Laboratory
// Author List:
//      Dave Brown 9 Sept. 2011
//------------------------------------------------------------------------
#ifndef DetUniformType_hh
#define DetUniformType_hh

#include "BTrk/DetectorModel/DetType.hh"
#include "BTrk/DetectorModel/DetMaterial.hh"

namespace mu2e {
class DetUniformType : public DetType {
  public:
    DetUniformType(const DetMaterial* mat) : _mat(mat) {}
    virtual ~DetUniformType() {}
// Trivial implementation of DetType interface
    virtual bool physicalMaterial(const TypeCoord*) const { return true; }
    virtual const DetMaterial& material(const TypeCoord*) const { return *_mat; }
  private:
    const DetMaterial* _mat;
  };
}  
#endif
