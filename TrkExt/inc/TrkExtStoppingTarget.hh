//
// StoppingTarget description for TrkExt
//
//
//  Original author MyeongJae Lee
//
//
#ifndef TrkExtStoppingTarget_HH
#define TrkExtStoppingTarget_HH

#include <string>
#include <vector>
#include "CLHEP/Vector/ThreeVector.h"
#include "TrkExt/inc/TrkExtShape.hh"
#include "TrkExt/inc/TrkExtMaterial.hh"

namespace mu2e {


  class foil_data_type {
  public:
    foil_data_type (double x1, double x2, double x3, double x4) :
      rout(x1),
      z0(x2),
      zc(x3),
      z1(x4) 
    {}
    double rout, z0, zc, z1;
  };



  class TrkExtStoppingTarget : public TrkExtShape, public TrkExtMaterial 
  {

  public:
    TrkExtStoppingTarget() ;
    ~TrkExtStoppingTarget() { }

    void initialize () ;
    bool contains (CLHEP::Hep3Vector& p) ;

  private:
    std::vector<foil_data_type> foil;
    int nfoil;
    std::string name;
    double rmax, zmin, zmax;

  };



} // end namespace mu2e


#endif
