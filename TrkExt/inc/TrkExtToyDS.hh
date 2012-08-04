//
// ToyDS description for TrkExt
//
//  $Id: TrkExtToyDS.hh,v 1.1 2012/08/04 00:22:10 mjlee Exp $
//  $Author: mjlee $
//  $Date: 2012/08/04 00:22:10 $
//
//  Original author MyeongJae Lee
//
//
#ifndef TrkExtToyDS_HH
#define TrkExtToyDS_HH

#include <string>
#include <vector>
#include "CLHEP/Vector/ThreeVector.h"
#include "TrkExt/inc/TrkExtShape.hh"
#include "TrkExt/inc/TrkExtMaterial.hh"

namespace mu2e {


  class TrkExtToyDS : public TrkExtShape, public TrkExtMaterial 
  {

  public:
    TrkExtToyDS() ;
    ~TrkExtToyDS() { }

    void initialize () ;
    bool contains (CLHEP::Hep3Vector& p) ;

  private:
    std::string name;
    double rin, rout, zmin, zmax;

  };



} // end namespace mu2e


#endif
