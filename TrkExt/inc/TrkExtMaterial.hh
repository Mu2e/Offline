//
// Material description for TrkExt
//
//
//  Original author MyeongJae Lee
//
//
#ifndef TrkExtMaterial_HH
#define TrkExtMaterial_HH

#include <string>
#include "CLHEP/Vector/ThreeVector.h"


namespace mu2e {


  class TrkExtMaterial {

  public:
    TrkExtMaterial() ;
    TrkExtMaterial(std::string n) ;
    ~TrkExtMaterial() { }

    double meanEnergyLoss (const CLHEP::Hep3Vector& p, double ds) ;
    double mostProbableEnergyLoss (const CLHEP::Hep3Vector& p, double ds) ;
    double scatteringAngle (const CLHEP::Hep3Vector& p, double ds) ;
    std::string name () { return _name; }

  private:
    
    static constexpr double _mec22 = 0.510998910 * 0.510998910;
    static constexpr double _mec2  = 0.510998910;

    std::string _name;
    double _dpmp[6];
    double _edpmp[6];
    double _dpsp[3];
    double _edpsp[3];
    double _thpar[3];
    double _rho;

    enum TrkExtMaterialList {
      Undefined = -1,
      Vac = 0,
      Al = 1,
      PE = 2
    } _matid;

  };



} // end namespace mu2e


#endif
