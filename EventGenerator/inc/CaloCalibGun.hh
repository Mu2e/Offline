#ifndef EventGenerator_CaloCalibGun_hh
#define EventGenerator_CaloCalibGun_hh
//
// Generate some number of DIO electrons.
//
// $Id: CaloCalibGun.hh,v 1.34 2013/07/22 18:57:42 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2013/07/22 18:57:42 $
//
//
// ====================================================================
//
// IMPORTANT NOTE:
//
//    _ehi MUST BE initialized before any of the CLHEP::Rand* variables
//
// ====================================================================

// C++ includes
#include <memory>

// Framework includes
#include "art/Framework/Principal/Run.h"

// Mu2e includes
#include "EventGenerator/inc/GeneratorBase.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"

// CLHEP includes
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Random/RandFlat.h"

// ROOT includes
#include "TTree.h"

class TTree;

namespace art {
  class Run;
}


namespace mu2e {

  class SimpleConfig;
  class DetectorSystem;

  class CaloCalibGun: public GeneratorBase {
  public:

    CaloCalibGun(CLHEP::HepRandomEngine& engine, art::Run& run, const SimpleConfig& config);
    virtual ~CaloCalibGun();

    virtual void generate( GenParticleCollection& );

  private:

    double _mean;
    double _energy;
    double _cosmin;
    double _cosmax;
    double _phimin;
    double _phimax;
    double _tmin;
    double _tmax;

    CLHEP::RandFlat     _randFlat;
    CLHEP::RandPoissonQ _randPoissonQ;
    RandomUnitSphere    _randomUnitSphere;

    const DetectorSystem  *_detSys;
    const DiskCalorimeter *_cal;
    int                    _nPipes;
    double                 _pipeRadius;
    std::vector<double>    _pipeTorRadius;
    std::vector<double>    _randomRad;
    CLHEP::Hep3Vector      _zPipeCenter;


    bool _doHistograms;
    TTree* _Ntupe;
    float _genErg;
    float _genTime;
    float _genCos;
    //float _genPhi;
    //float _genRad;
    float _genPosX;
    float _genPosY;
    float _genPosZ;
  };


} // end namespace mu2e,

#endif /* EventGenerator_CaloCalibGun_hh */
