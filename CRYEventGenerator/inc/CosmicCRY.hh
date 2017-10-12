#ifndef EventGenerator_CosmicCRY_hh
#define EventGenerator_CosmicCRY_hh

// Cosmic rays generator using CRY

#include <vector>

#include "EventGenerator/inc/GeneratorBase.hh"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoissonQ.h"

class TH1D;
class TH2D;
class TTree;
namespace art{
  class Run;
}

class CRYSetup;
class CRYGenerator;

namespace mu2e {

  class CosmicCRY: public GeneratorBase {

  public:
    CosmicCRY( art::Run& run, const SimpleConfig& config );
    virtual ~CosmicCRY();

    virtual void generate( GenParticleCollection&  );

  private:

    int   _verbose;
    bool  _doHistograms;
    TH2D *_hStartXZ;
    TH1D *_hStartY;
    // TH1D *_hStartPlane;
    TH1D *_hStartE;
    TH1D *_hStartTheta;
    TH1D *_hStartPhi;
    TTree *_tCry;

    double _muEMin;   // min and max values of muon energy (MeV)
    double _muEMax;
    double _muCosThMin; // min and max zenith angles
    double _muCosThMax;
    double _muPhiMin;
    double _muPhiMax;

    // CRY output options
    bool _returnMuons;
    bool _returnNeutrons;
    bool _returnProtons;
    bool _returnGammas;
    bool _returnElectrons;
    bool _returnPions;
    bool _returnKaons;

    // CRY input parameters:
    // - the date that we want to simulate flux, month-day-year
    // - latitude of the detector site, 41.8 for Fermilab
    // - altitude in meter, CRY accepts 3 values: 0, 2100, 11300, default to 0
    // - sub box length
    int _month;
    int _day;
    int _year;
    double _latitude;
    int _altitude;
    double _subboxLength;

    std::string _setupString;
    std::string _cryDataPath;

    CRYSetup * _crySetup;
    CRYGenerator * _cryGen;

    int pdgId;
    double ke0;
    double px0;
    double py0;
    double pz0;
    double ptot0;
    double x0;
    double y0;
    double z0;
    double theta0;
    double phi0;
  };  // CosmicCRY

}
#endif /* end of include guard: EventGenerator_CosmicCRY_hh */
