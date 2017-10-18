#ifndef CRYEventGenerator_CosmicCRY_hh
#define CRYEventGenerator_CosmicCRY_hh

// Cosmic rays generator using CRY

#include <vector>

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoissonQ.h"

class TH1D;
class TH2D;
class TH1I;
class TTree;
namespace art{
  class Run;
}

class CRYSetup;
class CRYGenerator;
class GenParticleCollection;

namespace mu2e {

  class CosmicCRY{

  public:
    CosmicCRY( art::Run& run, const SimpleConfig& config );
    virtual ~CosmicCRY();

    virtual void generate( GenParticleCollection&  );

  private:

    int  _verbose;
    bool _doHistograms;
    bool _saveTree;
    TH2D *_hStartXZ;
    TH1D *_hStartY;
    // TH1D *_hStartPlane;
    TH1D *_hStartE;
    TH1D *_hStartTheta;
    TH1D *_hStartPhi;
    TH1D *_hPtot;
    TH1D *_hPyOverPtot;
    TH1D *_hTime;
    TH1D *_hParticleType;
    TH1D *_hDensityOverR;
    TH1D *_hNegMuKE;
    TH1D *_hPosMuKE;
    TH2D *_hPtypeKE;
    TH1D *_hnSecondaries;
    TTree *_tCryPrimary;
    TTree *_tCrySecondaries;

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

    enum RefPointChoice {UNDEFINED, TRACKER, EXTMONFNAL, CALO, CUSTOMIZED};
    enum DirectionChoice {ALL, POSITIVE_X, NEGATIVE_X, POSITIVE_Z, NEGATIVE_Z,
      PHI_RANGE};
    RefPointChoice    _refPointChoice;
    DirectionChoice   _directionChoice;
    CLHEP::Hep3Vector _cosmicReferencePointInMu2e;
    bool _vertical;
    bool _dontProjectToSurface;

    int _evtId0;
    int _pdgId0;
    double _ke0;
    int _nSecondaries;

    static const int _maxNSecondaries = 300;
    int _pdgId1[_maxNSecondaries];
    double _x1[_maxNSecondaries];
    double _y1[_maxNSecondaries];
    double _z1[_maxNSecondaries];
    double _t1[_maxNSecondaries];
    double _ke1[_maxNSecondaries];
    double _px1[_maxNSecondaries];
    double _py1[_maxNSecondaries];
    double _pz1[_maxNSecondaries];
    double _ptot1[_maxNSecondaries];
    double _theta1[_maxNSecondaries];
    double _phi1[_maxNSecondaries];

    void makeTrees();
    void createSetupString();
  };  // CosmicCRY

}
#endif /* end of include guard: CRYEventGenerator_CosmicCRY_hh */
