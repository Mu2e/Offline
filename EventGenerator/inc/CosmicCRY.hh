#ifndef CRYEventGenerator_CosmicCRY_hh
#define CRYEventGenerator_CosmicCRY_hh

// Cosmic rays generator using CRY

#include <vector>

// #include "CLHEP/Random/RandFlat.h"
// #include "CLHEP/Random/RandPoissonQ.h"
// #include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Random/RandEngine.h"

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
    CosmicCRY(art::Run& run, const SimpleConfig& config,
        CLHEP::HepRandomEngine& engine);
    virtual ~CosmicCRY();

    virtual void generate( GenParticleCollection&  );

  private:

    int  _verbose;
    bool _doHistograms;
    bool _saveTree;
    TH2D *_hXZ;
    TH1D *_hY;
    TH1D *_hE;
    TH1D *_hTheta;
    TH1D *_hPhi;
    TH1D *_hPtot;
    TH1D *_hPyOverPtot;
    TH1D *_hTime;
    TH1D *_hLiveTime;
    TH1D *_hNegMuKE;
    TH1D *_hPosMuKE;
    TH2D *_hPtypeKE;
    TH1D *_hNSecondaries;
    TH2D *_hSecondPtotVsPrimKE;
    TH2D *_hShowerRadiusVsPrimKE;
    TH2D *_hNSecondariesVsPrimKE;

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
    std::shared_ptr<CRYGenerator> _cryGen;

    double _refY0;
    std::string _refPointChoice;
    std::string _directionChoice;
    CLHEP::Hep3Vector _cosmicReferencePointInMu2e;
    bool _vertical;

    bool _projectToEnvelope;

    bool _geomInfoObtained;
    double _envXmin;
    double _envXmax;
    double _envYmin;
    double _envYmax;
    double _envZmin;
    double _envZmax;

    std::vector<CLHEP::Hep3Vector> _envIntersections;
    void calIntersections(CLHEP::Hep3Vector orig, CLHEP::Hep3Vector dir);
    bool pointInBox(double x, double y, double x0, double y0, double x1, double z1);

    int _evtId0;
    int _pdgId0;
    double _ke0;
    double _t0;
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
    double calMEC(std::vector<CLHEP::Hep2Vector> points); // minimum enclosing circle
  };  // CosmicCRY

}
#endif /* end of include guard: CRYEventGenerator_CosmicCRY_hh */
