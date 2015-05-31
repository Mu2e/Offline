#ifndef WLSPrimaryGeneratorAction_h
#define WLSPrimaryGeneratorAction_h 1

#include "globals.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4PhysicsOrderedFreeVector.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4ThreeVector.hh"

class G4ParticleGun;
class G4Event;

class WLSPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
    WLSPrimaryGeneratorAction();

  public:

    WLSPrimaryGeneratorAction(int mode, int numberOfPhotons=-1, int simType=-1, int startBin=-1, bool verbose=false);
    ~WLSPrimaryGeneratorAction();

    void BuildEmissionSpectrum();
    bool SetNextBins();
    int  GeneratePhotonsInScintillator(G4Event *anEvent, int generatedPhotons);
    int  GenerateCerenkovPhotonsInFiber(G4Event *anEvent, int generatedPhotons);
    void GeneratePrimaries(G4Event*);

  private:

    CLHEP::HepRandomEngine*    _randomEngine;
    G4ParticleGun*             _particleGun;

    int                        _mode, _numberOfPhotons, _simType, _currentBin;
    bool                       _verbose;
    G4PhysicsOrderedFreeVector _emissionIntegral[2];
    G4MaterialPropertyVector*  _rindexFiber;
    double                     _cerenkovEnergyMinScintillator, _cerenkovEnergyMaxScintillator;
    double                     _cerenkovEnergyMinFiber, _cerenkovEnergyMaxFiber;
    double                     _maxRIndexScintillator, _maxRIndexFiber;
    double                     _yieldRatio;
    bool                       _first;

    double                     _minBinX, _minBinY, _minBinZ, _minBinBeta, _minBinTheta, _minBinPhi, _minBinR;
    double                     _maxBinX, _maxBinY, _maxBinZ, _maxBinBeta, _maxBinTheta, _maxBinPhi, _maxBinR;
};

#endif
