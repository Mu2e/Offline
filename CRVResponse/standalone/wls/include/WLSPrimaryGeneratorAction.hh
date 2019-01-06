#ifndef WLSPrimaryGeneratorAction_h
#define WLSPrimaryGeneratorAction_h 1

#include "globals.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4PhysicsOrderedFreeVector.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4ThreeVector.hh"
#include "WLSSteppingAction.hh"

class G4ParticleGun;
class G4Event;

class WLSPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
    WLSPrimaryGeneratorAction();

  public:

    WLSPrimaryGeneratorAction(WLSSteppingAction::simulationMode mode, 
                              int numberOfPhotons=-1, int simType=-1, int startBin=-1, bool verbose=false, double posY=0.0, double posZ=0.0);
    ~WLSPrimaryGeneratorAction();                                                                   //posY and posZ only used for mode 0 and 1

    void BuildEmissionSpectrum();
    bool SetNextBins();
    int  GeneratePhotonsInScintillator(G4Event *anEvent, int generatedPhotons);
    int  GenerateCerenkovPhotonsInFiber(G4Event *anEvent, int generatedPhotons);
    void GeneratePrimaries(G4Event*);

  private:

    CLHEP::HepRandomEngine*    _randomEngine;
    G4ParticleGun*             _particleGun;

    WLSSteppingAction::simulationMode _mode;

    int                        _numberOfPhotons, _simType, _currentBin;
    bool                       _verbose;
    G4PhysicsOrderedFreeVector _emissionIntegral;
    G4MaterialPropertyVector  *_rindexScintillator, *_rindexFiber;
    double                     _cerenkovEnergyMinScintillator, _cerenkovEnergyMaxScintillator;
    double                     _cerenkovEnergyMinFiber, _cerenkovEnergyMaxFiber;
    double                     _maxRIndexScintillator, _maxRIndexFiber;
    double                     _scintillationRiseTime, _scintillationDecayTime;
    bool                       _first;

    double                     _minBinX, _minBinY, _minBinZ, _minBinBeta, _minBinTheta, _minBinPhi, _minBinR;
    double                     _maxBinX, _maxBinY, _maxBinZ, _maxBinBeta, _maxBinTheta, _maxBinPhi, _maxBinR;

    double                     _posY, _posZ;  //only used for mode 0 and 1
};

#endif
