#ifndef WLSPrimaryGeneratorAction_h
#define WLSPrimaryGeneratorAction_h 1

#include "globals.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4PhysicsOrderedFreeVector.hh"
#include "G4ThreeVector.hh"

class G4ParticleGun;
class G4Event;

class WLSPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:

    WLSPrimaryGeneratorAction(int mode);
    ~WLSPrimaryGeneratorAction();

    void GeneratePrimaries(G4Event*);
    void BuildEmissionSpectrum();
    G4ThreeVector GetOptPhotonStartPoint();
    void SetOptPhotonPolar();
    void SetOptPhotonEnergy(int spectrum);
    void SetBins(int binx, int biny, int binz);

  private:

    CLHEP::HepRandomEngine*           _randomEngine;
    G4ParticleGun*             _particleGun;

    int                        _mode;
    G4PhysicsOrderedFreeVector _emissionIntegral[2];
    double                     _yieldRatio;
    double                     _cerenkovEnergyMin, _cerenkovEnergyMax;
    bool                       _first;
    int                        _binx, _biny, _binz;
    bool                       _hasBins;
};

#endif
