#include "G4ios.hh"
#include "G4Event.hh"

#include "G4ParticleGun.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

#include "G4PhysicsTable.hh"

#include "Randomize.hh"

#include "WLSPrimaryGeneratorAction.hh"
#include "WLSEventAction.hh"
#include "WLSSteppingAction.hh"
#include "WLSDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"

#include "TH3D.h"

#include "CLHEP/Random/Random.h"
#include "CLHEP/Random/RandGaussQ.h"

#include <stdexcept>

WLSPrimaryGeneratorAction::WLSPrimaryGeneratorAction(int mode) : _mode(mode), _first(true), _hasBins(false)
{
  _randomEngine = CLHEP::HepRandom::getTheEngine();
  _particleGun  = new G4ParticleGun(1);
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  if(_mode==-1)
  {
    G4ParticleDefinition* particle = particleTable->FindParticle("opticalphoton");
    _particleGun->SetParticleDefinition(particle);
  }
  if(_mode==0)
  {
    G4ParticleDefinition* particle = particleTable->FindParticle("proton");
//    G4ParticleDefinition* particle = particleTable->FindParticle("neutron");
//    G4ParticleDefinition* particle = particleTable->FindParticle("gamma");
//    G4ParticleDefinition* particle = particleTable->FindParticle("mu-");
    _particleGun->SetParticleDefinition(particle);
    _particleGun->SetParticleMomentumDirection(G4ThreeVector(1., 0., 0.));
//    _particleGun->SetParticleMomentumDirection(G4ThreeVector(1.0, -0.9, -5.0));
    _particleGun->SetParticleEnergy(120.*GeV);
    _particleGun->SetParticleTime(0);
  }  
}

WLSPrimaryGeneratorAction::~WLSPrimaryGeneratorAction()
{
  delete _particleGun;
}

void WLSPrimaryGeneratorAction::BuildEmissionSpectrum()
{
  G4Material* material = G4Material::GetMaterial("Polystyrene",true);
  G4MaterialPropertiesTable* materialPropertiesTable = material->GetMaterialPropertiesTable();

//scintillation
  _yieldRatio = materialPropertiesTable->GetConstProperty("YIELDRATIO");

  std::string tag[2]={"FASTCOMPONENT","SLOWCOMPONENT"};
  for(int iTag=0; iTag<2; iTag++)
  {
    const G4MaterialPropertyVector &component = *materialPropertiesTable->GetProperty(tag[iTag].c_str());

    double currentEnergy = component.Energy(0);
    double currentDifferentialProb = component[0];
    double currentIntegratedProb = 0.0;
    _emissionIntegral[iTag].InsertValues(currentEnergy, currentIntegratedProb);

    for(size_t i=1; i<component.GetVectorLength(); i++)
    {
      double prevEnergy  = currentEnergy;
      double prevDifferentialProb = currentDifferentialProb;
      double prevIntegratedProb = currentIntegratedProb;

      currentEnergy = component.Energy(i);
      currentDifferentialProb = component[i];
      currentIntegratedProb = 0.5 * (prevDifferentialProb + currentDifferentialProb);  //TODO: Why?
      currentIntegratedProb = prevIntegratedProb + (currentEnergy - prevEnergy) * currentIntegratedProb;
      _emissionIntegral[iTag].InsertValues(currentEnergy, currentIntegratedProb);
    }
  }

//Cerenkov
  G4MaterialPropertyVector *rindex = materialPropertiesTable->GetProperty("RINDEX");
  _cerenkovEnergyMin = rindex->GetMinLowEdgeEnergy();
  _cerenkovEnergyMax = rindex->GetMaxLowEdgeEnergy();
}

void WLSPrimaryGeneratorAction::SetBins(int binx, int biny, int binz)
{
  _hasBins=true;
  _binx=binx;
  _biny=biny;
  _binz=binz;
}

G4ThreeVector WLSPrimaryGeneratorAction::GetOptPhotonStartPoint()
{
  if(!_hasBins) throw std::logic_error("Bins need to be set for mode -1");

  G4ThreeVector point = WLSEventAction::Instance()->GetHistBinCenter(_binx, _biny, _binz);

  return(point);
}

void WLSPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
//GEANT simulation of opt photons / propagate / hits to build maps
  if(_mode==-1)
  {
    if(_first)
    {
      _first=false;

      G4Material* material = G4Material::GetMaterial("PMMA",true);
      G4MaterialPropertiesTable* materialPropertiesTable = material->GetMaterialPropertiesTable();
      materialPropertiesTable->AddConstProperty("WLSTIMECONSTANT", 0.0);

      BuildEmissionSpectrum();
    }

//    int generatedPhotons=100000;  //takes about 5 min
int generatedPhotons=100;  //takes about 5 min
    G4ThreeVector startCenter = GetOptPhotonStartPoint();
    double binWidthX = WLSEventAction::Instance()->GetHistBinWidthX(_binx);
    double binWidthY = WLSEventAction::Instance()->GetHistBinWidthY(_biny);
    double binWidthZ = WLSEventAction::Instance()->GetHistBinWidthZ(_binz);

    double fiberRadius = WLSDetectorConstruction::Instance()->GetClad2Radius();
    double fiberSeparation = WLSDetectorConstruction::Instance()->GetFiberSeparation();
    if(anEvent->GetEventID()>=2) //Cerenkov inside of fiber
    {
//      generatedPhotons=10000;
generatedPhotons=100;
      startCenter.setX(0);
      startCenter.setY(anEvent->GetEventID()==2?-fiberSeparation/2.0:fiberSeparation/2.0); //two fibers
      binWidthX=fiberRadius*2.0;
      binWidthY=fiberRadius*2.0;
    }

    WLSEventAction::Instance()->SetOptPhotonStart(startCenter);

    G4Navigator* navigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
    G4VPhysicalVolume *scintillator = WLSDetectorConstruction::Instance()->GetScintillatorVolume();

    int actualNumberOfGeneratedPhotons=0;
    for(int i=0; i<generatedPhotons; i++)
    {
      double dx=(G4UniformRand()-0.5)*binWidthX;
      double dy=(G4UniformRand()-0.5)*binWidthY;
      double dz=(G4UniformRand()-0.5)*binWidthZ;
      G4ThreeVector delta(dx,dy,dz);
      G4ThreeVector start = startCenter + delta;
      _particleGun->SetParticlePosition(start);

      if(anEvent->GetEventID()<2) //scintillation and Cerenkov only in scintillator
      {
        if(navigator->LocateGlobalPointAndSetup(start,NULL,false,true)!=scintillator) continue;
      }
      else  //Cerenkov in fiber
      {
        if(sqrt(dx*dx+dy*dy)>fiberRadius) continue;
      }
      actualNumberOfGeneratedPhotons++;

      double theta = acos(1.0-2.0*G4UniformRand());      //theta of emitted photon
      double phi   = CLHEP::twopi*G4UniformRand();       //phi of the emitted photon
      G4ThreeVector photonMomentum;
      photonMomentum.setRThetaPhi(1.0, theta, phi);
      _particleGun->SetParticleMomentumDirection(photonMomentum);

      SetOptPhotonPolar();
      SetOptPhotonEnergy(anEvent->GetEventID());  //energy for either the scintillation or the cerenkov spektrum
      _particleGun->SetParticleTime(0);
      _particleGun->GeneratePrimaryVertex(anEvent); 
    }
    WLSEventAction::Instance()->SetGeneratedOptPhotons(actualNumberOfGeneratedPhotons);
  }

//GEANT simulation of proton / energy dep / opt photons / propagate / hits
  if(_mode==0)
  {
    double beamsize=1.0*mm;
    double x0 = -1.5*cm;
//    double y0 = CLHEP::RandGaussQ::shoot(_randomEngine,0.0*cm,beamsize);
    double y0 = CLHEP::RandGaussQ::shoot(_randomEngine,1.3*cm,beamsize);
//    double y0 = CLHEP::RandGaussQ::shoot(_randomEngine,1.8*cm,beamsize); 
    double z0 = CLHEP::RandGaussQ::shoot(_randomEngine,20.0*cm,beamsize);
  
    _particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
    _particleGun->GeneratePrimaryVertex(anEvent);
  }
}

void WLSPrimaryGeneratorAction::SetOptPhotonPolar()
{
  G4ThreeVector normal (0., 1., 0.);
  G4ThreeVector kphoton = _particleGun->GetParticleMomentumDirection();
  G4ThreeVector product = normal.cross(kphoton);
  G4double modul2       = product*product;

  G4ThreeVector e_perpend (0., 0., 1.);
  if (modul2 > 0.) e_perpend = (1./std::sqrt(modul2))*product;
  G4ThreeVector e_paralle    = e_perpend.cross(kphoton);

  double angle   = CLHEP::twopi*G4UniformRand(); 
  G4ThreeVector polar = std::cos(angle)*e_paralle + std::sin(angle)*e_perpend;
  _particleGun->SetParticlePolarization(polar);
}

void WLSPrimaryGeneratorAction::SetOptPhotonEnergy(int table)
{
  double energy=0;
  if(table==0) //scintillation
  {
    int iTag=1;  //slow
    if(G4UniformRand()<=_yieldRatio) iTag=0; //fast

    double integratedProbMax = _emissionIntegral[iTag].GetMaxValue();
    double integratedProb = G4UniformRand()*integratedProbMax;
    energy = _emissionIntegral[iTag].GetEnergy(integratedProb);
  }
  else //cerenkov
  {
    energy = _cerenkovEnergyMin + G4UniformRand()*(_cerenkovEnergyMax-_cerenkovEnergyMin);
  }

  _particleGun->SetParticleEnergy(energy);
}
