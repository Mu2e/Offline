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

#include <stdexcept>

WLSPrimaryGeneratorAction::WLSPrimaryGeneratorAction(int mode, int numberOfPhotons, int simType, int startBin, bool verbose, double posY, double posZ) : 
                                                                                            _mode(mode), 
                                                                                            _numberOfPhotons(numberOfPhotons),
                                                                                            _simType(simType),
                                                                                            _currentBin(startBin-1), 
                                                                                            _verbose(verbose),
                                                                                            _first(true),
                                                                                            _posY(posY), _posZ(posZ)
{
  _randomEngine = CLHEP::HepRandom::getTheEngine();
  _particleGun  = new G4ParticleGun(1);
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  if(_mode==-1)
  {
    G4ParticleDefinition* particle = particleTable->FindParticle("opticalphoton");
    _particleGun->SetParticleDefinition(particle);
  }
  if(_mode==0 || _mode==1)
  {
//    G4ParticleDefinition* particle = particleTable->FindParticle("proton");
//    G4ParticleDefinition* particle = particleTable->FindParticle("neutron");
//    G4ParticleDefinition* particle = particleTable->FindParticle("gamma");
    G4ParticleDefinition* particle = particleTable->FindParticle("mu-"); 
//    G4ParticleDefinition* particle = particleTable->FindParticle("opticalphoton");  //WLS fiber test
    _particleGun->SetParticleDefinition(particle);
    _particleGun->SetParticleMomentumDirection(G4ThreeVector(1., 0., 0.));
    _particleGun->SetParticleEnergy(120.*GeV);
//    _particleGun->SetParticleEnergy(3.7*eV);  //WLS fiber test
    _particleGun->SetParticleTime(0);
  }  
}

WLSPrimaryGeneratorAction::~WLSPrimaryGeneratorAction()
{
  delete _particleGun;
}

void WLSPrimaryGeneratorAction::BuildEmissionSpectrum()
{
  G4Material* scintillator = G4Material::GetMaterial("Polystyrene",true);
  G4MaterialPropertiesTable* scintillatorPropertiesTable = scintillator->GetMaterialPropertiesTable();

//scintillation
  _yieldRatio = scintillatorPropertiesTable->GetConstProperty("YIELDRATIO");

  std::string tag[2]={"FASTCOMPONENT","SLOWCOMPONENT"};
  for(int iTag=0; iTag<2; iTag++)
  {
    const G4MaterialPropertyVector &component = *scintillatorPropertiesTable->GetProperty(tag[iTag].c_str());

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
  G4MaterialPropertyVector *rindexScintillator = scintillatorPropertiesTable->GetProperty("RINDEX");
  _cerenkovEnergyMinScintillator = rindexScintillator->GetMinLowEdgeEnergy();
  _cerenkovEnergyMaxScintillator = rindexScintillator->GetMaxLowEdgeEnergy();
  _maxRIndexScintillator = rindexScintillator->GetMaxValue();

  G4Material *fiber = G4Material::GetMaterial("PMMA",true); //fiber
  G4MaterialPropertiesTable* fiberPropertiesTable = fiber->GetMaterialPropertiesTable();
  _rindexFiber = fiberPropertiesTable->GetProperty("RINDEX");
  _cerenkovEnergyMinFiber = _rindexFiber->GetMinLowEdgeEnergy();
  _cerenkovEnergyMaxFiber = _rindexFiber->GetMaxLowEdgeEnergy();
  _maxRIndexFiber = _rindexFiber->GetMaxValue();
}

bool WLSPrimaryGeneratorAction::SetNextBins()
{
  const std::vector<double> &xBins     = WLSDetectorConstruction::Instance()->GetXBins();
  const std::vector<double> &yBins     = WLSDetectorConstruction::Instance()->GetYBins();
  const std::vector<double> &zBins     = WLSDetectorConstruction::Instance()->GetZBins();
  const std::vector<double> &betaBins  = WLSDetectorConstruction::Instance()->GetBetaBins();
  const std::vector<double> &thetaBins = WLSDetectorConstruction::Instance()->GetThetaBins();
  const std::vector<double> &phiBins   = WLSDetectorConstruction::Instance()->GetPhiBins();
  const std::vector<double> &rBins     = WLSDetectorConstruction::Instance()->GetRBins();

  int nXBins=xBins.size()-1;   //e.g. 3 bins need 4 entries in the vector (for 3 bin boundaries)
  int nYBins=yBins.size()-1;
  int nZBins=zBins.size()-1;
  int nBetaBins=betaBins.size()-1;
  int nThetaBins=thetaBins.size()-1;
  int nPhiBins=phiBins.size()-1;
  int nRBins=rBins.size()-1;

//bin# = 40*120*xBin + 120*yBin + zBin
  _currentBin++;
  if(_verbose) std::cout<<"Sim Type: "<<_simType<<"    Current Bin: "<<_currentBin<<"     Number Of Photons: "<<_numberOfPhotons<<std::endl;

  switch(_simType)
  {
    case 0: //scintillation in scintillator
    case 1: //Cerenkov in scintillator
            if(_currentBin<0 || _currentBin>=nZBins*nYBins*nXBins) return false;

            {
              int xbin = (_currentBin / (nZBins*nYBins)) % nXBins;
              int ybin = (_currentBin / nZBins) % nYBins;
              int zbin = _currentBin % nZBins;
              _minBinX=xBins[xbin]/mm;
              _minBinY=yBins[ybin]/mm;
              _minBinZ=zBins[zbin]/mm;
              _maxBinX=xBins[xbin+1]/mm;
              _maxBinY=yBins[ybin+1]/mm;
              _maxBinZ=zBins[zbin+1]/mm;
              if(_verbose)
              {
                std::cout<<"X: "<<_minBinX<<" ... "<<_maxBinX<<"    bin# "<<xbin<<"/"<<nXBins<<std::endl;
                std::cout<<"Y: "<<_minBinY<<" ... "<<_maxBinY<<"    bin# "<<ybin<<"/"<<nYBins<<std::endl;
                std::cout<<"Z: "<<_minBinZ<<" ... "<<_maxBinZ<<"    bin# "<<zbin<<"/"<<nZBins<<std::endl;
              }
            }
            break;
    case 2: //Cerenkov in fiber 0
    case 3: //Cerenkov in fiber 1
            if(_currentBin<0 || _currentBin>=nZBins*nRBins*nPhiBins*nThetaBins*nBetaBins) return false;

            {
              int betabin = (_currentBin / (nZBins*nRBins*nPhiBins*nThetaBins)) % nBetaBins;
              int thetabin = (_currentBin / (nZBins*nRBins*nPhiBins)) % nThetaBins;
              int phibin = (_currentBin / (nZBins*nRBins)) % nPhiBins;
              int rbin = (_currentBin / nZBins) % nRBins;
              int zbin = _currentBin % nZBins;
              _minBinBeta=betaBins[betabin];
              _minBinTheta=thetaBins[thetabin];
              _minBinPhi=phiBins[phibin];
              _minBinR=rBins[rbin]/mm;
              _minBinZ=zBins[zbin]/mm;
              _maxBinBeta=betaBins[betabin+1];
              _maxBinTheta=thetaBins[thetabin+1];
              _maxBinPhi=phiBins[phibin+1];
              _maxBinR=rBins[rbin+1]/mm;
              _maxBinZ=zBins[zbin+1]/mm;
              if(_verbose)
              {
                std::cout<<"Beta: "<<_minBinBeta<<" ... "<<_maxBinBeta<<std::endl;
                std::cout<<"Theta: "<<_minBinTheta<<" ... "<<_maxBinTheta<<std::endl;
                std::cout<<"Phi: "<<_minBinPhi<<" ... "<<_maxBinPhi<<std::endl;
                std::cout<<"R: "<<_minBinR<<" ... "<<_maxBinR<<std::endl;
                std::cout<<"Z: "<<_minBinZ<<" ... "<<_maxBinZ<<std::endl;
              }
            }
            break;
  }

  WLSEventAction::Instance()->SetStartZ(((_maxBinZ+_minBinZ)/2.0)*mm);
  return(true);
}

//both scintillation photons and (in approximation) also Cerenkov photons
int WLSPrimaryGeneratorAction::GeneratePhotonsInScintillator(G4Event *anEvent, int generatedPhotons) 
{
  G4ThreeVector startPosition0(_minBinX,_minBinY,_minBinZ);
  G4Navigator* navigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  G4VPhysicalVolume *scintillator = WLSDetectorConstruction::Instance()->GetScintillatorVolume();

  double dx=_maxBinX-_minBinX;
  double dy=_maxBinY-_minBinY;
  double dz=_maxBinZ-_minBinZ;

  int actualNumberOfGeneratedPhotons=0;
  for(int i=0; i<generatedPhotons; i++) 
  {
//pick photon start position

    double deltaX=G4UniformRand()*dx;
    double deltaY=G4UniformRand()*dy;
    double deltaZ=G4UniformRand()*dz;
    G4ThreeVector deltaPos(deltaX,deltaY,deltaZ);
    G4ThreeVector startPosition = startPosition0+deltaPos;

    if(navigator->LocateGlobalPointAndSetup(startPosition,NULL,false,true)!=scintillator) continue;
    actualNumberOfGeneratedPhotons++;

//pick random photon energy

    double photonEnergy;
    if(_simType==0)  //scintillation
    {
      int iTag=1;  //slow
      if(G4UniformRand()<=_yieldRatio) iTag=0; //fast

      double integratedProbMax = _emissionIntegral[iTag].GetMaxValue();
      double integratedProb = G4UniformRand()*integratedProbMax;
      photonEnergy = _emissionIntegral[iTag].GetEnergy(integratedProb);
    }
    else  //cerenkov
    {
      photonEnergy = _cerenkovEnergyMinScintillator + G4UniformRand()*(_cerenkovEnergyMaxScintillator-_cerenkovEnergyMinScintillator);
    }

//pick random direction (uniformly distributed)

    double cost = 1. - 2.*G4UniformRand();
    double sint = std::sqrt((1.-cost)*(1.+cost));

    double phi = twopi*G4UniformRand();
    double sinp = std::sin(phi);
    double cosp = std::cos(phi);

    double px = sint*cosp;
    double py = sint*sinp;
    double pz = cost;
    G4ParticleMomentum photonMomentum(px, py, pz);

//pick random polarization (uniformly distributed)

    double sx = cost*cosp;
    double sy = cost*sinp; 
    double sz = -sint;

    G4ThreeVector photonPolarization(sx, sy, sz);
    G4ThreeVector perp = photonMomentum.cross(photonPolarization);

    phi = twopi*G4UniformRand();
    sinp = std::sin(phi);
    cosp = std::cos(phi);
    photonPolarization = cosp * photonPolarization + sinp * perp;
    photonPolarization = photonPolarization.unit();

//generate photon
    _particleGun->SetParticleMomentumDirection(photonMomentum);
    _particleGun->SetParticlePolarization(photonPolarization);
    _particleGun->SetParticleEnergy(photonEnergy);
    _particleGun->SetParticleTime(0);
    _particleGun->SetParticlePosition(startPosition);
    _particleGun->GeneratePrimaryVertex(anEvent); 
  }
  return actualNumberOfGeneratedPhotons;
}

int WLSPrimaryGeneratorAction::GenerateCerenkovPhotonsInFiber(G4Event *anEvent, int generatedPhotons)
{
  double dEnergy = _cerenkovEnergyMaxFiber - _cerenkovEnergyMinFiber;
  double dBeta = _maxBinBeta - _minBinBeta;
  double dTheta = _maxBinTheta - _minBinTheta;
  double dPhi = _maxBinPhi - _minBinPhi;

  if(_maxBinBeta*_maxRIndexFiber<=1.0) return 0;  //Cerenkov photons can't be produced 

  double fiberSeparation = WLSDetectorConstruction::Instance()->GetFiberSeparation()/mm;
  G4ThreeVector startPosition(_minBinR,_simType==2?-fiberSeparation/2.0:fiberSeparation/2.0,_minBinZ);
  double dR = _maxBinR - _minBinR;
  double dZ = _maxBinZ - _minBinZ;

  for(int i=0; i<generatedPhotons; i++) 
  {
    double beta=_minBinBeta+dBeta*G4UniformRand();
    double maxCos = 1. / (beta*_maxRIndexFiber); 
    if(maxCos>=1.0) {i--; continue;}
    double maxSin2 = (1.0 - maxCos) * (1.0 + maxCos);

//pick random photon energy, which leads to random rIndex, which leads to random angle
//could be simplified, since all RIndeces are equal

    double photonEnergy,cosTheta,sin2Theta;
    do   
    {
      photonEnergy = _cerenkovEnergyMinFiber + dEnergy*G4UniformRand(); 
      double rIndex = _rindexFiber->Value(photonEnergy);
      cosTheta = 1. / (beta*rIndex);  
      if(cosTheta>=1.0) continue;
      sin2Theta = (1.0 - cosTheta)*(1.0 + cosTheta);
    } while(G4UniformRand()*maxSin2 > sin2Theta);  //TODO: Why this weighting with respect to sin2Theta???

//random direction of photon on cone surface defined by Theta 
//(in coordinate system with primary particle direction aligned with the z axis)

    double phi = twopi*G4UniformRand();
    double sinPhi = std::sin(phi);
    double cosPhi = std::cos(phi);

    double sinTheta = std::sqrt(sin2Theta); 
    double px = sinTheta*cosPhi;
    double py = sinTheta*sinPhi;
    double pz = cosTheta;

    G4ParticleMomentum photonMomentum(px, py, pz);

//rotate momentum direction back to global reference system 
    
    G4ThreeVector particleMomentum;
    double particleTheta=_minBinTheta + G4UniformRand()*dTheta;
    double particlePhi=_minBinPhi + G4UniformRand()*dPhi;
    particleMomentum.setRThetaPhi(1.0, particleTheta, particlePhi);
    photonMomentum.rotateUz(particleMomentum);

//polarization of new photon 

    double sx = cosTheta*cosPhi;
    double sy = cosTheta*sinPhi; 
    double sz = -sinTheta;

    G4ThreeVector photonPolarization(sx, sy, sz);

//rotate back to global reference system 

    photonPolarization.rotateUz(particleMomentum);

//start position

    double deltaX=G4UniformRand()*dR;  //TODO: the probability should be weighted with respect to r
    double deltaZ=G4UniformRand()*dZ;
    G4ThreeVector deltaPos(deltaX,0,deltaZ);
	
//generate photon
    _particleGun->SetParticleMomentumDirection(photonMomentum);
    _particleGun->SetParticlePolarization(photonPolarization);
    _particleGun->SetParticleEnergy(photonEnergy);
    _particleGun->SetParticleTime(0);
    _particleGun->SetParticlePosition(startPosition+deltaPos);
    _particleGun->GeneratePrimaryVertex(anEvent); 
  }

  return generatedPhotons;
}

void WLSPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
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

    if(!SetNextBins()) return;

    if(_simType<2)
    {
      int generatedPhotons=_numberOfPhotons;
      int actualNumberOfGeneratedPhotons = GeneratePhotonsInScintillator(anEvent, generatedPhotons);
      WLSEventAction::Instance()->SetGeneratedOptPhotons(actualNumberOfGeneratedPhotons);
    }
    else
    {
      int generatedPhotons=_numberOfPhotons;
      int actualNumberOfGeneratedPhotons = GenerateCerenkovPhotonsInFiber(anEvent, generatedPhotons);
      WLSEventAction::Instance()->SetGeneratedOptPhotons(actualNumberOfGeneratedPhotons);
    }
  }

  if(_mode==0 || _mode==1)
  {
//    double beamsize=2.0*mm;
//    double beamsize=1.0*mm;
    double beamsize=0.0*mm;
    double x0 = -10*mm;
    double y0 = G4RandGauss::shoot(_randomEngine,_posY*mm,beamsize);    //0 is at center //-13mm is at fiber 0
    double barlength = WLSDetectorConstruction::Instance()->GetBarLength()/mm;
    double z0 = G4RandGauss::shoot(_randomEngine,(-barlength/2.0+_posZ)*mm,beamsize);  //pos mm from left side of counter
    _particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));

    _particleGun->GeneratePrimaryVertex(anEvent);  //not used for WLS fiber test

// for WLS fiber test
/*
    for(int ii=0; ii<1000; ii++) 
    {
      double cost = 1. - 2.*G4UniformRand();
      double sint = std::sqrt((1.-cost)*(1.+cost));
      double phi = twopi*G4UniformRand();
      double sinp = std::sin(phi);
      double cosp = std::cos(phi);
      double sx = cost*cosp;
      double sy = cost*sinp; 
      double sz = -sint;
      G4ThreeVector photonPolarization(sx, sy, sz);
      G4ThreeVector photonMomentum = _particleGun->GetParticleMomentumDirection();
      G4ThreeVector perp = photonMomentum.cross(photonPolarization);
      phi = twopi*G4UniformRand();
      sinp = std::sin(phi);
      cosp = std::cos(phi);
      photonPolarization = cosp * photonPolarization + sinp * perp;
      photonPolarization = photonPolarization.unit();
      _particleGun->SetParticlePolarization(photonPolarization);
      _particleGun->GeneratePrimaryVertex(anEvent);
    }
*/
  }
}

