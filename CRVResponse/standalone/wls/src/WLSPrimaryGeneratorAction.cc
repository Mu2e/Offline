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

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include <stdexcept>

WLSPrimaryGeneratorAction::WLSPrimaryGeneratorAction(WLSSteppingAction::simulationMode mode, 
                                                     int numberOfPhotons, int simType, int startBin, bool verbose, double posY, double posZ) : 
                                                                                            _mode(mode), 
                                                                                            _numberOfPhotons(numberOfPhotons),
                                                                                            _simType(simType),
                                                                                            _currentBin(startBin), 
                                                                                            _verbose(verbose),
                                                                                            _first(true),
                                                                                            _posY(posY), _posZ(posZ)
{
  _randomEngine = CLHEP::HepRandom::getTheEngine();
  _particleGun  = new G4ParticleGun(1);
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  if(_mode==WLSSteppingAction::CreateLookupTables)
  {
    G4ParticleDefinition* particle = particleTable->FindParticle("opticalphoton");
    _particleGun->SetParticleDefinition(particle);
  }
  if(_mode==WLSSteppingAction::UseGeantOnly || _mode==WLSSteppingAction::UseGeantAndLookupTables)
  {
#ifndef FIBERTEST
    G4ParticleDefinition* particle = particleTable->FindParticle("proton"); 
    _particleGun->SetParticleEnergy(120.*GeV);
#else
#pragma message "USING FIBERTEST"
    G4ParticleDefinition* particle = particleTable->FindParticle("opticalphoton");
    _particleGun->SetParticleEnergy(2.88*eV); 
#endif
    _particleGun->SetParticleDefinition(particle);
    _particleGun->SetParticleMomentumDirection(G4ThreeVector(1., 0., 0.));
    _particleGun->SetParticleTime(0);
  }  
}

WLSPrimaryGeneratorAction::~WLSPrimaryGeneratorAction()
{
  delete _particleGun;
}

void WLSPrimaryGeneratorAction::BuildEmissionSpectrum()
{
  G4Material* scintillator = G4Material::GetMaterial("PolystyreneScint",true);
  G4MaterialPropertiesTable* scintillatorPropertiesTable = scintillator->GetMaterialPropertiesTable();

  _scintillationRiseTime = scintillatorPropertiesTable->GetConstProperty("FASTSCINTILLATIONRISETIME");
  _scintillationDecayTime = scintillatorPropertiesTable->GetConstProperty("FASTTIMECONSTANT");

//scintillation: build the emission integral
//from G4Scintillation::BuildThePhysicsTable()
  const G4MaterialPropertyVector &component = *scintillatorPropertiesTable->GetProperty("FASTCOMPONENT");  //ignore slow component

  double currentEnergy = component.Energy(0);
  double currentDifferentialProb = component[0];
  double currentIntegratedProb = 0.0;
  _emissionIntegral.InsertValues(currentEnergy, currentIntegratedProb);

  for(size_t i=1; i<component.GetVectorLength(); i++)
  {
    double prevEnergy  = currentEnergy;
    double prevDifferentialProb = currentDifferentialProb;
    double prevIntegratedProb = currentIntegratedProb;

    currentEnergy = component.Energy(i);
    currentDifferentialProb = component[i];
    currentIntegratedProb = 0.5 * (prevDifferentialProb + currentDifferentialProb);
    currentIntegratedProb = prevIntegratedProb + (currentEnergy - prevEnergy) * currentIntegratedProb;
    _emissionIntegral.InsertValues(currentEnergy, currentIntegratedProb);
  }

//Cerenkov: get the properties to generate the Cerenkov photons in scintillator
  _rindexScintillator = scintillatorPropertiesTable->GetProperty("RINDEX");
  _cerenkovEnergyMinScintillator = _rindexScintillator->GetMinLowEdgeEnergy();
  _cerenkovEnergyMaxScintillator = _rindexScintillator->GetMaxLowEdgeEnergy();
  _maxRIndexScintillator = (*_rindexScintillator)[0];
  for(size_t i=1; i<_rindexScintillator->GetVectorLength(); i++)
  {
    double currentRindex=(*_rindexScintillator)[i];
    if(currentRindex>_maxRIndexScintillator) _maxRIndexScintillator=currentRindex;
  }

//Cerenkov: get the properties to generate the Cerenkov photons in fiber (assume constant rindex)
  G4Material *fiber = G4Material::GetMaterial("PolystyreneFiber",true); //fiber
  G4MaterialPropertiesTable* fiberPropertiesTable = fiber->GetMaterialPropertiesTable();
  _rindexFiber = fiberPropertiesTable->GetProperty("RINDEX");
  _cerenkovEnergyMinFiber = _rindexFiber->GetMinLowEdgeEnergy();
  _cerenkovEnergyMaxFiber = _rindexFiber->GetMaxLowEdgeEnergy();
  _maxRIndexFiber = (*_rindexFiber)[0];
  for(size_t i=1; i<_rindexFiber->GetVectorLength(); i++)
  {
    double currentRindex=(*_rindexFiber)[i];
    if(currentRindex>_maxRIndexFiber) _maxRIndexFiber=currentRindex;
  }
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
  if(_verbose) std::cout<<"Sim Type: "<<_simType<<"    Current Bin: "<<_currentBin<<"     Number Of Photons: "<<_numberOfPhotons<<std::endl;
  WLSEventAction::Instance()->SetBinNumber(_currentBin);

  switch(_simType)
  {
    case 0: //scintillation in scintillator
            if(_currentBin<0 || _currentBin>=nZBins*nYBins*nXBins) return false;
            {
              int xbin = (_currentBin / (nZBins*nYBins)) % nXBins;
              int ybin = (_currentBin / nZBins) % nYBins;
              int zbin = _currentBin % nZBins;
              _minBinX=xBins[xbin]/mm;
              _minBinY=yBins[ybin]/mm;
              _minBinZ=zBins[zbin]/mm;
              _maxBinX=xBins[xbin+1]/mm;  //this is not a problem, since there is always one more entry than e.g. nXBins
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
    case 1: //Cerenkov in scintillator
            if(_currentBin<0 || _currentBin>=nZBins*nYBins*nXBins*nBetaBins) return false;
            {
              int xbin = (_currentBin / (nBetaBins*nZBins*nYBins)) % nXBins;
              int ybin = (_currentBin / (nBetaBins*nZBins)) % nYBins;
              int zbin = (_currentBin / nBetaBins) % nZBins;
              int betabin = _currentBin % nBetaBins;
              _minBinX=xBins[xbin]/mm;
              _minBinY=yBins[ybin]/mm;
              _minBinZ=zBins[zbin]/mm;
              _maxBinX=xBins[xbin+1]/mm;  //this is not a problem, since there is always one more entry than e.g. nXBins
              _maxBinY=yBins[ybin+1]/mm;
              _maxBinZ=zBins[zbin+1]/mm;
              _minBinBeta=betaBins[betabin];
              _maxBinBeta=betaBins[betabin+1];
              if(_verbose)
              {
                std::cout<<"X: "<<_minBinX<<" ... "<<_maxBinX<<"    bin# "<<xbin<<"/"<<nXBins<<std::endl;
                std::cout<<"Y: "<<_minBinY<<" ... "<<_maxBinY<<"    bin# "<<ybin<<"/"<<nYBins<<std::endl;
                std::cout<<"Z: "<<_minBinZ<<" ... "<<_maxBinZ<<"    bin# "<<zbin<<"/"<<nZBins<<std::endl;
                std::cout<<"Beta: "<<_minBinBeta<<" ... "<<_maxBinBeta<<std::endl;
              }
            }
            break;
    case 2: //Cerenkov in fiber (only the fiber at +y is simulated)
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

  _currentBin++;
  return(true);
}

//scintillation photons and Cerenkov photons in scintillator
int WLSPrimaryGeneratorAction::GeneratePhotonsInScintillator(G4Event *anEvent, int generatedPhotons) 
{
  G4ThreeVector startPosition0(_minBinX,_minBinY,_minBinZ);
  G4Navigator* navigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  G4VPhysicalVolume *scintillator = WLSDetectorConstruction::Instance()->GetScintillatorVolume();

  double dx=_maxBinX-_minBinX;
  double dy=_maxBinY-_minBinY;
  double dz=_maxBinZ-_minBinZ;
  double dBetaInverse=1.0/_minBinBeta-1.0/_maxBinBeta; //only relevant for Cerenkov  
  double dEnergy = _cerenkovEnergyMaxScintillator-_cerenkovEnergyMinScintillator;  //only relevant for Cerenkov

  if(_simType==1 && _maxBinBeta*_maxRIndexFiber<=1.0) return 0;  //Cerenkov photons can't be produced 

  int actualNumberOfGeneratedPhotons=0;
  for(int i=0; i<generatedPhotons; i++) 
  {
    double startTime=0;

//pick photon start position
    double deltaX=G4UniformRand()*dx;
    double deltaY=G4UniformRand()*dy;
    double deltaZ=G4UniformRand()*dz;
    G4ThreeVector deltaPos(deltaX,deltaY,deltaZ);
    G4ThreeVector startPosition = startPosition0+deltaPos;

    //check whether the random position is really inside of the scintillator
    if(navigator->LocateGlobalPointAndSetup(startPosition,NULL,false,true)!=scintillator) continue;
    actualNumberOfGeneratedPhotons++;

//pick random photon energy

    double photonEnergy;
    if(_simType==0)  //scintillation
    {
      //from G4Scintillation::PostStepDoIt()
      double integratedProbMax = _emissionIntegral.GetMaxValue();
      double integratedProb = G4UniformRand()*integratedProbMax;
      photonEnergy = _emissionIntegral.GetEnergy(integratedProb);

      //from G4Scintillation::sample_time()
      while(1)
      {
        double ran1 = G4UniformRand();
        double ran2 = G4UniformRand();
        startTime = -1.0*_scintillationDecayTime*std::log(1-ran1);
        double tempvar = 1.0-std::exp(-1.0*startTime/_scintillationRiseTime);
        if(ran2 <= tempvar) break;
      }
    }
    else  //cerenkov
    {
      double BetaInverse=1.0/_maxBinBeta + G4UniformRand()*dBetaInverse;   
      if(BetaInverse>=_maxRIndexScintillator) {i--; continue;}

      //from G4Cerenkov::PostStepDoIt()
      double maxCos = BetaInverse / _maxRIndexScintillator; 
      double maxSin2 = (1.0 - maxCos) * (1.0 + maxCos);
      double sampledRI, cosTheta, sin2Theta;
      do
      {
          photonEnergy = _cerenkovEnergyMinScintillator + G4UniformRand()*dEnergy; 
          sampledRI = _rindexScintillator->Value(photonEnergy);
          cosTheta = BetaInverse / sampledRI;  
          sin2Theta = (1.0 - cosTheta)*(1.0 + cosTheta);
      } while (G4UniformRand()*maxSin2 > sin2Theta);
    }

//As a simplification, a uniformly distributed direction and polarization is used from the scintillation and Cerenkov photons.
//While this is correct for the scintillation, it is wrong for the Cerenkov photons. However, the correct calculation
//would require the direction and speed of the charged particle, which would require more dimensions in the lookup table.
//The number of Cerenkov photons in the scintillator is small compared to the number of scintillator photons, 
//and most Cerenkov photons are absorbed and wavelength shifted right away with a uniformly distributed direction and polarization. 

//pick random direction (uniformly distributed) from G4Scintillation::PostStepDoIt()

    double cost = 1. - 2.*G4UniformRand();
    double sint = std::sqrt((1.-cost)*(1.+cost));

    double phi = CLHEP::twopi*G4UniformRand();
    double sinp = std::sin(phi);
    double cosp = std::cos(phi);

    double px = sint*cosp;
    double py = sint*sinp;
    double pz = cost;
    G4ParticleMomentum photonMomentum(px, py, pz);

//pick random polarization (uniformly distributed) from G4Scintillation::PostStepDoIt()

    double sx = cost*cosp;
    double sy = cost*sinp; 
    double sz = -sint;

    G4ThreeVector photonPolarization(sx, sy, sz);
    G4ThreeVector perp = photonMomentum.cross(photonPolarization);

    phi = CLHEP::twopi*G4UniformRand();
    sinp = std::sin(phi);
    cosp = std::cos(phi);
    photonPolarization = cosp * photonPolarization + sinp * perp;
    photonPolarization = photonPolarization.unit();

//generate photon
    _particleGun->SetParticleMomentumDirection(photonMomentum);
    _particleGun->SetParticlePolarization(photonPolarization);
    _particleGun->SetParticleEnergy(photonEnergy);
    _particleGun->SetParticleTime(startTime);
    _particleGun->SetParticlePosition(startPosition);
    _particleGun->GeneratePrimaryVertex(anEvent); 
  }
  return actualNumberOfGeneratedPhotons;
}

//Cerenkov photons in fiber
int WLSPrimaryGeneratorAction::GenerateCerenkovPhotonsInFiber(G4Event *anEvent, int generatedPhotons)
{
  double dEnergy = _cerenkovEnergyMaxFiber - _cerenkovEnergyMinFiber;
  double dBetaInverse=1.0/_minBinBeta-1.0/_maxBinBeta;
  double dTheta = _maxBinTheta - _minBinTheta;
  double dPhi = _maxBinPhi - _minBinPhi;

  if(_maxBinBeta*_maxRIndexFiber<=1.0) return 0;  //Cerenkov photons can't be produced 

  double fiberSeparation = WLSDetectorConstruction::Instance()->GetFiberSeparation()/mm;
  G4ThreeVector startPosition(_minBinR,fiberSeparation/2.0,_minBinZ);  //take only the fiber at +y

  for(int i=0; i<generatedPhotons; i++) 
  {
    //from G4Cerenkov::PostStepDoIt() 
    //as a simplification, the different index of refraction in the tiny sections of the cladding around the core is ignored.
    double BetaInverse=1.0/_maxBinBeta + G4UniformRand()*dBetaInverse;   
    if(BetaInverse>=_maxRIndexFiber) {i--; continue;}

    double maxCos = BetaInverse / _maxRIndexFiber; 
    double maxSin2 = (1.0 - maxCos) * (1.0 + maxCos);
    double sampledRI, cosTheta, sin2Theta;
    double photonEnergy=0;
    do   
    {
      photonEnergy = _cerenkovEnergyMinFiber + dEnergy*G4UniformRand(); 
      sampledRI = _rindexFiber->Value(photonEnergy);
      cosTheta = BetaInverse / sampledRI;  
      sin2Theta = (1.0 - cosTheta)*(1.0 + cosTheta);
    } while(G4UniformRand()*maxSin2 > sin2Theta);

//random direction of photon on cone surface defined by Theta 
//(in coordinate system with primary particle direction aligned with the z axis)

    double phi = CLHEP::twopi*G4UniformRand();
    double sinPhi = std::sin(phi);
    double cosPhi = std::cos(phi);

    double sinTheta = std::sqrt(sin2Theta); 
    double px = sinTheta*cosPhi;
    double py = sinTheta*sinPhi;
    double pz = cosTheta;

    G4ParticleMomentum photonMomentum(px, py, pz);

//randomly generate a direction (in the global reference system) of a primary particle 
//which causes the Cerenkov photon
//(theta=0 is in +z direction, phi=0 is in +x direction, phi=pi/2 is in +y direction)
    
    G4ThreeVector particleMomentum;
    double particleTheta=_minBinTheta + G4UniformRand()*dTheta;
    double particlePhi=_minBinPhi + G4UniformRand()*dPhi;
    particleMomentum.setRThetaPhi(1.0, particleTheta, particlePhi);

//rotate photon momentum direction back to global reference system 

    photonMomentum.rotateUz(particleMomentum);

//polarization of photon 

    double sx = cosTheta*cosPhi;
    double sy = cosTheta*sinPhi; 
    double sz = -sinTheta;

    G4ThreeVector photonPolarization(sx, sy, sz);

//rotate photon polarization back to global reference system 

    photonPolarization.rotateUz(particleMomentum);

//start position 
//pick a random point within the shell of x=_minBinR and x=_maxBinR
//at fixed y (at fiber center)
//deltaX is generated in such a way that one has a constant probability of hitting a specific point in the shell
//A = 2*pi*int(r_1,r_2,r*dr)
//F(r) = int(r_1,r,r*dr)/int(r_1,r_2,r*dr)
//F(r) = ( r^2 - r_1^2 ) / ( r_2^2 - r_1^2 )
//uniformly distributed random number u = F(r)
//u * ( r_2^2 - r_1^2 ) + r_1^2 =  r^2
    double deltaX = sqrt((_maxBinR*_maxBinR-_minBinR*_minBinR)*G4UniformRand() + _minBinR*_minBinR);
    double deltaZ=G4UniformRand()*(_maxBinZ - _minBinZ);
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
  if(_mode==WLSSteppingAction::CreateLookupTables)
  {
    if(_first)
    {
      _first=false;

      G4Material* material = G4Material::GetMaterial("PolystyreneFiber",true);
      G4MaterialPropertiesTable* materialPropertiesTable = material->GetMaterialPropertiesTable();
      materialPropertiesTable->AddConstProperty("WLSTIMECONSTANT", 0.0);  //set WLS decay time to 0
                                                                          //so that the true travel time can be determined

      BuildEmissionSpectrum();
    }

    if(!SetNextBins()) return;

    int generatedPhotons=0;
    if(_simType<2) generatedPhotons = GeneratePhotonsInScintillator(anEvent, _numberOfPhotons);
    else generatedPhotons = GenerateCerenkovPhotonsInFiber(anEvent, _numberOfPhotons);
    WLSEventAction::Instance()->SetGeneratedOptPhotons(generatedPhotons);
  }

  if(_mode==WLSSteppingAction::UseGeantOnly || _mode==WLSSteppingAction::UseGeantAndLookupTables)
  {
//    double beamsize=0.5*mm;  //resolution of wire chamber
    double beamsize=0.0*mm;
    double x0 = -10*mm;
    double y0 = G4RandGauss::shoot(_randomEngine,_posY*mm,beamsize);    //0 is at center //-13mm is at fiber 0
    double scintillatorHalfLength = WLSDetectorConstruction::Instance()->GetScintillatorHalfLength()/mm;
    double z0 = G4RandGauss::shoot(_randomEngine,(-scintillatorHalfLength+_posZ)*mm,beamsize);  //pos mm from left side of counter
    _particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));

#ifndef FIBERTEST
    _particleGun->GeneratePrimaryVertex(anEvent);  //not used for WLS fiber test
#else
    for(int ii=0; ii<4000; ii++) 
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
#endif
  }
}

