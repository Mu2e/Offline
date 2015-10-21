#include "G4Run.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"

#include "WLSEventAction.hh"
#include "WLSSteppingAction.hh"

#include "MakeCrvPhotonArrivals.hh"

#include "G4ProcessManager.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4RunManager.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"

#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include <sstream>

#include <TH3D.h>
#include <TNtuple.h>

#include "G4LossTableManager.hh"
#include "G4NistManager.hh"
#include "G4EmSaturation.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"

WLSSteppingAction* WLSSteppingAction::_fgInstance = NULL;

WLSSteppingAction::WLSSteppingAction(int mode, const std::string &lookupFileName) : _mode(mode), _engine(0), _randFlat(_engine)
{
  _fgInstance = this;

  if(_mode==0)
  {
    _crvPhotonArrivals = std::unique_ptr<MakeCrvPhotonArrivals>(new MakeCrvPhotonArrivals(_randFlat));
    _crvPhotonArrivals->LoadLookupTable(lookupFileName);
  }

  Reset();
}

WLSSteppingAction::~WLSSteppingAction()
{
}

void WLSSteppingAction::UserSteppingAction(const G4Step* theStep)
{
  G4OpBoundaryProcessStatus theStatus = Undefined;

  G4ProcessManager* OpManager = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();

  if (OpManager) 
  {
     G4int MAXofPostStepLoops = OpManager->GetPostStepProcessVector()->entries();
     G4ProcessVector* fPostStepDoItVector = OpManager->GetPostStepProcessVector(typeDoIt);

     for ( G4int i=0; i<MAXofPostStepLoops; i++) 
     {
         G4VProcess* fCurrentProcess = (*fPostStepDoItVector)[i];
         G4OpBoundaryProcess *opProcess = dynamic_cast<G4OpBoundaryProcess*>(fCurrentProcess);
         if (opProcess) { theStatus = opProcess->GetStatus(); break;}
     }
  }

  if(theStatus==Detection)
  {
     G4VPhysicalVolume* thePostPV = theStep->GetPostStepPoint()->GetPhysicalVolume();
     if(thePostPV)
     {
         if(thePostPV->GetName()=="PhotonDet")
         {
           _arrivalTimes[0][thePostPV->GetCopyNo()].push_back(theStep->GetPostStepPoint()->GetGlobalTime());
           if(_mode==-1)
           {
             int numberOfFiberEmissions=0;
             int trackID = theStep->GetTrack()->GetTrackID();
             while(1)
             {
                std::map<int,int>::const_iterator wlsIter=_wlsTracks.find(trackID);
                if(wlsIter!=_wlsTracks.end())
                {
                  numberOfFiberEmissions++;
                  trackID=wlsIter->second;  //parent ID
                }
                else break;
             }
             _fiberEmissions[thePostPV->GetCopyNo()].push_back(numberOfFiberEmissions);
           }
         }
     }
  }

  if(_mode==-1)
  {
     if(theStep->GetTrack()->GetCreatorProcess()!=NULL)
     {
       if(theStep->GetTrack()->GetCreatorProcess()->GetProcessName()=="OpWLS")
       {
         int trackID=theStep->GetTrack()->GetTrackID();
         int parentID=theStep->GetTrack()->GetParentID();
         _wlsTracks[trackID]=parentID;
       }
     }
  }

  if(_mode==0)
  {
    const G4ThreeVector &p1 = theStep->GetPreStepPoint()->GetPosition();
    const G4ThreeVector &p2 = theStep->GetPostStepPoint()->GetPosition();
    const double t1 = theStep->GetPreStepPoint()->GetGlobalTime();
    const double t2 = theStep->GetPostStepPoint()->GetGlobalTime();
    int PDGcode = theStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
    double beta = (theStep->GetPreStepPoint()->GetBeta() + theStep->GetPostStepPoint()->GetBeta())/2.0;
    double charge = theStep->GetTrack()->GetParticleDefinition()->GetPDGCharge();
    double energyDepositedTotal= theStep->GetTotalEnergyDeposit();
    double energyDepositedNonIonizing = theStep->GetNonIonizingEnergyDeposit();

    static bool first=true;
    if(first)
    {
      first=false;

      G4Material* scintillator = G4Material::GetMaterial("Polystyrene",true);
      G4MaterialPropertiesTable* scintillatorPropertiesTable = scintillator->GetMaterialPropertiesTable();
      double scintillationYield = scintillatorPropertiesTable->GetConstProperty("SCINTILLATIONYIELD");
      double scintillatorDecayTimeFast = scintillatorPropertiesTable->GetConstProperty("FASTTIMECONSTANT");
      double scintillatorDecayTimeSlow = scintillatorPropertiesTable->GetConstProperty("SLOWTIMECONSTANT");

      G4Material* fiber = G4Material::GetMaterial("PMMA",true);
      G4MaterialPropertiesTable* fiberPropertiesTable = fiber->GetMaterialPropertiesTable();
      double fiberDecayTime = fiberPropertiesTable->GetConstProperty("WLSTIMECONSTANT");

      _crvPhotonArrivals->SetScintillationYield(scintillationYield);
      _crvPhotonArrivals->SetScintillatorDecayTimeFast(scintillatorDecayTimeFast);
      _crvPhotonArrivals->SetScintillatorDecayTimeSlow(scintillatorDecayTimeSlow);
      _crvPhotonArrivals->SetFiberDecayTime(fiberDecayTime);
    }

    if(PDGcode!=0)  //ignore optical photons
    {
      _crvPhotonArrivals->MakePhotons(p1, p2, t1, t2,  
                            PDGcode, beta, charge,
                            energyDepositedTotal,
                            energyDepositedNonIonizing);
 
      for(int SiPM=0; SiPM<4; SiPM++)
      {
        std::vector<double> times=_crvPhotonArrivals->GetArrivalTimes(SiPM);
        _arrivalTimes[1][SiPM].insert(_arrivalTimes[1][SiPM].end(),times.begin(),times.end());
      }
//      Test(theStep, PDGcode);
    }

  }
}

const std::vector<double> &WLSSteppingAction::GetArrivalTimes(int i, int SiPM)
{
  return _arrivalTimes[i][SiPM];
}

const std::vector<int> &WLSSteppingAction::GetFiberEmissions(int SiPM)
{
  return _fiberEmissions[SiPM];
}

void WLSSteppingAction::Reset()
{
  for(int i=0; i<2; i++)
  for(int SiPM=0; SiPM<4; SiPM++)
  {
    _arrivalTimes[i][SiPM].clear();
    _arrivalTimes[i][SiPM].reserve(10000);
    if(i==0) 
    {
     _fiberEmissions[SiPM].clear();
     _fiberEmissions[SiPM].reserve(10000);
    }
  }

  _wlsTracks.clear();
}

void WLSSteppingAction::Test(const G4Step *theStep, int PDGcode)
{
  std::cout<<"Original/Visible Energy Deposition (G4): "<<theStep->GetTotalEnergyDeposit()<<"/"<<G4LossTableManager::Instance()->EmSaturation()->VisibleEnergyDeposition(theStep)<<"   PDGcode: "<<PDGcode<<std::endl;

  G4Material* Polystyrene = G4Material::GetMaterial("Polystyrene",true);
  double BirksConstant = Polystyrene->GetIonisation()->GetBirksConstant();

  std::cout<<"ELECTRON RANGE"<<std::endl;
  for(double e=0.001*eV; e<100.0*TeV; e*=1.5)
  {
    std::cout<<theStep->GetTrack()->GetMaterialCutsCouple()->GetMaterial()->GetName();
    std::cout<<"  Energy: "<<e;
    std::cout<<"  Range: "<<G4LossTableManager::Instance()->GetRange(G4Electron::Electron(), e, theStep->GetTrack()->GetMaterialCutsCouple());
    std::cout<<"  Error: "<<BirksConstant*e/G4LossTableManager::Instance()->GetRange(G4Electron::Electron(), e, theStep->GetTrack()->GetMaterialCutsCouple());
    std::cout<<"  Fit: "<<BirksConstant*(27.0*exp(-0.247*pow(fabs(log(e)+8.2),1.6))+0.177);
    std::cout<<std::endl;
  }

  std::cout<<"PROTON RANGE"<<std::endl;
  double ratio = 0;
  double chargeSq = 0; 
  double norm = 0.0;
  const G4ElementVector* theElementVector = Polystyrene->GetElementVector();
  const double* theAtomNumDensityVector = Polystyrene->GetVecNbOfAtomsPerVolume();
  size_t nelm = Polystyrene->GetNumberOfElements();
  for(size_t i=0; i<nelm; ++i) 
  {
    const G4Element* elm = (*theElementVector)[i];
    double Z = elm->GetZ();
    double w = Z*Z*theAtomNumDensityVector[i];
    ratio += w/G4NistManager::Instance()->GetAtomicMassAmu(G4int(Z));
    chargeSq = Z*Z*w;
    norm += w;
  }
  ratio *= CLHEP::proton_mass_c2/norm;
  chargeSq /= norm;
  for(double e=1.0*eV; e<1.0*TeV; e*=2.0)
  {
    std::cout<<theStep->GetTrack()->GetMaterialCutsCouple()->GetMaterial()->GetName();
    std::cout<<"  Energy: "<<e;
    std::cout<<"  Range: "<<G4LossTableManager::Instance()->GetRange(G4Proton::Proton(), e*ratio, theStep->GetTrack()->GetMaterialCutsCouple());
    std::cout<<"  Error: "<<BirksConstant*e/(G4LossTableManager::Instance()->GetRange(G4Proton::Proton(), e*ratio, theStep->GetTrack()->GetMaterialCutsCouple())/chargeSq);
    std::cout<<std::endl;
  }
}

