#include "G4Run.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"

#include "WLSDetectorConstruction.hh"
#include "WLSEventAction.hh"
#include "WLSSteppingAction.hh"

#include "MakeCrvPhotons.hh"

#include "G4ProcessManager.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4RunManager.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"

#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include <sstream>

//#define PHOTONTEST
#ifdef PHOTONTEST
#include <TNtuple.h>
#endif

#include "G4LossTableManager.hh"
#include "G4NistManager.hh"
#include "G4EmSaturation.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"

WLSSteppingAction* WLSSteppingAction::_fgInstance = NULL;

WLSSteppingAction::WLSSteppingAction(simulationMode mode, const std::string &lookupFileName) : 
                                     _mode(mode), _engine(0), _randFlat(_engine), _randGaussQ(_engine), _randPoissonQ(_engine)
                                     //lookupFileName only used for simulationMode::UseGeantAndLookupTables
{
  _fgInstance = this;

  //load lookup tables
  if(_mode==UseGeantAndLookupTables)
  {
    _crvPhotons = std::unique_ptr<mu2eCrv::MakeCrvPhotons>(new mu2eCrv::MakeCrvPhotons(_randFlat, _randGaussQ, _randPoissonQ));
    _crvPhotons->LoadLookupTable(lookupFileName);
  }

#ifdef PHOTONTEST
  _ntuple = new TNtuple("CRVPhotons","CRVPhotons","SiPM:Energy:Length:StartZ:x:y:time:angle");
#endif
  Reset();
}

WLSSteppingAction::~WLSSteppingAction()
{
#ifdef PHOTONTEST
  _ntuple->SaveAs("CRVPhotons.root");
  delete _ntuple;
#endif
}

void WLSSteppingAction::UserSteppingAction(const G4Step* theStep)
{
  G4OpBoundaryProcessStatus theStatus = G4OpBoundaryProcessStatus::Undefined;

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

  //create a list of fiber photons and their parents to find the number of re-emissions
  if(theStep->GetTrack()->GetCreatorProcess()!=NULL)
  {
    if(theStep->GetTrack()->GetCreatorProcess()->GetProcessName()=="OpWLSY11")
    {
      int trackID=theStep->GetTrack()->GetTrackID();
      int parentID=theStep->GetTrack()->GetParentID();
      _wlsTrackParents[trackID]=parentID;
      if(_wlsTrackParents.find(parentID)==_wlsTrackParents.end()) _tracksGettingAbsorbedInFiber.insert(parentID);
    }
  }


  G4VPhysicalVolume* thePostPV = theStep->GetPostStepPoint()->GetPhysicalVolume();
  if(thePostPV)
  {
    if(theStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding()==0)
    {
      std::string currentVolume=thePostPV->GetName();
      std::string creationVolume=theStep->GetTrack()->GetLogicalVolumeAtVertex()->GetName();
      int trackID=theStep->GetTrack()->GetTrackID();
      if((currentVolume=="WLSFiber" || currentVolume.compare(0,4,"Clad")==0) && creationVolume=="Scintillator") _tracksHittingFiber.insert(trackID);
    }
    
/*
    double trueStepLength=theStep->GetStepLength();
    double stepLength=(theStep->GetPreStepPoint()->GetPosition()-theStep->GetPostStepPoint()->GetPosition()).mag();
    double diff=(trueStepLength-stepLength)/stepLength;
    if(diff>0.01)
    {
      std::cout<<"-----------  "<<trueStepLength<<"   "<<stepLength<<"     "<<diff<<"   ";
      std::cout<<"             "<<G4LossTableManager::Instance()->EmSaturation()->VisibleEnergyDepositionAtAStep(theStep)<<"   ";
      std::cout<<theStep->GetTrack()->GetParticleDefinition()->GetParticleName()<<std::endl;
    }
*/

/*
//    if(theStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding()!=0)
    {
      std::cout<<theStep->GetTrack()->GetTrackID()<<"  "<<theStep->GetTrack()->GetParentID()<<"   ";
      std::cout<<theStep->GetTrack()->GetParticleDefinition()->GetParticleName()<<"   ";
      std::cout<<theStep->GetTrack()->GetTrackLength()<<"   "<<theStep->GetTrack()->GetTotalEnergy()<<"   ";
      std::cout<<theStep->GetPreStepPoint()->GetPosition()<<"   ";
      std::cout<<theStep->GetPostStepPoint()->GetPosition()<<"   ";
      std::cout<<"       ";
      if(theStep->GetPreStepPoint()->GetMaterial()!=NULL)
      {
        std::cout<<theStep->GetPreStepPoint()->GetMaterial()->GetName()<<"  ";
      }
      if(theStep->GetPostStepPoint()->GetMaterial()!=NULL)
      {
        std::cout<<theStep->GetPostStepPoint()->GetMaterial()->GetName()<<"  ";
      }
      if(theStep->GetTrack()->GetCreatorProcess()!=NULL)
      {
        std::cout<<theStep->GetTrack()->GetCreatorProcess()->GetProcessName()<<"  ";
      }
      if(theStep->GetTrack()->GetLogicalVolumeAtVertex()!=NULL)
      {
        std::cout<<theStep->GetTrack()->GetLogicalVolumeAtVertex()->GetName();
      }
      std::cout<<std::endl;
    }
*/

    if(theStatus==Detection)
    {
         if(thePostPV->GetName()=="SiPM")
         {
//std::cout<<"DETECTION  "<<thePostPV->GetCopyNo()<<std::endl;
           //a photon reached a SiPM
#ifdef PHOTONTEST
           _ntuple->Fill((float)(thePostPV->GetCopyNo()),
                                 theStep->GetTrack()->GetTotalEnergy(),
                                 theStep->GetTrack()->GetTrackLength(),
                                 theStep->GetTrack()->GetVertexPosition().z(),
                                 theStep->GetPostStepPoint()->GetPosition().x(),
                                 theStep->GetPostStepPoint()->GetPosition().y(), //WLS fiber test
                                 theStep->GetPostStepPoint()->GetGlobalTime(),
                                 theStep->GetPostStepPoint()->GetMomentum().unit().dot(CLHEP::Hep3Vector(0.0,0.0,1.0)));
#endif

           //run through the list of parent photons to find the number of re-emissions
           int numberOfFiberEmissions=0;
           int trackID = theStep->GetTrack()->GetTrackID();
           while(1)
           {
             std::map<int,int>::const_iterator wlsIter=_wlsTrackParents.find(trackID);
             if(wlsIter!=_wlsTrackParents.end())
             {
               numberOfFiberEmissions++;
               trackID=wlsIter->second;  //parent ID
             }
             else break;
           }
           _photonInfo[thePostPV->GetCopyNo()].emplace_back(theStep->GetPostStepPoint()->GetGlobalTime(),numberOfFiberEmissions);
           if(numberOfFiberEmissions==0) _zeroFiberEmissions++;
         }
    }
  }

  //if a lookup table is used
  if(_mode==UseGeantAndLookupTables)
  {
    const G4ThreeVector &p1 = theStep->GetPreStepPoint()->GetPosition();
    const G4ThreeVector &p2 = theStep->GetPostStepPoint()->GetPosition();
    const double t1 = theStep->GetPreStepPoint()->GetGlobalTime();
    const double t2 = theStep->GetPostStepPoint()->GetGlobalTime();
    int PDGcode = theStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
    double beta = (theStep->GetPreStepPoint()->GetBeta() + theStep->GetPostStepPoint()->GetBeta())/2.0;
    double charge = theStep->GetTrack()->GetParticleDefinition()->GetPDGCharge();
    double visibleEnergyDeposited = G4LossTableManager::Instance()->EmSaturation()->VisibleEnergyDepositionAtAStep(theStep);
    double trueStepLength = theStep->GetStepLength();  //may be longer than (p1-p2).mag() due to scattering

    static bool first=true;
    if(first)
    {
      first=false;

      //these constants are extracted from the G4Material
      //in a real run, they would be provided by the fcl file
      G4Material* scintillator = G4Material::GetMaterial("PolystyreneScint",true);
      G4MaterialPropertiesTable* scintillatorPropertiesTable = scintillator->GetMaterialPropertiesTable();
      double scintillationYield = scintillatorPropertiesTable->GetConstProperty("SCINTILLATIONYIELD");
      _crvPhotons->SetScintillationYield(scintillationYield);
    }

    if(PDGcode!=0)  //ignore optical photons
    {
     int reflector = WLSDetectorConstruction::Instance()->GetReflectorOption();
      _crvPhotons->MakePhotons(p1, p2, t1, t2,  
                            beta, charge,
                            visibleEnergyDeposited,
                            trueStepLength,reflector);
 
      for(int SiPM=0; SiPM<4; SiPM++)
      {
        const std::vector<double> &times=_crvPhotons->GetArrivalTimes(SiPM);
        _arrivalTimesFromLookupTables[SiPM].insert(_arrivalTimesFromLookupTables[SiPM].end(),times.begin(),times.end());
      }
    }
  }
}

const std::vector<WLSSteppingAction::PhotonInfo> &WLSSteppingAction::GetPhotonInfo(int SiPM)
{
  return _photonInfo[SiPM];
}

const std::vector<double> &WLSSteppingAction::GetArrivalTimesFromLookupTables(int SiPM)
{
  return _arrivalTimesFromLookupTables[SiPM];
}

void WLSSteppingAction::Reset()
{
  for(int SiPM=0; SiPM<4; SiPM++)
  {
    _photonInfo[SiPM].clear();
    _photonInfo[SiPM].reserve(10000);
    _arrivalTimesFromLookupTables[SiPM].clear();
    _arrivalTimesFromLookupTables[SiPM].reserve(10000);
  }

  _wlsTrackParents.clear();
  _tracksGettingAbsorbedInFiber.clear();
  _tracksHittingFiber.clear();
  _zeroFiberEmissions=0;
}

void WLSSteppingAction::PrintFiberStats()
{
  std::cout<<"Full GEANT4:    Tracks hitting fiber: "<<_tracksHittingFiber.size()<<"    Tracks getting absorbed in fiber: "<<_tracksGettingAbsorbedInFiber.size();
  std::cout<<"       Photons detected which have not been wavelength shifted in fiber: "<<_zeroFiberEmissions<<std::endl;
}

