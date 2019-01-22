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

WLSSteppingAction::WLSSteppingAction(simulationMode mode, const std::string &lookupFileName, const std::string &visibleEnergyAdjustmentFileName) : 
                                                         _mode(mode), _engine(0), _randFlat(_engine), _randGaussQ(_engine), _randPoissonQ(_engine)
                                                                                 //lookupFileName and visibleEnergyAdjustmentFileName
                                                                                 //only used for simulationMode::UseGeantAndLookupTables
{
  _fgInstance = this;

  //load lookup tables
  if(_mode==UseGeantAndLookupTables)
  {
    _crvPhotons = std::unique_ptr<mu2eCrv::MakeCrvPhotons>(new mu2eCrv::MakeCrvPhotons(_randFlat, _randGaussQ, _randPoissonQ));
    _crvPhotons->LoadLookupTable(lookupFileName);
    _crvPhotons->LoadVisibleEnergyAdjustmentTable(visibleEnergyAdjustmentFileName);
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
    }
  }


  G4VPhysicalVolume* thePostPV = theStep->GetPostStepPoint()->GetPhysicalVolume();
  if(thePostPV)
  {
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
        std::cout<<theStep->GetTrack()->GetCreatorProcess()->GetProcessName();
      }
      std::cout<<std::endl;
    }
*/

    if(theStatus==Detection)
    {
         if(thePostPV->GetName()=="PhotonDet")
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
           _fiberEmissions[thePostPV->GetCopyNo()].push_back(numberOfFiberEmissions);
           _arrivalTimes[thePostPV->GetCopyNo()].push_back(theStep->GetPostStepPoint()->GetGlobalTime());
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
    double energyDepositedTotal= theStep->GetTotalEnergyDeposit();
    double energyDepositedNonIonizing = theStep->GetNonIonizingEnergyDeposit();
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
                            PDGcode, beta, charge,
                            energyDepositedTotal,
                            energyDepositedNonIonizing,
                            trueStepLength,0,reflector);
 
      for(int SiPM=0; SiPM<4; SiPM++)
      {
        const std::vector<double> &times=_crvPhotons->GetArrivalTimes(SiPM);
        _arrivalTimesFromLookupTables[SiPM].insert(_arrivalTimesFromLookupTables[SiPM].end(),times.begin(),times.end());
      }
    }
  }

//  ShowVisibleEnergyTable(theStep);

}

const std::vector<double> &WLSSteppingAction::GetArrivalTimes(int SiPM)
{
  return _arrivalTimes[SiPM];
}

const std::vector<double> &WLSSteppingAction::GetArrivalTimesFromLookupTables(int SiPM)
{
  return _arrivalTimesFromLookupTables[SiPM];
}

const std::vector<int> &WLSSteppingAction::GetFiberEmissions(int SiPM)
{
  return _fiberEmissions[SiPM];
}

void WLSSteppingAction::Reset()
{
  for(int SiPM=0; SiPM<4; SiPM++)
  {
    _arrivalTimes[SiPM].clear();
    _arrivalTimes[SiPM].reserve(10000);
    _arrivalTimesFromLookupTables[SiPM].clear();
    _arrivalTimesFromLookupTables[SiPM].reserve(10000);
    _fiberEmissions[SiPM].clear();
    _fiberEmissions[SiPM].reserve(10000);
  }

  _wlsTrackParents.clear();
}

void WLSSteppingAction::ShowVisibleEnergyTable(const G4Step *theStep)
{
  if(theStep->GetTotalEnergyDeposit()==0) return;

  G4Material* material = const_cast<G4Material*>(theStep->GetTrack()->GetMaterialCutsCouple()->GetMaterial());
  double BirksConstant = material->GetIonisation()->GetBirksConstant();
  std::cout<<material->GetName()<<"  Birks Constant: "<<BirksConstant<<std::endl;

  std::cout<<"PDGcode: "<<theStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding()<<std::endl;
  std::cout<<"Original Energy Deposition (G4): "<<theStep->GetTotalEnergyDeposit()<<std::endl;
  std::cout<<"Original Nonionizting Energy Deposition (G4): "<<theStep->GetNonIonizingEnergyDeposit()<<std::endl;
  std::cout<<"Visible Energy Deposition (G4): "<<G4LossTableManager::Instance()->EmSaturation()->VisibleEnergyDepositionAtAStep(theStep)<<std::endl;
  std::cout<<"Step Length: "<<theStep->GetStepLength()<<std::endl;
  const G4ThreeVector &p1 = theStep->GetPreStepPoint()->GetPosition();
  const G4ThreeVector &p2 = theStep->GetPostStepPoint()->GetPosition();
  std::cout<<"             "<<(p1-p2).mag()<<std::endl;

  std::cout<<"ELECTRON RANGE"<<std::endl;
  for(double e=0.001*eV; e<1.0*TeV; e*=1.2)
  {
    std::cout<<material->GetName();
    std::cout<<"  Energy: "<<e;
    std::cout<<"  Range: "<<G4LossTableManager::Instance()->GetRange(G4Electron::Electron(), e, theStep->GetTrack()->GetMaterialCutsCouple());
    std::cout<<"  Energy/Range: "<<e/G4LossTableManager::Instance()->GetRange(G4Electron::Electron(), e, theStep->GetTrack()->GetMaterialCutsCouple());
    std::cout<<std::endl;
  }

  std::cout<<"PROTON RANGE"<<std::endl;
  double ratio = 0;
  double chargeSq = 0; 
  double norm = 0.0;
  const G4ElementVector* theElementVector = material->GetElementVector();
  const double* theAtomNumDensityVector = material->GetVecNbOfAtomsPerVolume();
  size_t nelm = material->GetNumberOfElements();
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
  for(double e=0.001*eV; e<1.0*TeV; e*=1.2)
  {
    std::cout<<material->GetName();
    std::cout<<"  Energy: "<<e;
    std::cout<<"  Range: "<<G4LossTableManager::Instance()->GetRange(G4Proton::Proton(), e*ratio, theStep->GetTrack()->GetMaterialCutsCouple());
    std::cout<<"  Energy/(Range/chargeSq): "<<e/(G4LossTableManager::Instance()->GetRange(G4Proton::Proton(), e*ratio, theStep->GetTrack()->GetMaterialCutsCouple())/chargeSq);
    std::cout<<std::endl;
  }
}

