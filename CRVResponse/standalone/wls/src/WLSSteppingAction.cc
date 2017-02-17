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

WLSSteppingAction::WLSSteppingAction(int mode, const std::string &lookupFileName, const std::string &visibleEnergyAdjustmentFileName) : 
                                                                       _mode(mode), _engine(0), _randFlat(_engine), _randGaussQ(_engine), _randPoissonQ(_engine)
{
  _fgInstance = this;

  //load lookup tables
  if(_mode==1)
  {
    _crvPhotonArrivals = std::unique_ptr<mu2eCrv::MakeCrvPhotonArrivals>(new mu2eCrv::MakeCrvPhotonArrivals(_randFlat, _randGaussQ, _randPoissonQ));
    _crvPhotonArrivals->LoadLookupTable(lookupFileName);
    _crvPhotonArrivals->LoadVisibleEnergyAdjustmentTable(visibleEnergyAdjustmentFileName);
  }

  _ntuple = new TNtuple("CRVPhotons","CRVPhotons","SiPM:Energy:Length:StartZ"); //WLS fiber test
  Reset();
}

WLSSteppingAction::~WLSSteppingAction()
{
  _ntuple->SaveAs("CRVPhotons.root");  //WLS fiber test
  delete _ntuple;
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
           _ntuple->Fill((float)(thePostPV->GetCopyNo()),
                                 theStep->GetTrack()->GetTotalEnergy(),
                                 theStep->GetTrack()->GetTrackLength(),
                                 theStep->GetTrack()->GetVertexPosition().z()); //WLS fiber test
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

  //creating a list of fiber photons and their parents to find the number of re-emissions, while lookup tables are created
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

  //if a lookup table is used
  if(_mode==1)
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

      //these constances are extracted from the G4Material
      //in a real run, they would be provided by the fcl file
      G4Material* scintillator = G4Material::GetMaterial("Polystyrene",true);
      G4MaterialPropertiesTable* scintillatorPropertiesTable = scintillator->GetMaterialPropertiesTable();
      double scintillationYield = scintillatorPropertiesTable->GetConstProperty("SCINTILLATIONYIELD");
      double scintillatorBirksConstant = scintillator->GetIonisation()->GetBirksConstant();
      double scintillatorRatioFastSlow = scintillatorPropertiesTable->GetConstProperty("YIELDRATIO");
      double scintillatorDecayTimeFast = scintillatorPropertiesTable->GetConstProperty("FASTTIMECONSTANT");
      double scintillatorDecayTimeSlow = scintillatorPropertiesTable->GetConstProperty("SLOWTIMECONSTANT");

      G4Material* fiber = G4Material::GetMaterial("PMMA",true);
      G4MaterialPropertiesTable* fiberPropertiesTable = fiber->GetMaterialPropertiesTable();
      double fiberDecayTime = fiberPropertiesTable->GetConstProperty("WLSTIMECONSTANT");

      _crvPhotonArrivals->SetScintillationYield(scintillationYield);
      _crvPhotonArrivals->SetScintillatorBirksConstant(scintillatorBirksConstant);
      _crvPhotonArrivals->SetScintillatorRatioFastSlow(scintillatorRatioFastSlow);
      _crvPhotonArrivals->SetScintillatorDecayTimeFast(scintillatorDecayTimeFast);
      _crvPhotonArrivals->SetScintillatorDecayTimeSlow(scintillatorDecayTimeSlow);
      _crvPhotonArrivals->SetFiberDecayTime(fiberDecayTime);
    }

    if(PDGcode!=0)  //ignore optical photons
    {
      _crvPhotonArrivals->MakePhotons(p1, p2, t1, t2,  
                            PDGcode, beta, charge,
                            energyDepositedTotal,
                            energyDepositedNonIonizing,
                            trueStepLength);
 
      for(int SiPM=0; SiPM<4; SiPM++)
      {
        std::vector<double> times=_crvPhotonArrivals->GetArrivalTimes(SiPM);
        _arrivalTimes[1][SiPM].insert(_arrivalTimes[1][SiPM].end(),times.begin(),times.end());
      }
    }
  }

//  ShowVisibleEnergyTable(theStep);

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

void WLSSteppingAction::ShowVisibleEnergyTable(const G4Step *theStep)
{
  if(theStep->GetTotalEnergyDeposit()==0) return;

  G4Material* material = const_cast<G4Material*>(theStep->GetTrack()->GetMaterialCutsCouple()->GetMaterial());
  double BirksConstant = material->GetIonisation()->GetBirksConstant();
  std::cout<<material->GetName()<<"  Birks Constant: "<<BirksConstant<<std::endl;

  std::cout<<"PDGcode: "<<theStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding()<<std::endl;
  std::cout<<"Original Energy Deposition (G4): "<<theStep->GetTotalEnergyDeposit()<<std::endl;
  std::cout<<"Original Nonionizting Energy Deposition (G4): "<<theStep->GetNonIonizingEnergyDeposit()<<std::endl;
  std::cout<<"Visible Energy Deposition (G4): "<<G4LossTableManager::Instance()->EmSaturation()->VisibleEnergyDeposition(theStep)<<std::endl;
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

