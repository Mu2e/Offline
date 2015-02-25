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

#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include <sstream>

#include <TH3D.h>
#include <TNtuple.h>

/*
#include "G4LossTableManager.hh"
#include "G4EmSaturation.hh"
*/

WLSSteppingAction* WLSSteppingAction::_fgInstance = NULL;

WLSSteppingAction::WLSSteppingAction(int mode, const std::string &lookupFileName) : _mode(mode), _engine(0), _randFlat(_engine)
{
  _fgInstance = this;

  if(_mode==0)
  {
    _crvPhotonArrivals = std::unique_ptr<MakeCrvPhotonArrivals>(new MakeCrvPhotonArrivals(_randFlat));
    _crvPhotonArrivals->LoadLookupTable(lookupFileName);
  }
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
     G4String thePostPVname = " ";
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
      _crvPhotonArrivals->SetScintillationYield(_scintillationYield);
      _crvPhotonArrivals->SetScintillatorDecayTimeFast(_scintillatorDecayTimeFast);
      _crvPhotonArrivals->SetScintillatorDecayTimeSlow(_scintillatorDecayTimeSlow);
      _crvPhotonArrivals->SetFiberDecayTime(_fiberDecayTime);
    }

    if(PDGcode!=0)  //ignore optical photons
    {
//std::cout<<"G4 : "<<G4LossTableManager::Instance()->EmSaturation()->VisibleEnergyDeposition(theStep)<<std::endl;
      _crvPhotonArrivals->MakePhotons(p1, p2, t1, t2,  
                            PDGcode, beta, charge,
                            energyDepositedTotal,
                            energyDepositedNonIonizing);
 
      for(int SiPM=0; SiPM<4; SiPM++)
      {
        std::vector<double> times=_crvPhotonArrivals->GetArrivalTimes(SiPM);
        _arrivalTimes[1][SiPM].assign(times.begin(),times.end());
      }
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
    if(i==0) _fiberEmissions[SiPM].clear();
  }

  _wlsTracks.clear();

  if(_mode==0)
  {
    _crvPhotonArrivals->Reset();
  }
}

