// Called at every G4 step.
//
// Original author Rob Kutschke
//

// C++ includes
#include <cstdio>
#include <cmath>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "G4Step.hh"
#include "G4Threading.hh"

// Mu2e includes
#include "Mu2eG4/inc/Mu2eG4SteppingAction.hh"
#include "Mu2eG4/inc/Mu2eG4UserHelpers.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "Mu2eG4/inc/getPhysicalVolumeOrThrow.hh"
#include "Mu2eG4/inc/SimParticleHelper.hh"
#include "Mu2eG4/inc/UserTrackInformation.hh"
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"
#include "Mu2eG4/inc/Mu2eG4ResourceLimits.hh"
#include "Mu2eG4/inc/Mu2eG4TrajectoryControl.hh"

using namespace std;

namespace mu2e {

  Mu2eG4SteppingAction::Mu2eG4SteppingAction(const fhicl::ParameterSet& pset,
                                             const std::vector<double>& timeVDtimes,
                                             IMu2eG4Cut& steppingCuts,
                                             IMu2eG4Cut& commonCuts,
                                             const Mu2eG4TrajectoryControl& trajectoryControl,
                                             const Mu2eG4ResourceLimits& lim) :
    pset_(pset),

    steppingCuts_(&steppingCuts),
    commonCuts_(&commonCuts),

    mu2elimits_(&lim),

    numTrackSteps_(),
    numKilledTracks_(),
    stepLimitKillerVerbose_(pset.get<bool>("debug.stepLimitKillerVerbose")),

    // Things related to time virtual detector
    tvd_time_(timeVDtimes),
    tvd_collection_(nullptr),
    tvd_warning_printed_(false),

    trajectoryControl_(&trajectoryControl),

    // Default values for parameters that are optional in the run time configuration.
    _debugEventList(pset.get<std::vector<int> >("debug.eventList", std::vector<int>())),
    _debugTrackList(pset.get<std::vector<int> >("debug.trackList", std::vector<int>())),

    _spHelper()
  {
    if( (!tvd_time_.empty()) && (G4Threading::G4GetThreadId() <= 0) ) {
        G4cout << "Time virtual detector is enabled. Particles are recorded at";
        for( unsigned int i=0; i<tvd_time_.size(); ++i ) G4cout << " " << tvd_time_[i];
            G4cout << " ns" << G4endl;
    }
  }//end ctor

  // A helper function to manage the printout.
  void Mu2eG4SteppingAction::printit( G4String const& s,
                                      G4int id,
                                      G4ThreeVector const& pos,
                                      G4ThreeVector const& mom,
                                      double localTime,
                                      double globalTime ){

      // It is easier to line up printout in columns with printf than with cout.
      // fixme use G4cout to be compatible with MT
      printf ( "%-8s %4d %15.4f %15.4f %15.4f %15.4f %15.4f %15.4f %15.4f %13.4f %13.4f\n",
              s.data(), id,
              pos.x(), pos.y(), pos.z(),
              mom.x(), mom.y(), mom.z(),
              mom.mag(),
              localTime, globalTime);
 }

  void Mu2eG4SteppingAction::beginRun(PhysicsProcessInfo* processInfo,
                                      CLHEP::Hep3Vector const& mu2eOrigin ){
      
    _processInfo    = processInfo;
    _mu2eOrigin     =  mu2eOrigin;
  }

  void Mu2eG4SteppingAction::finishConstruction() {
      
    // We have to wait until G4 geometry is constructed
    // to get phys volume pointers that are used in the
    // volume to cut value map.
      for(const auto& spec: trajectoryControl_->perVolumeMinDistance()) {
          auto vol = getPhysicalVolumeOrThrow(spec.first);
          mcTrajectoryVolumePtDistances_[vol] = spec.second;
      }
  }
    
  void Mu2eG4SteppingAction::BeginOfTrack() {
      numTrackSteps_ = 0;
      const auto oldSize = _trajectory.size();
      _trajectory.clear();
      _trajectory.reserve(oldSize + oldSize/8);
  }

  void Mu2eG4SteppingAction::EndOfTrack() {
  }

    
  void Mu2eG4SteppingAction::BeginOfEvent(StepPointMCCollection& outputHits,
                                          const SimParticleHelper& spHelper) {
    numKilledTracks_ = 0;
    tvd_collection_  = &outputHits;
    tvd_warning_printed_ = false;
    _spHelper    = &spHelper;
  }

    
  void Mu2eG4SteppingAction::UserSteppingAction(const G4Step* step){

      numTrackSteps_++;
      
      G4Track* track = step->GetTrack();
      
      // Pre and post stepping points.
      G4StepPoint const* prept  = step->GetPreStepPoint();
      G4StepPoint const* postpt = step->GetPostStepPoint();
      
      /// FIXME: Optimization: if the current track (not step!) has
      /// initial momentum below the "g4.mcTrajectoryMomentumCut"
      /// threshold, we can skip adding trajectory points completely.
      /// NB: there is also "g4.saveTrajectoryMomentumCut".
      /// Do we need to look at both of them?
      //
      // Determine whether we add the current step to MCTrajectory
      const double mcTrajCurrentCut = mcTrajectoryMinDistanceCut(prept->GetPhysicalVolume());
      // In some cases we know to accept the point even without computing the distance
      bool computeMCTrajDistance = (!_trajectory.empty()) && (mcTrajCurrentCut > 0.);
      if(!computeMCTrajDistance || (((prept->GetPosition() - _mu2eOrigin) -  _trajectory.back().pos()).mag() >= mcTrajCurrentCut)) {
          _trajectory.emplace_back ( prept->GetPosition() - _mu2eOrigin,
                                    prept->GetGlobalTime(),
                                    prept->GetKineticEnergy()
                                    );
      }

      // Save hits in time virtual detector
      for( unsigned int i=0; i<tvd_time_.size(); ++i ) {
          if( prept->GetGlobalTime()<=tvd_time_[i] && postpt->GetGlobalTime()>tvd_time_[i] ) {
              addTimeVDHit(step,i+1);
          }
      }
      
      if(steppingCuts_->steppingActionCut(step)) {
          killTrack(track, ProcessCode::mu2eKillerVolume, fStopAndKill);
      } else if(commonCuts_->steppingActionCut(step)) {
          killTrack(track, ProcessCode::mu2eKillerVolume, fStopAndKill);
      } else if(killTooManySteps(track)) {
          killTrack( track, ProcessCode::mu2eMaxSteps, fStopAndKill);
      }
      
      //----------------------------------------------------------------
      // Do we want to do make debug printout for this event?
      if ( !_debugEventList.inList() ) return;
      
      // Get information about this track.
      G4int id = track->GetTrackID();

      // If no tracks are listed, then printout for all tracks.
      // If some tracks are listed, then printout only for those tracks.
      if ( _debugTrackList.size() > 0 ){
          if ( !_debugTrackList.inList(id) ) return;
      }
      
      const G4DynamicParticle*  pParticle = track->GetDynamicParticle();
      double theKEnergy  = pParticle->GetKineticEnergy();
      const G4ThreeVector& theMomentumDirection = pParticle->GetMomentumDirection();
      G4int prec = G4cout.precision(15);
      G4cout << __func__ << " KE " << theKEnergy
             << " Momentum direction " << theMomentumDirection
             << G4endl;
      G4cout.precision(prec);

      // Momentum at the the pre point.
      G4ThreeVector const& mom = prept->GetMomentum();
      
      // On the last step on a track the post step point does not have an
      // associated physical volume. So we need to protect against that.
      G4String preVolume, postVolume;
      G4int preCopy(-1), postCopy(-1);
      
      // Get the names if they are defined.
      if ( prept->GetPhysicalVolume() ){
          preVolume = prept->GetPhysicalVolume()->GetName();
          preCopy   = prept->GetPhysicalVolume()->GetCopyNo();
      }
      
      if ( postpt->GetPhysicalVolume() ){
          postVolume = postpt->GetPhysicalVolume()->GetName();
          postCopy   = postpt->GetPhysicalVolume()->GetCopyNo();
      }
      
      // Status report. // fixme use G4cout to be compatible with MT
      printf ( "Step number: %d\n", numTrackSteps_ );
      
      printit ( "Pre: ", id,
               prept->GetPosition(),
               prept->GetMomentum(),
               prept->GetLocalTime(),
               prept->GetGlobalTime()
               );
      
      printit ( "Step:", id,
               track->GetPosition(),
               track->GetMomentum(),
               track->GetLocalTime(),
               track->GetGlobalTime()
               );
      
      printit ( "Post: ", id,
               postpt->GetPosition(),
               postpt->GetMomentum(),
               postpt->GetLocalTime(),
               postpt->GetGlobalTime()
               );
      fflush(stdout);
      
      G4cout << "Pre  Volume and copy: " << preVolume  << " " << preCopy  << G4endl;
      G4cout << "Post Volume and copy: " << postVolume << " " << postCopy << G4endl;
      
      printf ( "\n");
      fflush(stdout);
  
  } // end Mu2eG4SteppingAction


  // Kill tracks that take too many steps.
  bool Mu2eG4SteppingAction::killTooManySteps( const G4Track* const track ){
      if(numTrackSteps_ <= mu2elimits_->maxStepsPerTrack()) {
          return false;
      }
      
      if ( stepLimitKillerVerbose_ ) {
        G4cout << __func__ << " WARNING: kill particle in "
               << track->GetStep()->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetName()
               << " due to large number of steps." << G4endl;
        Mu2eG4UserHelpers::printKilledTrackInfo(track);
      }
      
      ++numKilledTracks_;
      return true;
  }

  // Record why the track is to be killed, then kill it.
  void Mu2eG4SteppingAction::killTrack( G4Track* track, ProcessCode::enum_type code, G4TrackStatus status ){
      
      // Get user track informaton object from the track.
      G4VUserTrackInformation* info = track->GetUserInformation();
      UserTrackInformation* tinfo   = static_cast<UserTrackInformation*>(info);
      
      // Record why the track was killed.
      tinfo->setProcessCode(ProcessCode(code));
      
      // Kill the track
      track->SetTrackStatus(status);
  }

  G4bool Mu2eG4SteppingAction::addTimeVDHit(const G4Step* aStep, int id){

      if(tvd_collection_->size() >= mu2elimits_->maxStepPointCollectionSize()) {
          if(!tvd_warning_printed_) {
              tvd_warning_printed_ = true;
              mf::LogWarning("G4") << "Maximum number of entries reached in time virtual detector: "
              << tvd_collection_->size() << endl;
          }
          return false;
      }
      
      // Which process caused this step to end?
      ProcessCode endCode(_processInfo->
                          findAndCount(Mu2eG4UserHelpers::findStepStoppingProcessName(aStep)));
      
      // The point's coordinates are saved in the mu2e coordinate system.
      tvd_collection_->
      push_back(StepPointMC(_spHelper->particlePtr(aStep->GetTrack()),
                            id,
                            0,
                            0,
                            aStep->GetPostStepPoint()->GetGlobalTime(),
                            aStep->GetPostStepPoint()->GetProperTime(),
                            aStep->GetPostStepPoint()->GetPosition() - _mu2eOrigin,
                            aStep->GetPostStepPoint()->GetMomentum(),
                            aStep->GetStepLength(),
                            endCode
                            ));
      return true;
  }


  std::vector<MCTrajectoryPoint> const& Mu2eG4SteppingAction::trajectory() {
      return _trajectory;
  }


  void  Mu2eG4SteppingAction::swapTrajectory(std::vector<MCTrajectoryPoint>& trajectory) {
      std::swap( trajectory, _trajectory);
  }

  double Mu2eG4SteppingAction::mcTrajectoryMinDistanceCut(const G4VPhysicalVolume* vol) const {
      
      const auto it = mcTrajectoryVolumePtDistances_.find(vol);
      return (it != mcTrajectoryVolumePtDistances_.end()) ? it->second : trajectoryControl_->defaultMinPointDistance();
  }

} // end namespace mu2e
