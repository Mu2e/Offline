// Called at every G4 step.
//
// Original author Rob Kutschke
//

// C++ includes
#include <cstdio>
#include <cmath>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

// Mu2e includes
#include "Mu2eG4/inc/Mu2eG4SteppingAction.hh"
#include "Mu2eG4/inc/Mu2eG4UserHelpers.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "MCDataProducts/inc/PDGCode.hh"
#include "Mu2eG4/inc/getPhysicalVolumeOrThrow.hh"
#include "Mu2eG4/inc/SimParticleHelper.hh"

// G4 includes
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4SteppingManager.hh"
#include "G4String.hh"

using namespace std;

namespace mu2e {

  Mu2eG4SteppingAction::Mu2eG4SteppingAction(const fhicl::ParameterSet& pset, IMu2eG4Cut& cuts) :
    pset_(pset),

    steppingCuts_(&cuts),

    maxStepsPerTrack_(pset.get<int>("ResourceLimits.maxStepsPerTrack")),
    numTrackSteps_(),
    numKilledTracks_(),
    stepLimitKillerVerbose_(pset.get<bool>("debug.stepLimitKillerVerbose")),

    stepPointCollectionSizeLimit_(pset.get<StepPointMCCollection::size_type>("ResourceLimits.maxStepPointCollectionSize")),

    // Things related to time virtual detector
    tvd_time_(pset.get<vector<double> >("TimeVD", vector<double>())),
    tvd_collection_(nullptr),
    tvd_warning_printed_(false),

    mcTrajectoryDefaultMinPointDistance_(pset.get<double>("TrajectoryControl.defaultMinPointDistance")),

    // Default values for parameters that are optional in the run time configuration.
    _debugEventList(pset.get<std::vector<int> >("debug.eventList", std::vector<int>())),
    _debugTrackList(pset.get<std::vector<int> >("debug.trackList", std::vector<int>())),

    _spHelper()
  {
    if(maxStepsPerTrack_ > 0) {
      cout << "Limit maximum number of steps per track in Mu2eG4SteppingAction to "
           << maxStepsPerTrack_ << endl;
    }

    if(!tvd_time_.empty()) {
      cout << "Time virtual detector is enabled. Particles are recorded at";
      for( unsigned int i=0; i<tvd_time_.size(); ++i ) cout << " " << tvd_time_[i];
      cout << " ns" << endl;
    }
  }

  // A helper function to manage the printout.
  void Mu2eG4SteppingAction::printit( G4String const& s,
                                      G4int id,
                                      G4ThreeVector const& pos,
                                      G4ThreeVector const& mom,
                                      double localTime,
                                      double globalTime ){

    // It is easier to line up printout in columns with printf than with cout.
    printf ( "%-8s %4d %15.4f %15.4f %15.4f %15.4f %15.4f %15.4f %15.4f %13.4f %13.4f\n",
             s.data(), id,
             pos.x(), pos.y(), pos.z(),
             mom.x(), mom.y(), mom.z(),
             mom.mag(),
             localTime, globalTime);
  }

  void Mu2eG4SteppingAction::beginRun(PhysicsProcessInfo& processInfo,
                                      CLHEP::Hep3Vector const& mu2eOrigin ){
    _processInfo    = &processInfo;
    _mu2eOrigin     =  mu2eOrigin;
  }

  void Mu2eG4SteppingAction::finishConstruction() {
    // We have to wait until G4 geometry is constructed
    // to get phys volume pointers that are used in the
    // volume to cut value map.
    std::cout<<"Mu2eG4SteppingAction: mcTrajectoryDefaultMinPointDistance = "<<mcTrajectoryDefaultMinPointDistance_<<std::endl;
    const fhicl::ParameterSet& volumeCutsPS{pset_.get<fhicl::ParameterSet>("TrajectoryControl.perVolumeMinDistance")};
    const std::vector<std::string> volnames{volumeCutsPS.get_keys()};
    for(const auto& k: volnames) {
      auto vol = getPhysicalVolumeOrThrow(k);
      double cut = mcTrajectoryVolumePtDistances_[vol] = volumeCutsPS.get<double>(k);
      std::cout<<"Mu2eG4SteppingAction: mcTrajectory distance cut = "
               <<cut<<" for volume "<<k<<std::endl;
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
    if(!computeMCTrajDistance || ((prept->GetPosition() -  _trajectory.back().vect()).mag() >= mcTrajCurrentCut)) {
      _trajectory.emplace_back ( prept->GetPosition(), prept->GetGlobalTime() );
    }

    // Get kinetic energy at the begin of the step
    const double preStepEK = prept->GetKineticEnergy();

    G4VUserTrackInformation* info = track->GetUserInformation();
    UserTrackInformation* tinfo   = static_cast<UserTrackInformation*>(info);

    tinfo->setStepInfo(preStepEK, numTrackSteps_);

    // Save hits in time virtual detector
    for( unsigned int i=0; i<tvd_time_.size(); ++i ) {
      if( prept->GetGlobalTime()<=tvd_time_[i] && postpt->GetGlobalTime()>tvd_time_[i] ) {
        addTimeVDHit(step,i+1);
      }
    }

    if(steppingCuts_->steppingActionCut(step)) {
      // FIXME: do we need to differentiate stopping codes?
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

    // Position and momentum at the the pre point.
    G4ThreeVector const& pos = prept->GetPosition();
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

    // Status report.
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

    cout << "Pre  Volume and copy: " << preVolume  << " " << preCopy  << endl;
    cout << "Post Volume and copy: " << postVolume << " " << postCopy << endl;

    printf ( "\n");
    fflush(stdout);

  } // end Mu2eG4SteppingAction


  // Kill tracks that take too many steps.
  bool Mu2eG4SteppingAction::killTooManySteps( const G4Track* track ){

    if(maxStepsPerTrack_ > 0 && numTrackSteps_ <= maxStepsPerTrack_ ) {
      return false;
    }

    if ( stepLimitKillerVerbose_ ) {
      cout << "Mu2eG4SteppingAction: kill particle pdg="
           << track->GetDefinition()->GetPDGEncoding()
           << " due to large number of steps." << endl;
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

    if( tvd_collection_->size() >= stepPointCollectionSizeLimit_ ) {
      if(!tvd_warning_printed_) {
        tvd_warning_printed_ = true;
        mf::LogWarning("G4") << "Maximum number of entries reached in time virtual detector "
                             << stepPointCollectionSizeLimit_ << endl;
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

  std::vector<CLHEP::HepLorentzVector> const& Mu2eG4SteppingAction::trajectory() {
    return _trajectory;
  }

  void  Mu2eG4SteppingAction::swapTrajectory(std::vector<CLHEP::HepLorentzVector>& trajectory) {
    std::swap( trajectory, _trajectory);
  }

  double Mu2eG4SteppingAction::mcTrajectoryMinDistanceCut(const G4VPhysicalVolume* vol) const {
    const auto it = mcTrajectoryVolumePtDistances_.find(vol);
    return (it != mcTrajectoryVolumePtDistances_.end()) ? it->second : mcTrajectoryDefaultMinPointDistance_;
  }

  //================================================================
  void Mu2eG4SteppingAction::checkConfigRelics(const SimpleConfig& config) {
    static const std::vector<std::string> keys = {
      "g4.steppingActionStepsSizeLimit",
      "g4.killLowEKine",
      "g4SteppingAction.killInTheseVolumes",
      "g4SteppingAction.killerVerbose",
      "g4SteppingAction.killInHallAir",
      "g4.eKineMin",
      "g4.killLowEKinePDG",
      "g4.eKineMinPDG",
      "g4.steppingActionEventList",
      "g4.steppingActionTrackList",
      "g4.steppingActionMaxSteps",
      "g4.steppingActionMaxGlobalTime",
      "g4.steppingActionTimeVD",
      "g4.mcTrajectoryVolumes",
      "g4.mcTrajectoryVolumePtDistances",
      "g4.mcTrajectoryDefaultMinPointDistance"
    };

    std::string present;
    for(const auto k: keys) {
      if(config.hasName(k)) {
        present += k+" ";
      }
    }
    if(!present.empty()) {
      throw cet::exception("CONFIG")<<"Mu2eG4SteppingAction: Please use fcl to configure Mu2eG4_module. "
                                    <<"Detected obsolete SimpleConfig parameters: "<<present;
    }
  }

} // end namespace mu2e
