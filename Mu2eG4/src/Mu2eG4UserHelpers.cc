//
// A collection of Geant4 user helper functions
// initially extracted from the TrackingAction
//
// $Id: Mu2eG4UserHelpers.cc,v 1.1 2012/12/19 23:37:13 genser Exp $
// $Author: genser $
// $Date: 2012/12/19 23:37:13 $
//
// Original author K. L. Genser based on Rob's TrackingAction
//
// Notes:

// c++ includes
#include <algorithm>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

// Mu2e includes.

#include "Mu2eG4/inc/Mu2eG4UserHelpers.hh"
#include "Mu2eG4/inc/UserTrackInformation.hh"
#include "MCDataProducts/inc/ProcessCode.hh"

// G4 includes

#include "G4Track.hh"
#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4VProcess.hh"
#include "G4RunManager.hh"
#include "G4TrackingManager.hh"

using namespace std;

namespace mu2e {

  namespace Mu2eG4UserHelpers {

    // Enable/disable storing of trajectories based on several considerations
    void controlTrajectorySaving(G4Track const* trk, int _sizeLimit, int _currentSize){

      // Do not add the trajectory if the corresponding SimParticle is missing.
      //    if( _sizeLimit>0 && _currentSize>_sizeLimit ) return;

      bool keep =  _sizeLimit>0 && _currentSize>_sizeLimit && saveThisTrajectory(trk);
      G4TrackingManager* trkmgr = G4EventManager::GetEventManager()->GetTrackingManager();
      if ( keep ) {
        trkmgr->SetStoreTrajectory(true);
      } else{
        trkmgr->SetStoreTrajectory(false);
      }

    }

    // The per track decision on whether or not to store trajectories.
    bool saveThisTrajectory(G4Track const* trk){

      // This is a guess at what might be useful.  Feel free to improve it.
      // We might want to change the momentum cut depending on which volume
      // the track starts in.
      CLHEP::Hep3Vector const& mom = trk->GetMomentum();
      bool keep = ( mom.mag() > 50.*CLHEP::MeV ); // fixme get it from config and make it static 

      return keep;
    }

    //Retrieve kinetic energy at the beginnig of the last step from UserTrackInfo
    double getPreLastStepKE(G4Track const* trk) {

      G4VUserTrackInformation* info = trk->GetUserInformation();
      UserTrackInformation const* tinfo   = dynamic_cast<UserTrackInformation*>(info);
      return tinfo->preLastStepKE();

    }

    // Get the number of G4 steps the track is made of
    int getNSteps(G4Track const* trk) {

      G4VUserTrackInformation* info = trk->GetUserInformation();
      UserTrackInformation const* tinfo   = dynamic_cast<UserTrackInformation*>(info);
      return tinfo->nSteps();

    }

    // Find the name of the code for the process that created this track.
    ProcessCode findCreationCode(G4Track const* trk){
      G4VProcess const* process = trk->GetCreatorProcess();

      // If there is no creator process, then the G4Track was created by PrimaryGenerator action.
      if ( process == 0 ){
        return ProcessCode::mu2ePrimary;
      }

      // Extract the name from the process and look up the code.
      string name = process->GetProcessName(); // may need to protect it
      return ProcessCode::findByName(name);
    }

    // Find the name of the process that stopped this track.
    G4String findStoppingProcessName(G4Track const* trk){

      // First check to see if Mu2e code killed this track.
      G4VUserTrackInformation* info = trk->GetUserInformation();
      UserTrackInformation const* tinfo   = dynamic_cast<UserTrackInformation*>(info);

      if ( tinfo->isForced() ){
        return tinfo->code().name();
      }

      // Otherwise, G4 killed this track.
      G4VProcess const* process = trk->GetStep()->GetPostStepPoint()->GetProcessDefinedStep();

      if (process) {
        return process->GetProcessName();
      } else {

        { // forcing mf own scope to prevent output interleaving; 
          // does not seem to work anyway
          mf::LogWarning("G4")
            << "ProcessDefinedStep NotSpecified for " 
            << trk->GetParticleDefinition()->GetParticleName()
            << endl;
        }

        static bool printItOnce = true;
      
        if (printItOnce) {

          // printItOnce = false;

          G4VProcess const* process = trk->GetStep()->GetPreStepPoint()->GetProcessDefinedStep();
          G4String pname = (process) ? process->GetProcessName() : "NotSpecified";

          cout << __func__ << " ProcessDefinedStep NotSpecified for " 
               << trk->GetParticleDefinition()->GetParticleName()
               << ", PostStepPoint StepStatus: "
               << trk->GetStep()->GetPostStepPoint()->GetStepStatus()
               << ", PreStepPoint  StepStatus: "
               << trk->GetStep()->GetPreStepPoint()->GetStepStatus()
               << ", PreStepPoint ProcessDefinedStep: "
               << pname
               << endl;

          int id       = trk->GetTrackID();
          int parentId = trk->GetParentID();

          cout << __func__ << " Track info:" 
               << " ID: " 
               << id
               << ", ParentID: "
               << parentId
               << ", CurrentStepNumber: "
               << trk->GetCurrentStepNumber()
               << ", TrackStatus: "
               << trk->GetTrackStatus()
               << ", TrackCreationCode: "
               << findCreationCode(trk)
               << endl;

          // the code commented out below is looking into the event past
          // and would require passing a lot of data around

          //         if ( parentId != 0 ){
          //           map_type::iterator i(_transientMap.find(key_type(parentId)));
          //           if ( i == _transientMap.end() ){
          //             throw cet::exception("RANGE")
          //               << "Could not find parent SimParticle in findStoppingProcess.  id: "
          //               << id
          //               << "\n";
          //           }

          //           int pdgId = i->second.pdgId();
          //           cout << __func__ << " Parent info:" 
          //                << " ID: " 
          //                << parentId
          //                << ", PDGiD: "
          //                << pdgId
          //                << ", name: "
          //                << G4ParticleTable::GetParticleTable()->FindParticle(pdgId)->GetParticleName()
          //                << ", endG4Status: "
          //                << i->second.endG4Status()
          //                << ", stoppingCode: "
          //                << i->second.stoppingCode()
          //                << endl;

        }

      }

      return G4String("NotSpecified");
    }

    void printTrackInfo(G4Track const* trk, std::string const& text,
                        map_type const& _transientMap,
                        art::CPUTimer const& _timer,
                        CLHEP::Hep3Vector const& _mu2eOrigin,
                        bool isEnd) {

      const G4Event* event = G4RunManager::GetRunManager()->GetCurrentEvent();

      // Get some properties of the tracks.
      G4VPhysicalVolume* pvol = trk->GetVolume();
      G4String volName = (pvol !=0) ?
        pvol->GetName(): "Unknown Volume";

      G4ParticleDefinition* pdef = trk->GetDefinition();
      G4String partName = (pdef !=0) ?
        pdef->GetParticleName() : "Unknown Particle";

      int id       = trk->GetTrackID();
      int parentId = trk->GetParentID();

      cout << text
           << setw(5) << event->GetEventID()   << " "
           << setw(4) << id                    << " "
           << setw(4) << parentId              << " "
           << setw(8) << partName              << " | "
           << trk->GetPosition()-_mu2eOrigin   << " "
           << trk->GetMomentum()               << " "
           << trk->GetKineticEnergy()          << " "
           << volName                          << " ";

      if ( isEnd ){
        cout << trk->GetProperTime() <<  " | ";
        map_type::const_iterator i(_transientMap.find(key_type(id)));
        if ( i != _transientMap.end() ){
          SimParticle const& particle = i->second;
          cout << particle.startGlobalTime() <<  " ";
        } else {
          cout << -1. <<  " ";
        }

        cout << trk->GetGlobalTime() << " | ";
        cout << _timer.cpuTime() << " "
             << _timer.realTime()
             << endl;
      }

      cout << endl;

    }

    bool checkCrossReferences( bool doPrint, bool doThrow, map_type const& _transientMap ){

      // Start by assuming we are ok; any error will turn this to false.
      bool ok(true);

      // Loop over all simulated particles.
      for ( map_type::const_iterator i=_transientMap.begin();
            i!=_transientMap.end(); ++i ){

        // The next particle to look at.
        SimParticle const& sim = i->second;

        key_type const & simid = sim.id();

        // Check that daughters point to the mother.
        std::vector<key_type> const& dau = sim.daughterIds();

        for ( std::vector<key_type>::const_iterator j=dau.begin();
              j!=dau.end(); ++j) {

          key_type parentId;

          map_type::const_iterator fdi = _transientMap.find(*j);
          bool daugterFound = fdi !=_transientMap.end(); 
          if (daugterFound) {
            parentId = (fdi->second).parentId();
          }

          if ( !daugterFound || parentId != simid ){

            // Daughter does not point back to the parent.
            ok = false;
            if ( doPrint ){
              mf::LogError("G4")
                << "StudyTrackingAction::checkCrossReferences: daughter does not point back to parent.\n";
            }
            if ( doThrow ){
              throw cet::exception("MU2EG4")
                << "StudyTrackingAction::checkCrossReferences: daughter does not point back to parent.\n";
            }
          }

        }

        // Check that this particle is in the list of its parent's daughters.
        if ( sim.hasParent() ){
          key_type parentId = sim.parentId();

          map_type::const_iterator fpi = _transientMap.find(parentId);
          bool parentFound = fpi != _transientMap.end();
        
          if ( !parentFound ){
            ok = false;
            if ( doPrint ){
              mf::LogError("G4")
                << "StudyTrackingAction::checkCrossReferences: parent not in list.\n";
            }
            if ( doThrow ){
              throw cet::exception("MU2EG4")
                << "StudyTrackingAction::checkCrossReferences: parent not in list.\n";
            }
          } else {

            std::vector<key_type> const& mdau = (fpi->second).daughterIds();
            bool inList(false);

            if (find(mdau.begin(), mdau.end(), simid)!=mdau.end()) {
              inList = true;
            }

            if ( !inList ){
              ok = false;
              if ( doPrint ){
                mf::LogError("G4")
                  << "StudyTrackingAction::checkCrossReferences: daughter is not found amoung parent's daughters.\n";
              }
              if ( doThrow ){
                throw cet::exception("MU2EG4")
                  << "StudyTrackingAction::checkCrossReferences: daughter is not found amoung parent's daughters.\n";
              }
            }

          }

        }

      }

      return ok;

    }

  }

}
