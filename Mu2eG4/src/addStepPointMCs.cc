//
// Add StepPointMC objects to the event.
//
// $Id: addStepPointMCs.cc,v 1.1 2010/04/07 23:19:57 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/04/07 23:19:57 $
//
// Original author Rob Kutschke
//
// Notes
// 1) At the creation of the geometry manager there was a check to
//    ensure that at most one tracker is created.
//

// Mu2e includes
#include "Mu2eG4/inc/addStepPointMCs.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "LTrackerGeom/inc/LTracker.hh"
#include "ITrackerGeom/inc/ITracker.hh"
#include "Mu2eG4/inc/StepPointG4.hh"

// G4 includes
#include "G4Event.hh"
#include "G4SDManager.hh"

namespace mu2e{

  // Public entry point: pick an tracker specific thing to do. See note 1).
  void addStepPointMCs( const G4Event* g4event, StepPointMCCollection& hits  ){

    edm::Service<GeometryService> geom;

    if ( geom->hasElement<LTracker>() ){
      addL( g4event, hits );
    }
    else if ( geom->hasElement<ITracker>() ) {
      addI( g4event, hits );
    }
  }

  // For the LTracker.
  void addL( const G4Event* g4event, StepPointMCCollection& outputHits  ){
    
    // G4 Hit collections for this event.
    G4HCofThisEvent* hce = g4event->GetHCofThisEvent();
    
    // Get the collection ID for the Sensitive Layer hits.
    // Magic name needs to move to a proper home.
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    G4int colId = SDman->GetCollectionID("StepPointG4Collection");
    
    if ( colId >= 0 && hce != 0 ){
      
      StepPointG4Collection* hits = static_cast<StepPointG4Collection*>(hce->GetHC(colId));
      G4int nHits = hits->entries();
        
      for (G4int i=0;i<nHits;i++) {
        StepPointG4* h = (*hits)[i];
        outputHits.push_back( h->hit() );
      }
    }
  } // end addL

  // For the ITracker.
  void addI( const G4Event* g4event, StepPointMCCollection& outputHits  ){

    // G4 Hit collections for this event.
    G4HCofThisEvent* hce = g4event->GetHCofThisEvent();

    // Get the collection ID for the Sensitive Layer hits.
    G4SDManager* SDman = G4SDManager::GetSDMpointer();

    GeomHandle<ITracker> itracker;
    std::stringstream SDname;
    std::string sSDname;
    std::string::size_type loc;
    for (int iSlr=0; iSlr < itracker->nSuperLayers(); ++iSlr) {
      for (int iRng = 0; iRng < itracker->nRing(); ++iRng) {
        SDname.str(std::string());
        SDname<< "StepPointG4Collection_";
        SDname<<"wvolS"<< std::setw(2) << iSlr << "R" << std::setw(2) <<iRng;
        sSDname=SDname.str();
        loc=sSDname.find( " ", 0 );
        while ( loc != std::string::npos ){
          sSDname.replace(loc,1,"0");
          loc = sSDname.find( " ", 0 );
        }
        //    			std::cout<<sSDname<<endl;
        G4int colId = SDman->GetCollectionID(sSDname.c_str());
        if ( colId >= 0 && hce != 0 ){
          StepPointG4Collection* hits = static_cast<StepPointG4Collection*>(hce->GetHC(colId));
          G4int nHits = hits->entries();

          for (G4int i=0;i<nHits;i++) {
            StepPointG4* h = (*hits)[i];
            outputHits.push_back( h->hit() );
          }

        }
        SDname.str(std::string());
        SDname<< "StepPointG4Collection_";
        SDname<<"gvolS"<< std::setw(2) << iSlr << "R" << std::setw(2) <<iRng;
        sSDname=SDname.str();
        loc = sSDname.find( " ", 0 );
        while ( loc != std::string::npos ){
          sSDname.replace(loc,1,"0");
          loc = sSDname.find( " ", 0 );
        }
        //    			std::cout<<sSDname<<endl;
        colId = SDman->GetCollectionID(sSDname.c_str());
        if ( colId >= 0 && hce != 0 ){
          StepPointG4Collection* hits = static_cast<StepPointG4Collection*>(hce->GetHC(colId));
          G4int nHits = hits->entries();

          for (G4int i=0;i<nHits;i++) {
            StepPointG4* h = (*hits)[i];
            outputHits.push_back( h->hit() );
          }
        }
      }
    }

  } // end addI

}


