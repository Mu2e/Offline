//
// A class to hold hits created by G4 in straw detectors.
// These hits are made by the sensitive detector class for straws, StrawSD.
// 
// $Id: StrawG4Hit.cc,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//

// Mu2e incldues
#include "Mu2eG4/inc/StrawG4Hit.hh"


// G4 includes 
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"



namespace mu2e {

  G4Allocator<StrawG4Hit> StrawG4HitAllocator;

  // I don't like this but that is the way G4 requires it to be defined.
  // Should compare equality of content not if they have the same address.
  G4int StrawG4Hit::operator==(const StrawG4Hit& right) const
  {
    return (this==&right) ? 1 : 0;
  }
  

  void StrawG4Hit::Draw(){
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    if(pVVisManager)
      {
	G4Circle circle(_hit.position());
	circle.SetScreenSize(2.);
	circle.SetFillStyle(G4Circle::filled);
	G4Colour colour(1.,0.,0.);
	G4VisAttributes attribs(colour);
	circle.SetVisAttributes(attribs);
	pVVisManager->Draw(circle);
      }
  }

  void StrawG4Hit::Print(){

    if ( _hit.trackId() != 1 ){

      // There is no units category for momentum, so use energy.
      G4cout << "  trackId: "        << _hit.trackId()
	     << "  strawIndex: "     << _hit.strawIndex()
	     << "  energy deposit: " << G4BestUnit(_hit.eDep(),"Energy")
	     << "  position: "       << G4BestUnit(_hit.position(),"Length") 
	     << "  momentum: "       << G4BestUnit(_hit.momentum(),"Energy") 
	     << "  time: "           << G4BestUnit(_hit.time(),"Time") 
	     << G4endl;
    }
  }

} // namespace mu2e
