//
// A class to hold hits created by G4 in most sensitive detectors.
// 
// $Id: StepPointG4.cc,v 1.2 2010/04/04 20:35:13 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/04/04 20:35:13 $
//
// Original author Rob Kutschke
//

// Mu2e incldues
#include "Mu2eG4/inc/StepPointG4.hh"

// G4 includes 
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

namespace mu2e {

  G4Allocator<StepPointG4> StepPointG4Allocator;

  // G4 defines this operator to compare by address, not to compare content.
  G4int StepPointG4::operator==(const StepPointG4& right) const
  {
    return (this==&right) ? 1 : 0;
  }
  

  void StepPointG4::Draw(){
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

  void StepPointG4::Print(){

    if ( _hit.trackId() != 1 ){

      // There is no units category for momentum, so use energy.
      G4cout << "  trackId: "        << _hit.trackId()
             << "  volumeId: "       << _hit.volumeId()
             << "  energy deposit: " << G4BestUnit(_hit.eDep(),"Energy")
             << "  position: "       << G4BestUnit(_hit.position(),"Length") 
             << "  momentum: "       << G4BestUnit(_hit.momentum(),"Energy") 
             << "  time: "           << G4BestUnit(_hit.time(),"Time") 
             << "  step Length: "    << G4BestUnit(_hit.stepLength(),"Length")
             << G4endl;
    }
  }

} // namespace mu2e
