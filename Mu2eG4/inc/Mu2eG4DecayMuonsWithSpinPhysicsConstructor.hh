//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//---------------------------------------------------------------------------
//
// ClassName:   Mu2eG4DecayMuonsWithSpinPhysicsConstructor based on G4DecayPhysics & F05PhysicsList
//              a G4VPhysicsConstructor
// Author: KLG
//
//
//----------------------------------------------------------------------------
//

#ifndef Mu2eG4DecayMuonsWithSpinPhysicsConstructor_h
#define Mu2eG4DecayMuonsWithSpinPhysicsConstructor_h 1

#include "Geant4/globals.hh"
#include "Geant4/G4VPhysicsConstructor.hh"

#include "Geant4/G4DecayWithSpin.hh"

class Mu2eG4DecayMuonsWithSpinPhysicsConstructor : public G4VPhysicsConstructor
{
  public:
    Mu2eG4DecayMuonsWithSpinPhysicsConstructor(G4int ver = 1);
    Mu2eG4DecayMuonsWithSpinPhysicsConstructor(const G4String& name, G4int ver = 1);
    virtual ~Mu2eG4DecayMuonsWithSpinPhysicsConstructor() = default;
    Mu2eG4DecayMuonsWithSpinPhysicsConstructor(const Mu2eG4DecayMuonsWithSpinPhysicsConstructor &) = delete;
    Mu2eG4DecayMuonsWithSpinPhysicsConstructor & operator=(const Mu2eG4DecayMuonsWithSpinPhysicsConstructor &) = delete;
    Mu2eG4DecayMuonsWithSpinPhysicsConstructor & operator=( Mu2eG4DecayMuonsWithSpinPhysicsConstructor && ) = delete;
  public:
    // This method will be invoked in the Construct() method.
    // each particle type will be instantiated
  virtual void ConstructParticle();

    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type
  virtual void ConstructProcess();

private:

  G4int    verbose;
};


#endif
