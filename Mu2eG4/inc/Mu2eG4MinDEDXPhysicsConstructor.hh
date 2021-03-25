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
//---------------------------------------------------------------------------
//
// ClassName:   Mu2eG4MinDEDXPhysicsConstructor
//
// Author:      KLG based on G4EmStandardPhysics
//
//----------------------------------------------------------------------------
//
// This class constructs special limited EM physics
//

#ifndef Mu2eG4MinDEDXPhysicsConstructor_h
#define Mu2eG4MinDEDXPhysicsConstructor_h 1

#include "Geant4/G4VPhysicsConstructor.hh"
#include "Geant4/globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Mu2eG4MinDEDXPhysicsConstructor : public G4VPhysicsConstructor
{
public:
  Mu2eG4MinDEDXPhysicsConstructor(G4int ver = 0);

  virtual ~Mu2eG4MinDEDXPhysicsConstructor();

  virtual void ConstructParticle();
  virtual void ConstructProcess();

private:
  G4int  verbose;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
