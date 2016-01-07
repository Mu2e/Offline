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
// $Id$

#ifndef PBARDiffractiveExcitation_h
#define PBARDiffractiveExcitation_h 1
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      ---------------- PBARDiffractiveExcitation --------------
//             by Gunter Folger, October 1998.
//      diffractive Excitation used by strings models
//	Take a projectile and a target
//	excite the projectile and target
// Modified: KLG added Zhengyun's pbar related modifications (on the top of 9.6.p04)
// ------------------------------------------------------------

#include "globals.hh"
class G4VSplitableHadron;
class G4ExcitedString;
#include "G4FTFParameters.hh"
#include "G4ElasticHNScattering.hh"  // Uzhi 3.09.09
#include "G4ThreeVector.hh"

class PBARDiffractiveExcitation 
{

  public:

      PBARDiffractiveExcitation();
      virtual ~PBARDiffractiveExcitation();

      virtual G4bool ExciteParticipants (G4VSplitableHadron *aPartner, 
                                         G4VSplitableHadron * bPartner,
                                         G4FTFParameters *theParameters,
                                         G4ElasticHNScattering *theElastic) const;

      virtual void CreateStrings        (G4VSplitableHadron * aHadron, 
                                         G4bool isProjectile,
                                         G4ExcitedString * &FirstString, 
                                         G4ExcitedString * &SecondString,
                                         G4FTFParameters *theParameters) const;
  private:

      PBARDiffractiveExcitation(const PBARDiffractiveExcitation &right);
      
      G4ThreeVector GaussianPt(G4double  AveragePt2, G4double maxPtSquare) const;
      G4double ChooseP(G4double Pmin, G4double Pmax) const;
      G4double GetQuarkFractionOfKink(G4double zmin, G4double zmax) const;
      void UnpackMeson( G4int IdPDG, G4int &Q1, G4int &Q2) const;  // Uzhi 7.09.09
      void UnpackBaryon(G4int IdPDG, G4int &Q1, G4int &Q2, G4int &Q3) const; // Uzhi 7.09.09
      G4int NewNucleonId(G4int Q1, G4int Q2, G4int Q3) const; // Uzhi 7.09.09

      const PBARDiffractiveExcitation & operator=(const PBARDiffractiveExcitation &right);
      int operator==(const PBARDiffractiveExcitation &right) const;
      int operator!=(const PBARDiffractiveExcitation &right) const;

};

#endif
