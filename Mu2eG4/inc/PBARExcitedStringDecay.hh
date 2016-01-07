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
// Modified: KLG added Zhengyun's pbar related modifications (on the top of 9.6.p04)
//
#ifndef PBARExcitedStringDecay_h
#define PBARExcitedStringDecay_h 1

#include "globals.hh"
#include "G4VStringFragmentation.hh"
#include "G4ExcitedStringVector.hh"
#include "G4KineticTrackVector.hh"
#include "Mu2eG4/inc/PBARLundStringFragmentation.hh"

class PBARExcitedStringDecay: public G4VStringFragmentation 
{
  public:
      PBARExcitedStringDecay();
      PBARExcitedStringDecay(PBARVLongitudinalStringDecay * aStringDecay);
      virtual ~PBARExcitedStringDecay();

  private:
      PBARExcitedStringDecay(const PBARExcitedStringDecay &right);
      const PBARExcitedStringDecay & operator=(const PBARExcitedStringDecay &right);
      int operator==(const PBARExcitedStringDecay &right) const;
      int operator!=(const PBARExcitedStringDecay &right) const;

  public:

      virtual G4KineticTrackVector * FragmentStrings(const G4ExcitedStringVector * theStrings);

  private:
      G4KineticTrackVector * FragmentString(const G4ExcitedString &theString);
      G4bool EnergyAndMomentumCorrector(G4KineticTrackVector* Output, G4LorentzVector& TotalCollisionMom);   
  
      PBARVLongitudinalStringDecay * theStringDecay;

};

#endif


