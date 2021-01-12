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
// ClassName:   
//
// Author: 2007  tatsumi Koi, Gunter Folger
//   created from HadronPhysicsFTFP_BERT
// Modified:
//
//----------------------------------------------------------------------------
//
#ifndef HadronPhysicsShielding_MU2E02_h
#define HadronPhysicsShielding_MU2E02_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"

#include "G4PiKBuilder.hh"
#include "G4BertiniPiKBuilder.hh"
#include "G4FTFPPiKBuilder.hh"

#include "G4ProtonBuilder.hh"
#include "G4BertiniProtonBuilder.hh"
#include "G4FTFPNeutronBuilder.hh"
#include "G4FTFPProtonBuilder.hh"

#include "G4NeutronBuilder.hh"
#include "G4BertiniNeutronBuilder.hh"
#include "G4FTFPNeutronBuilder.hh"
#if G4VERSION<4099
#include "G4LEPNeutronBuilder.hh"
#endif
#include "G4NeutronHPBuilder.hh"

#include "G4HyperonFTFPBuilder.hh"
#include "G4AntiBarionBuilder.hh"
#include "G4FTFPAntiBarionBuilder.hh"

class HadronPhysicsShielding_MU2E02 : public G4VPhysicsConstructor
{
  public: 
    //HadronPhysicsShielding_MU2E02(G4int verbose =1,G4bool blend=false);
    HadronPhysicsShielding_MU2E02(G4int verbose =1);
    HadronPhysicsShielding_MU2E02(const G4String& name, G4bool quasiElastic=false);
    virtual ~HadronPhysicsShielding_MU2E02();

  public: 
    virtual void ConstructParticle();
    virtual void ConstructProcess();
    void UseLEND( G4String ss="" ){useLEND=true;evaluation=ss;};
    void UnuseLEND(){useLEND=false;};

  private:
    void CreateModels();
    
    G4NeutronBuilder * theNeutrons;
  // G4NeutronHPBuilder * theHPNeutron; // it is invoked for theNeutrons in cc file
    G4VNeutronBuilder * theLENeutron;
    G4BertiniNeutronBuilder * theBertiniNeutron;
    G4FTFPNeutronBuilder * theFTFPNeutron;
#if G4VERSION<4099
    G4LEPNeutronBuilder * theLEPNeutron;        //needed for capture&fission;
#else
    void* theLEPNeutron; // fake
#endif

    G4PiKBuilder * thePiK;
    G4BertiniPiKBuilder * theBertiniPiK;
    G4FTFPPiKBuilder * theFTFPPiK;
    
    G4ProtonBuilder * thePro;
    G4BertiniProtonBuilder * theBertiniPro;
    G4FTFPProtonBuilder * theFTFPPro;    

    G4HyperonFTFPBuilder * theHyperon;
    
    G4AntiBarionBuilder * theAntiBaryon;
    G4FTFPAntiBarionBuilder * theFTFPAntiBaryon;

    G4VCrossSectionDataSet * theChipsKaonMinus;
    G4VCrossSectionDataSet * theChipsKaonPlus;
    G4VCrossSectionDataSet * theChipsKaonZero;
 
    G4bool QuasiElastic;

    G4VCrossSectionDataSet * NeutronHPJENDLHEInelastic;

    G4bool useLEND;
    G4String evaluation;
};

#endif

