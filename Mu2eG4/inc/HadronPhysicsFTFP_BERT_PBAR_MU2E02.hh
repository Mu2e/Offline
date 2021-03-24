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
// $Id$
//
//---------------------------------------------------------------------------
//
// ClassName:   
//
// Author: 2007  Gunter Folger
//   created from HadronPhysicsFTFP
// Modified:
// 23.11.2005 G.Folger: migration to non static particles
// 08.06.2006 V.Ivanchenko: remove stopping
// 19.06.2008 G.Folger: change default for QE to NOT use Chips QE
// Modified: KLG added Zhengyun's pbar related modifications (on the top of 9.6.p04)
//
//----------------------------------------------------------------------------
//
#ifndef HadronPhysicsFTFP_BERT_PBAR_MU2E02_h
#define HadronPhysicsFTFP_BERT_PBAR_MU2E02_h 1

#include "Geant4/globals.hh"
#include "Geant4/G4ios.hh"

#include "Geant4/G4VPhysicsConstructor.hh"

#include "Geant4/G4PiKBuilder.hh"
#include "Geant4/G4BertiniPiKBuilder.hh"
#include "Geant4/G4FTFPPiKBuilder.hh"

#include "Geant4/G4ProtonBuilder.hh"
#include "Geant4/G4BertiniProtonBuilder.hh"
#include "Geant4/G4FTFPNeutronBuilder.hh"
#include "Mu2eG4/inc/PBARFTFPProtonBuilder.hh"

#include "Geant4/G4NeutronBuilder.hh"
#include "Geant4/G4BertiniNeutronBuilder.hh"
#include "Geant4/G4FTFPNeutronBuilder.hh"
#if G4VERSION<4099
#include "Geant4/G4LEPNeutronBuilder.hh"
#endif
#include "Geant4/G4HyperonFTFPBuilder.hh"
#include "Geant4/G4AntiBarionBuilder.hh"
#include "Geant4/G4FTFPAntiBarionBuilder.hh"

class HadronPhysicsFTFP_BERT_PBAR_MU2E02 : public G4VPhysicsConstructor
{
  public: 
    HadronPhysicsFTFP_BERT_PBAR_MU2E02(G4int verbose =1);
    HadronPhysicsFTFP_BERT_PBAR_MU2E02(const G4String& name, G4bool quasiElastic=false);
    virtual ~HadronPhysicsFTFP_BERT_PBAR_MU2E02();

  public: 
    virtual void ConstructParticle();
    virtual void ConstructProcess();

  private:
    void CreateModels();
    
    G4NeutronBuilder * theNeutrons;
    G4BertiniNeutronBuilder * theBertiniNeutron;
    G4FTFPNeutronBuilder * theFTFPNeutron;
#if G4VERSION<4099
    G4LEPNeutronBuilder * theLEPNeutron;        //needed for capture&fission
#endif 
    G4PiKBuilder * thePiK;
    G4BertiniPiKBuilder * theBertiniPiK;
    G4FTFPPiKBuilder * theFTFPPiK;
    
    G4ProtonBuilder * thePro;
    G4BertiniProtonBuilder * theBertiniPro;
    PBARFTFPProtonBuilder * theFTFPPro;    
    
    G4HyperonFTFPBuilder * theHyperon;
    
    G4AntiBarionBuilder * theAntiBaryon;
    G4FTFPAntiBarionBuilder * theFTFPAntiBaryon;

    G4bool QuasiElastic;
    G4VCrossSectionDataSet * ChipsKaonMinus;
    G4VCrossSectionDataSet * ChipsKaonPlus;
    G4VCrossSectionDataSet * ChipsKaonZero;
    G4VCrossSectionDataSet * BGGProton;
    G4VCrossSectionDataSet * BGGNeutron;
    
};

#endif

