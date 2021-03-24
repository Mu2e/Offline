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
// Author: 2007 Gunter Folger
//   created from HadronPhysicsFTFP
//
// Modified: KLG added Zhengyun's pbar related modifications (on the top of 9.6.p04)
//
//----------------------------------------------------------------------------
//
#if G4VERSION<4099
#include <iomanip>   

#include "Mu2eG4/inc/HadronPhysicsFTFP_BERT_PBAR_MU2E02.hh"

#include "Geant4/globals.hh"
#include "Geant4/G4ios.hh"
#include "Geant4/G4SystemOfUnits.hh"
#include "Geant4/G4ParticleDefinition.hh"
#include "Geant4/G4ParticleTable.hh"

#include "Geant4/G4MesonConstructor.hh"
#include "Geant4/G4BaryonConstructor.hh"
#include "Geant4/G4ShortLivedConstructor.hh"

#include "Geant4/G4ChipsKaonMinusInelasticXS.hh"
#include "Geant4/G4ChipsKaonPlusInelasticXS.hh"
#include "Geant4/G4ChipsKaonZeroInelasticXS.hh"
#include "Geant4/G4CrossSectionDataSetRegistry.hh"

#include "Geant4/G4PhysListUtil.hh"

// factory
#include "Geant4/G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(HadronPhysicsFTFP_BERT_PBAR_MU2E02);

HadronPhysicsFTFP_BERT_PBAR_MU2E02::HadronPhysicsFTFP_BERT_PBAR_MU2E02(G4int)
    :  G4VPhysicsConstructor("hInelastic FTFP_BERT_PBAR_MU2E02")
    , theNeutrons(0)
    , theBertiniNeutron(0)
    , theFTFPNeutron(0)
#if G4VERSION<4099
    , theLEPNeutron(0)
#endif
    , thePiK(0)
    , theBertiniPiK(0)
    , theFTFPPiK(0)
    , thePro(0)
    , theBertiniPro(0)
    , theFTFPPro(0)
    , theHyperon(0)
    , theAntiBaryon(0)
    , theFTFPAntiBaryon(0)
    , QuasiElastic(false)
    , ChipsKaonMinus(0)
    , ChipsKaonPlus(0)
    , ChipsKaonZero(0)
{}

HadronPhysicsFTFP_BERT_PBAR_MU2E02::HadronPhysicsFTFP_BERT_PBAR_MU2E02(const G4String& name, G4bool quasiElastic)
    :  G4VPhysicsConstructor(name) 
    , theNeutrons(0)
    , theBertiniNeutron(0)
    , theFTFPNeutron(0)
#if G4VERSION<4099
    , theLEPNeutron(0)
#endif
    , thePiK(0)
    , theBertiniPiK(0)
    , theFTFPPiK(0)
    , thePro(0)
    , theBertiniPro(0)
    , theFTFPPro(0)
    , theHyperon(0)
    , theAntiBaryon(0)
    , theFTFPAntiBaryon(0)
    , QuasiElastic(quasiElastic)
    , ChipsKaonMinus(0)
    , ChipsKaonPlus(0)
    , ChipsKaonZero(0)
{}

void HadronPhysicsFTFP_BERT_PBAR_MU2E02::CreateModels()
{

  theNeutrons=new G4NeutronBuilder;
  theFTFPNeutron=new G4FTFPNeutronBuilder(QuasiElastic);
  theNeutrons->RegisterMe(theFTFPNeutron);
  theNeutrons->RegisterMe(theBertiniNeutron=new G4BertiniNeutronBuilder);
  theBertiniNeutron->SetMinEnergy(0.0*GeV);
  theBertiniNeutron->SetMaxEnergy(5*GeV);
#if G4VERSION<4099
  theNeutrons->RegisterMe(theLEPNeutron=new G4LEPNeutronBuilder);
  theLEPNeutron->SetMinInelasticEnergy(0.0*eV);   // no inelastic from LEP
  theLEPNeutron->SetMaxInelasticEnergy(0.0*eV);  
#endif

  thePro=new G4ProtonBuilder;
  theFTFPPro=new PBARFTFPProtonBuilder(QuasiElastic);
  thePro->RegisterMe(theFTFPPro);
  thePro->RegisterMe(theBertiniPro=new G4BertiniProtonBuilder);
  theBertiniPro->SetMaxEnergy(5*GeV);

  thePiK=new G4PiKBuilder;
  theFTFPPiK=new G4FTFPPiKBuilder(QuasiElastic);
  thePiK->RegisterMe(theFTFPPiK);
  thePiK->RegisterMe(theBertiniPiK=new G4BertiniPiKBuilder);
  theBertiniPiK->SetMaxEnergy(5*GeV);
  
  theHyperon=new G4HyperonFTFPBuilder;
    
  theAntiBaryon=new G4AntiBarionBuilder;
  theAntiBaryon->RegisterMe(theFTFPAntiBaryon=new  G4FTFPAntiBarionBuilder(QuasiElastic));
}

HadronPhysicsFTFP_BERT_PBAR_MU2E02::~HadronPhysicsFTFP_BERT_PBAR_MU2E02()
{
  delete theNeutrons;
  delete theBertiniNeutron;
  delete theFTFPNeutron;
#if G4VERSION<4099
  delete theLEPNeutron;    
#endif

  delete thePiK;
  delete theBertiniPiK;
  delete theFTFPPiK;
    
  delete thePro;
  delete theBertiniPro;
  delete theFTFPPro;    
    
  delete theHyperon;
  delete theAntiBaryon;
  delete theFTFPAntiBaryon;
}

void HadronPhysicsFTFP_BERT_PBAR_MU2E02::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

#include "Geant4/G4ProcessManager.hh"
void HadronPhysicsFTFP_BERT_PBAR_MU2E02::ConstructProcess()
{
  CreateModels();
  theNeutrons->Build();
  thePro->Build();
  thePiK->Build();

  // use CHIPS cross sections also for Kaons
  ChipsKaonMinus = G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsKaonMinusInelasticXS::Default_Name());
  ChipsKaonPlus = G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsKaonPlusInelasticXS::Default_Name());
  ChipsKaonZero = G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsKaonZeroInelasticXS::Default_Name());
  //
    
  G4PhysListUtil::FindInelasticProcess(G4KaonMinus::KaonMinus())->AddDataSet(ChipsKaonMinus);
  G4PhysListUtil::FindInelasticProcess(G4KaonPlus::KaonPlus())->AddDataSet(ChipsKaonPlus);
  G4PhysListUtil::FindInelasticProcess(G4KaonZeroShort::KaonZeroShort())->AddDataSet(ChipsKaonZero );
  G4PhysListUtil::FindInelasticProcess(G4KaonZeroLong::KaonZeroLong())->AddDataSet(ChipsKaonZero );
    
  theHyperon->Build();
  theAntiBaryon->Build();
}
#endif
