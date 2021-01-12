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
// Author: 2010 Tatsumi Koi, Gunter Folger
//   created from HadronPhysicsFTFP_BERT
//
// Modified:
// 30 06 2013 K.Genser Mu2e version with BERT FTPF transition between 9.5 and 9.9 GeV
//
//----------------------------------------------------------------------------
//
#include <iomanip>

#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Mu2eG4/inc/HadronPhysicsShielding_MU2E00.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

#include "G4BGGNucleonInelasticXS.hh"
#include "G4NeutronHPBGGNucleonInelasticXS.hh"
#include "G4NeutronHPJENDLHEInelasticData.hh"
#include "G4NeutronHPInelasticData.hh"

#include "G4ChipsKaonMinusInelasticXS.hh"
#include "G4ChipsKaonPlusInelasticXS.hh"
#include "G4ChipsKaonZeroInelasticXS.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4PhysListUtil.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(HadronPhysicsShielding_MU2E00);

HadronPhysicsShielding_MU2E00::HadronPhysicsShielding_MU2E00( G4int )
    :  G4VPhysicsConstructor("hInelastic Shielding")
    , theNeutrons(0)
    , theLENeutron(0)
    , theBertiniNeutron(0)
    , theFTFPNeutron(0)
    , theLEPNeutron(0)
    , thePiK(0)
    , theBertiniPiK(0)
    , theFTFPPiK(0)
    , thePro(0)
    , theBertiniPro(0)
    , theFTFPPro(0)
#if G4VERSION<4099
    , theMiscCHIPS(0)
#endif
    , QuasiElastic(false)
    , theCHIPSInelastic(0)
    , BGGxsNeutron(0)
    , NeutronHPJENDLHEInelastic(0)
    , BGGxsProton(0)
    , useLEND(false)
    , evaluation()
{}

HadronPhysicsShielding_MU2E00::HadronPhysicsShielding_MU2E00(const G4String& name, G4bool quasiElastic)
    :  G4VPhysicsConstructor(name)
    , theNeutrons(0)
    , theLENeutron(0)
    , theBertiniNeutron(0)
    , theFTFPNeutron(0)
    , theLEPNeutron(0)
    , thePiK(0)
    , theBertiniPiK(0)
    , theFTFPPiK(0)
    , thePro(0)
    , theBertiniPro(0)
    , theFTFPPro(0)
#if G4VERSION<4099
    , theMiscCHIPS(0)
#endif
    , QuasiElastic(quasiElastic)
    , theCHIPSInelastic(0)
    , BGGxsNeutron(0)
    , NeutronHPJENDLHEInelastic(0)
    , BGGxsProton(0)
    , useLEND(false)
    , evaluation()
{}

#include "G4NeutronLENDBuilder.hh"
void HadronPhysicsShielding_MU2E00::CreateModels()
{

#if G4VERSION>4099
  mf::LogError("PHYS") << " This Mu2e Physics Constructor has not been certified for use with Geant4 v10+.";
  G4cout << "Error: This Mu2e Physics Constructor has not been certified for use with Geant4 v10+." << G4endl;
  throw cet::exception("BADINPUT")<<"This Mu2e Physics Constructor has not been certified for use with Geant4 v10+.\n";
#endif

  const G4double minFTFPEnergy         =  9.5*GeV;
  const G4double maxBertiniEnergy      =  9.9*GeV;
  const G4double minNonHPNeutronEnergy = 19.9*MeV;

  theNeutrons=new G4NeutronBuilder;
  theFTFPNeutron=new G4FTFPNeutronBuilder(QuasiElastic);
  theFTFPNeutron->SetMinEnergy(minFTFPEnergy);
  theNeutrons->RegisterMe(theFTFPNeutron);
  theNeutrons->RegisterMe(theBertiniNeutron=new G4BertiniNeutronBuilder);
  theBertiniNeutron->SetMinEnergy(minNonHPNeutronEnergy);
  theBertiniNeutron->SetMaxEnergy(maxBertiniEnergy);
#if G4VERSION<4099
  theNeutrons->RegisterMe(theLEPNeutron=new G4LEPNeutronBuilder);
  theLEPNeutron->SetMinEnergy(minNonHPNeutronEnergy);
  theLEPNeutron->SetMinInelasticEnergy(0.0*eV);   // no inelastic from LEP
  theLEPNeutron->SetMaxInelasticEnergy(0.0*eV);
  //theNeutrons->RegisterMe(theHPNeutron=new G4NeutronHPBuilder);
#endif
    if ( useLEND != true )
     theNeutrons->RegisterMe(theLENeutron=new G4NeutronHPBuilder);
  else
  {
     theNeutrons->RegisterMe(theLENeutron=new G4NeutronLENDBuilder(evaluation));
  }

  thePro=new G4ProtonBuilder;
  theFTFPPro=new G4FTFPProtonBuilder(QuasiElastic);
  theFTFPPro->SetMinEnergy(minFTFPEnergy);
  thePro->RegisterMe(theFTFPPro);
  thePro->RegisterMe(theBertiniPro=new G4BertiniProtonBuilder);
  theBertiniPro->SetMaxEnergy(maxBertiniEnergy);

  thePiK=new G4PiKBuilder;
  theFTFPPiK=new G4FTFPPiKBuilder(QuasiElastic);
  theFTFPPiK->SetMinEnergy(minFTFPEnergy);
  thePiK->RegisterMe(theFTFPPiK);
  thePiK->RegisterMe(theBertiniPiK=new G4BertiniPiKBuilder);
  theBertiniPiK->SetMaxEnergy(maxBertiniEnergy);

#if G4VERSION<4099
  theMiscCHIPS=new G4MiscCHIPSBuilder;
#endif

}

HadronPhysicsShielding_MU2E00::~HadronPhysicsShielding_MU2E00()
{
  delete theNeutrons;
  delete theBertiniNeutron;
  delete theFTFPNeutron;
  //delete theHPNeutron;
  delete theLENeutron;

  delete thePiK;
  delete theBertiniPiK;
  delete theFTFPPiK;

  delete thePro;
  delete theBertiniPro;
  delete theFTFPPro;

#if G4VERSION<4099
  delete theMiscCHIPS;
#endif
  delete theCHIPSInelastic;
  delete BGGxsNeutron;
  delete NeutronHPJENDLHEInelastic;
  delete BGGxsProton;
}

void HadronPhysicsShielding_MU2E00::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();
}

#include "G4ProcessManager.hh"
void HadronPhysicsShielding_MU2E00::ConstructProcess()
{
  CreateModels();

  //BGGxsNeutron=new  G4BGGNucleonInelasticXS(G4Neutron::Neutron());
  thePro->Build();
  theNeutrons->Build();

//  BGGxsNeutron=new  G4NeutronHPBGGNucleonInelasticXS(G4Neutron::Neutron());
//  FindInelasticProcess(G4Neutron::Neutron())->AddDataSet(BGGxsNeutron);
//

  G4PhysListUtil::FindInelasticProcess(G4Neutron::Neutron())->AddDataSet(new  G4BGGNucleonInelasticXS(G4Neutron::Neutron()));
  NeutronHPJENDLHEInelastic=new G4NeutronHPJENDLHEInelasticData;
  G4PhysListUtil::FindInelasticProcess(G4Neutron::Neutron())->AddDataSet(NeutronHPJENDLHEInelastic);
  G4PhysListUtil::FindInelasticProcess(G4Neutron::Neutron())->AddDataSet(new G4NeutronHPInelasticData);

//  BGGxsProton=new G4BGGNucleonInelasticXS(G4Proton::Proton());
//  G4PhysListUtil::FindInelasticProcess(G4Proton::Proton())->AddDataSet(BGGxsProton);

  thePiK->Build();
  // use CHIPS cross sections also for Kaons

  G4CrossSectionDataSetRegistry* registry = G4CrossSectionDataSetRegistry::Instance();

  G4PhysListUtil::FindInelasticProcess(G4KaonMinus::KaonMinus())->AddDataSet(registry->GetCrossSectionDataSet(G4ChipsKaonMinusInelasticXS::Default_Name()));
  G4PhysListUtil::FindInelasticProcess(G4KaonPlus::KaonPlus())->AddDataSet(registry->GetCrossSectionDataSet(G4ChipsKaonPlusInelasticXS::Default_Name()));
  G4PhysListUtil::FindInelasticProcess(G4KaonZeroShort::KaonZeroShort())->AddDataSet(registry->GetCrossSectionDataSet(G4ChipsKaonZeroInelasticXS::Default_Name()));
  G4PhysListUtil::FindInelasticProcess(G4KaonZeroLong::KaonZeroLong())->AddDataSet(registry->GetCrossSectionDataSet(G4ChipsKaonZeroInelasticXS::Default_Name()));

#if G4VERSION<4099
  theMiscCHIPS->Build();
#endif

}
