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
// ClassName:   PBARFTFPProtonBuilder
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 18.11.2010 G.Folger, use G4CrossSectionPairGG for relativistic rise of cross
//             section at high energies.
// 30.03.2009 V.Ivanchenko create cross section by new
// Modified: KLG added Zhengyun's pbar related modifications (on the top of 9.6.p04)
//
//----------------------------------------------------------------------------
//
#if G4VERSION<4099
#include "Mu2eG4/inc/PBARFTFPProtonBuilder.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4BGGNucleonInelasticXS.hh"

PBARFTFPProtonBuilder::
PBARFTFPProtonBuilder(G4bool quasiElastic) 
{
  theMin = 4.*GeV;
  theMax = 100.*TeV; 
  theModel = new G4TheoFSGenerator("FTFP");

  theStringModel = new PBARFTFModel();
  theStringDecay = new PBARExcitedStringDecay(theLund = new PBARLundStringFragmentation);
  theStringModel->SetFragmentationModel(theStringDecay);

  thePreEquilib = new G4PreCompoundModel(theHandler = new G4ExcitationHandler);
  theCascade = new G4GeneratorPrecompoundInterface(thePreEquilib);

  theModel->SetHighEnergyGenerator(theStringModel);
  if (quasiElastic)
  {
     theQuasiElastic=new G4QuasiElasticChannel;
     theModel->SetQuasiElasticChannel(theQuasiElastic);
  } else 
  {  theQuasiElastic=0;}  

  theModel->SetTransport(theCascade);
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(100*TeV);
}

void PBARFTFPProtonBuilder::
Build(G4ProtonInelasticProcess * aP)
{
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(theMax);
  aP->RegisterMe(theModel);
    
    aP->AddDataSet(new G4BGGNucleonInelasticXS(G4Proton::Proton()));
}

PBARFTFPProtonBuilder::
~PBARFTFPProtonBuilder() 
{
  delete theStringDecay;
  delete theStringModel;
  delete theModel;
  delete theCascade;
  if ( theQuasiElastic ) delete theQuasiElastic;
  //delete theHandler;
  delete theLund;
}

void PBARFTFPProtonBuilder::
Build(G4HadronElasticProcess * )
{
}

 // 2002 by J.P. Wellisch
#endif
