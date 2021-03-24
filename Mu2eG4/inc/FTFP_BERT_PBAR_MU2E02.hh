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
//   created from FTFP
//
// Modified: KLG added Zhengyun's pbar related modifications (on the top of 9.6.p04)
//
//----------------------------------------------------------------------------
//
#ifndef TFTFP_BERT_PBAR_MU2E02_h
#define TFTFP_BERT_PBAR_MU2E02_h 1

#include <CLHEP/Units/SystemOfUnits.h>

#include "Geant4/globals.hh"
#include "Geant4/G4VModularPhysicsList.hh"
#include "CompileTimeConstraints.hh"

template<class T>
class TFTFP_BERT_PBAR_MU2E02: public T
{
public:
  TFTFP_BERT_PBAR_MU2E02(G4int ver = 1);
  virtual ~TFTFP_BERT_PBAR_MU2E02();
  
public:
  // SetCuts() 
  virtual void SetCuts();

private:
  enum {ok = CompileTimeConstraints::IsA<T, G4VModularPhysicsList>::ok };
};
#include "Mu2eG4/inc/FTFP_BERT_PBAR_MU2E02.icc"
typedef TFTFP_BERT_PBAR_MU2E02<G4VModularPhysicsList> FTFP_BERT_PBAR_MU2E02;

#endif
