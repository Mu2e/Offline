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
// -----------------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, Maxim Komogorov, 1-Jul-1998
//               redesign  Gunter Folger, August/September 2001
// -----------------------------------------------------------------------------

#if G4VERSION<4099

#include "Mu2eG4/inc/PBARVLongitudinalStringDecay.hh"
#include "Geant4/G4PhysicalConstants.hh"
#include "Geant4/G4SystemOfUnits.hh"
#include "Geant4/G4ios.hh"
#include "Randomize.hh"
#include "Geant4/G4FragmentingString.hh"

#include "Geant4/G4ParticleDefinition.hh"
#include "Geant4/G4ParticleTypes.hh"
#include "Geant4/G4ParticleChange.hh"
#include "Geant4/G4VShortLivedParticle.hh"
#include "Geant4/G4ShortLivedConstructor.hh"
#include "Geant4/G4ParticleTable.hh"
#if G4VERSION<4099
#include "Geant4/G4ShortLivedTable.hh"
#endif
#include "Geant4/G4PhaseSpaceDecayChannel.hh"
#include "Geant4/G4VDecayChannel.hh"
#include "Geant4/G4DecayTable.hh"

#include "Geant4/G4DiQuarks.hh"
#include "Geant4/G4Quarks.hh"
#include "Geant4/G4Gluons.hh"

//------------------------debug switches
//#define DEBUG_LightFragmentationTest 1


//********************************************************************************
// Constructors

PBARVLongitudinalStringDecay::PBARVLongitudinalStringDecay()
{
   MassCut  = 0.35*GeV; 
   ClusterMass = 0.15*GeV;

   SmoothParam      = 0.9; 
   StringLoopInterrupt    = 1000;
   ClusterLoopInterrupt   =  500;

// Changable Parameters below.
   SigmaQT = 0.5 * GeV;  // 0.5 0.1
   
   StrangeSuppress  = 0.44;    //  27 % strange quarks produced, ie. u:d:s=1:1:0.27
   DiquarkSuppress  = 0.07;
   DiquarkBreakProb = 0.1;
   
   //... pspin_meson is probability to create vector meson 
   pspin_meson = 0.5;

   //... pspin_barion is probability to create 3/2 barion 
   pspin_barion = 0.5;

   //... vectorMesonMix[] is quark mixing parameters for vector mesons (Variable spin = 3)
   vectorMesonMix.resize(6);
   vectorMesonMix[0] = 0.5;
   vectorMesonMix[1] = 0.0;
   vectorMesonMix[2] = 0.5;
   vectorMesonMix[3] = 0.0;
   vectorMesonMix[4] = 1.0;
   vectorMesonMix[5] = 1.0; 

   //... scalarMesonMix[] is quark mixing parameters for scalar mesons (Variable spin=1)
   scalarMesonMix.resize(6);
   scalarMesonMix[0] = 0.5; 
   scalarMesonMix[1] = 0.25; 
   scalarMesonMix[2] = 0.5; 
   scalarMesonMix[3] = 0.25; 
   scalarMesonMix[4] = 1.0; 
   scalarMesonMix[5] = 0.5; 

// Parameters may be changed until the first fragmentation starts
   PastInitPhase=false;
   hadronizer = new G4HadronBuilder(pspin_meson,pspin_barion,
		   		scalarMesonMix,vectorMesonMix);
   Kappa = 1.0 * GeV/fermi;


}
   

PBARVLongitudinalStringDecay::~PBARVLongitudinalStringDecay()
   {
   delete hadronizer;
   }

//=============================================================================

// Operators

//const  & PBARVLongitudinalStringDecay::operator=(const PBARVLongitudinalStringDecay &)
//    {
//    }

//-----------------------------------------------------------------------------

int PBARVLongitudinalStringDecay::operator==(const PBARVLongitudinalStringDecay &) const
    {
	throw G4HadronicException(__FILE__, __LINE__, "PBARVLongitudinalStringDecay::operator== forbidden");
	return false;
    }

//-------------------------------------------------------------------------------------

int PBARVLongitudinalStringDecay::operator!=(const PBARVLongitudinalStringDecay &) const
    {
	throw G4HadronicException(__FILE__, __LINE__, "PBARVLongitudinalStringDecay::operator!= forbidden");
	return true;
    }

//***********************************************************************************

// For changing Mass Cut used for selection of very small mass strings
void PBARVLongitudinalStringDecay::SetMassCut(G4double aValue){MassCut=aValue;}

//-----------------------------------------------------------------------------

// For handling a string with very low mass

G4KineticTrackVector* PBARVLongitudinalStringDecay::LightFragmentationTest(const
		G4ExcitedString * const string)
{
   // Check string decay threshold
		
	G4KineticTrackVector * result=0;  // return 0 when string exceeds the mass cut
	
	pDefPair hadrons((G4ParticleDefinition *)0,(G4ParticleDefinition *)0);

	G4FragmentingString aString(*string);

	if ( sqr(FragmentationMass(&aString,0,&hadrons)+MassCut) < aString.Mass2()) {
		return 0;
	}

// The string mass is very low ---------------------------
	
	result=new G4KineticTrackVector;
        
	if ( hadrons.second ==0 )
	{
// Substitute string by light hadron, Note that Energy is not conserved here!

/*		 
#ifdef DEBUG_LightFragmentationTest
	       G4cout << "VlongSF Warning replacing string by single hadron " <<G4endl;
	       G4cout << hadrons.first->GetParticleName() 
		      << "string .. " << string->Get4Momentum() << " " 
		      << string->Get4Momentum().m() << G4endl;
#endif		      
*/
	       G4ThreeVector   Mom3 = string->Get4Momentum().vect();
	       G4LorentzVector Mom(Mom3, 
	       			   std::sqrt(Mom3.mag2() + 
                                             sqr(hadrons.first->GetPDGMass())));
               result->push_back(new G4KineticTrack(hadrons.first, 0, 
                                                  string->GetPosition(),
                                                          Mom));
	} else 
	{
//... string was qq--qqbar type: Build two stable hadrons,

#ifdef DEBUG_LightFragmentationTest
	       G4cout << "VlongSF Warning replacing qq-qqbar string by TWO hadrons " 
		      << hadrons.first->GetParticleName() << " / " 
		      << hadrons.second->GetParticleName()
		      << "string .. " << string->Get4Momentum() << " " 
		      << string->Get4Momentum().m() << G4endl;
#endif		      

	       G4LorentzVector  Mom1, Mom2;
	       Sample4Momentum(&Mom1, hadrons.first->GetPDGMass(), 
			       &Mom2,hadrons.second->GetPDGMass(),
			       string->Get4Momentum().mag());

	       result->push_back(new G4KineticTrack(hadrons.first, 0, 
                                                    string->GetPosition(), 
                                                            Mom1));
	       result->push_back(new G4KineticTrack(hadrons.second, 0, 
                                                    string->GetPosition(), 
                                                    Mom2));

               G4ThreeVector Velocity = string->Get4Momentum().boostVector();
               result->Boost(Velocity);          
	}

	return result;
	
}

//----------------------------------------------------------------------------------------

G4double PBARVLongitudinalStringDecay::FragmentationMass(
            const G4FragmentingString *	const string,
		Pcreate build, pDefPair * pdefs       )
{
	
        G4double mass;
        static G4bool NeedInit(true);
	static std::vector<double> nomix;
	static G4HadronBuilder * minMassHadronizer;
	if ( NeedInit ) 
	{
	   NeedInit = false;
	   nomix.resize(6);
	   for ( G4int i=0; i<6 ; i++ ) nomix[i]=0;

//	   minMassHadronizer=new G4HadronBuilder(pspin_meson,pspin_barion,nomix,nomix);
	   minMassHadronizer=hadronizer;
	}

	if ( build==0 ) build=&G4HadronBuilder::BuildLowSpin;

        G4ParticleDefinition *Hadron1, *Hadron2=0;

        if (!string->FourQuarkString() )
        {
           // spin 0 meson or spin 1/2 barion will be built

//G4cout<<"String Left Right "<<string->GetLeftParton()<<" "<<string->GetRightParton()<<G4endl;
           Hadron1 = (minMassHadronizer->*build)(string->GetLeftParton(),
			                         string->GetRightParton());
//G4cout<<"Hadron1 "<<Hadron1->GetParticleName()<<G4endl;
           mass= (Hadron1)->GetPDGMass();
        } else
        {
           //... string is qq--qqbar: Build two stable hadrons,
           //...    with extra uubar or ddbar quark pair
	   G4int iflc = (G4UniformRand() < 0.5)? 1 : 2;
	   if (string->GetLeftParton()->GetPDGEncoding() < 0) iflc = -iflc;

	   //... theSpin = 4; spin 3/2 baryons will be built
	   Hadron1 = (minMassHadronizer->*build)(string->GetLeftParton(),
                                                 FindParticle(iflc)       );
	   Hadron2 = (minMassHadronizer->*build)(string->GetRightParton(),
                                                 FindParticle(-iflc)      );
           mass = (Hadron1)->GetPDGMass() + (Hadron2)->GetPDGMass();
        }
	
	if ( pdefs != 0 ) 
	{ // need to return hadrons as well....
	   pdefs->first  = Hadron1;
	   pdefs->second = Hadron2;
	}
	   
        return mass;
}

//----------------------------------------------------------------------------

G4ParticleDefinition* PBARVLongitudinalStringDecay::FindParticle(G4int Encoding) 
   {
   G4ParticleDefinition* ptr = G4ParticleTable::GetParticleTable()->FindParticle(Encoding);
      if (ptr == NULL)
       {
       G4cout << "Particle with encoding "<<Encoding<<" does not exist!!!"<<G4endl;
       throw G4HadronicException(__FILE__, __LINE__, "Check your particle table");
       }
   return ptr;    
   }

//-----------------------------------------------------------------------------
//   virtual void Sample4Momentum(G4LorentzVector* Mom,     G4double Mass, 
//                                G4LorentzVector* AntiMom, G4double AntiMass, 
//                                G4double InitialMass)=0; 
//-----------------------------------------------------------------------------

//*********************************************************************************
//   For decision on continue or stop string fragmentation
//   virtual G4bool StopFragmenting(const G4FragmentingString  * const string)=0;
//   virtual G4bool IsFragmentable(const G4FragmentingString * const string)=0;

//   If a string can not fragment, make last break into 2 hadrons
//   virtual G4bool SplitLast(G4FragmentingString * string, 
//                            G4KineticTrackVector * LeftVector,
//                            G4KineticTrackVector * RightVector)=0;
//-----------------------------------------------------------------------------
//
//   If a string fragments, do the following
//
//   For transver of a string to its CMS frame
//-----------------------------------------------------------------------------

G4ExcitedString *PBARVLongitudinalStringDecay::CPExcited(const G4ExcitedString & in)
{
	G4Parton *Left=new G4Parton(*in.GetLeftParton());
	G4Parton *Right=new G4Parton(*in.GetRightParton());
	return new G4ExcitedString(Left,Right,in.GetDirection());
}

//-----------------------------------------------------------------------------

G4KineticTrack * PBARVLongitudinalStringDecay::Splitup(
		        G4FragmentingString *string, 
			G4FragmentingString *&newString)
{
//G4cout<<"Start SplitUP"<<G4endl;
       //... random choice of string end to use for creating the hadron (decay)   
       G4int SideOfDecay = (G4UniformRand() < 0.5)? 1: -1;
       if (SideOfDecay < 0)
       {
	  string->SetLeftPartonStable();
       } else
       {
          string->SetRightPartonStable();
       }

       G4ParticleDefinition *newStringEnd;
       G4ParticleDefinition * HadronDefinition;
       if (string->DecayIsQuark())
       {
       	   HadronDefinition= QuarkSplitup(string->GetDecayParton(), newStringEnd);
/*yzy
G4cout << "QuarkSplitup" << G4endl;
*/
       } else {
           HadronDefinition= DiQuarkSplitup(string->GetDecayParton(), newStringEnd);
/*yzy
G4cout << "DiQuarkSplitup" << G4endl;
*/
       }      
//*yzy
//G4cout<<"New had "<<HadronDefinition->GetParticleName()<<G4endl;
// create new String from old, ie. keep Left and Right order, but replace decay

       newString=new G4FragmentingString(*string,newStringEnd); // To store possible
                                                                // quark containt of new string
//G4cout<<"SplitEandP "<<G4endl;
       G4LorentzVector* HadronMomentum=SplitEandP(HadronDefinition, string, newString);

       delete newString; newString=0;                          
	
       G4KineticTrack * Hadron =0;
       if ( HadronMomentum != 0 ) {    

	   G4ThreeVector   Pos;
	   Hadron = new G4KineticTrack(HadronDefinition, 0,Pos, *HadronMomentum);
 
	   newString=new G4FragmentingString(*string,newStringEnd,
	   				HadronMomentum);

	   delete HadronMomentum;
       }      
//G4cout<<"End SplitUP"<<G4endl;
       return Hadron;
}

//--------------------------------------------------------------------------------------

G4ParticleDefinition *
		PBARVLongitudinalStringDecay::QuarkSplitup(G4ParticleDefinition*
		decay, G4ParticleDefinition *&created)
{
    G4int IsParticle=(decay->GetPDGEncoding()>0) ? -1 : +1; // if we have a quark, 
                                                            // we need antiquark 
                                                            // (or diquark)
    pDefPair QuarkPair = CreatePartonPair(IsParticle);
    created = QuarkPair.second;
    return hadronizer->Build(QuarkPair.first, decay);
    
}

//-----------------------------------------------------------------------------

G4ParticleDefinition *PBARVLongitudinalStringDecay::DiQuarkSplitup(
							G4ParticleDefinition* decay,
							G4ParticleDefinition *&created)
{
   //... can Diquark break or not? 
   if (G4UniformRand() < DiquarkBreakProb ){
   //... Diquark break

      G4int stableQuarkEncoding = decay->GetPDGEncoding()/1000;
      G4int decayQuarkEncoding = (decay->GetPDGEncoding()/100)%10;
      if (G4UniformRand() < 0.5)
         {
         G4int Swap = stableQuarkEncoding;
         stableQuarkEncoding = decayQuarkEncoding;
         decayQuarkEncoding = Swap;
         }

      G4int IsParticle=(decayQuarkEncoding>0) ? -1 : +1; 
			// if we have a quark, we need antiquark)
      pDefPair QuarkPair = CreatePartonPair(IsParticle,false);  // no diquarks wanted
      //... Build new Diquark
      G4int QuarkEncoding=QuarkPair.second->GetPDGEncoding();
      G4int i10  = std::max(std::abs(QuarkEncoding), std::abs(stableQuarkEncoding));
      G4int i20  = std::min(std::abs(QuarkEncoding), std::abs(stableQuarkEncoding));
      G4int spin = (i10 != i20 && G4UniformRand() <= 0.5)? 1 : 3;
      G4int NewDecayEncoding = -1*IsParticle*(i10 * 1000 + i20 * 100 + spin);
      created = FindParticle(NewDecayEncoding);
      G4ParticleDefinition * decayQuark=FindParticle(decayQuarkEncoding);
      G4ParticleDefinition * had=hadronizer->Build(QuarkPair.first, decayQuark);
      return had;
//      return hadronizer->Build(QuarkPair.first, decayQuark);
   
   } else {
   //... Diquark does not break
 
      G4int IsParticle=(decay->GetPDGEncoding()>0) ? +1 : -1; 
			// if we have a diquark, we need quark)
      pDefPair QuarkPair = CreatePartonPair(IsParticle,false);  // no diquarks wanted
      created = QuarkPair.second;

      G4ParticleDefinition * had=hadronizer->Build(QuarkPair.first, decay);
      return had;
//      return G4ParticleDefinition * had=hadronizer->Build(QuarkPair.first, decay);
   }
}

//-----------------------------------------------------------------------------

G4int PBARVLongitudinalStringDecay::SampleQuarkFlavor(void)
   {
   return (1 + (int)(G4UniformRand()/StrangeSuppress));
   }

//-----------------------------------------------------------------------------

PBARVLongitudinalStringDecay::pDefPair PBARVLongitudinalStringDecay::CreatePartonPair(G4int NeedParticle,G4bool AllowDiquarks)
{
//  NeedParticle = +1 for Particle, -1 for Antiparticle

    if ( AllowDiquarks && G4UniformRand() < DiquarkSuppress )
    {
      // Create a Diquark - AntiDiquark pair , first in pair is anti to IsParticle
      G4int q1  = SampleQuarkFlavor();
      G4int q2  = SampleQuarkFlavor();
      G4int spin = (q1 != q2 && G4UniformRand() <= 0.5)? 1 : 3;
                                     //   convention: quark with higher PDG number is first
      G4int PDGcode = (std::max(q1,q2) * 1000 + std::min(q1,q2) * 100 + spin) * NeedParticle;
      return pDefPair (FindParticle(-PDGcode),FindParticle(PDGcode));
      

    } else {
      // Create a Quark - AntiQuark pair, first in pair  IsParticle
      G4int PDGcode=SampleQuarkFlavor()*NeedParticle;
      return pDefPair (FindParticle(PDGcode),FindParticle(-PDGcode));
    }

}

//-----------------------------------------------------------------------------
G4ThreeVector PBARVLongitudinalStringDecay::SampleQuarkPt(G4double ptMax)
   {
   G4double Pt;
   if ( ptMax < 0 ) {
      // sample full gaussian
      Pt = -std::log(G4UniformRand());
   } else {
      // sample in limited range
      Pt = -std::log(CLHEP::RandFlat::shoot(std::exp(-sqr(ptMax)/sqr(SigmaQT)), 1.));
   }
   Pt = SigmaQT * std::sqrt(Pt);
   G4double phi = 2.*pi*G4UniformRand();
   return G4ThreeVector(Pt * std::cos(phi),Pt * std::sin(phi),0);
   }

//******************************************************************************

void PBARVLongitudinalStringDecay::CalculateHadronTimePosition(G4double theInitialStringMass, G4KineticTrackVector* Hadrons)
   {

//   `yo-yo` formation time
//   const G4double kappa = 1.0 * GeV/fermi/4.;      
   G4double kappa = GetStringTensionParameter();
   for(size_t c1 = 0; c1 < Hadrons->size(); c1++)
      {
      G4double SumPz = 0; 
      G4double SumE  = 0;
      for(size_t c2 = 0; c2 < c1; c2++)
         {
         SumPz += Hadrons->operator[](c2)->Get4Momentum().pz();
         SumE  += Hadrons->operator[](c2)->Get4Momentum().e();   
         } 
      G4double HadronE  = Hadrons->operator[](c1)->Get4Momentum().e();
      G4double HadronPz = Hadrons->operator[](c1)->Get4Momentum().pz();
      Hadrons->operator[](c1)->SetFormationTime(
(theInitialStringMass - 2.*SumPz + HadronE - HadronPz)/(2.*kappa)/c_light); 

      G4ThreeVector aPosition(0, 0,     
(theInitialStringMass - 2.*SumE  - HadronE + HadronPz)/(2.*kappa));
      Hadrons->operator[](c1)->SetPosition(aPosition);

      } 
   }

//-----------------------------------------------------------------------------

void PBARVLongitudinalStringDecay::SetSigmaTransverseMomentum(G4double aValue)
{
	if ( PastInitPhase ) {
		throw G4HadronicException(__FILE__, __LINE__, "4VLongitudinalStringDecay::SetSigmaTransverseMomentum after FragmentString() not allowed");
	} else {
		SigmaQT = aValue;
	}
}

//----------------------------------------------------------------------------------------------------------

void PBARVLongitudinalStringDecay::SetStrangenessSuppression(G4double aValue)
{
	if ( PastInitPhase ) {
		throw G4HadronicException(__FILE__, __LINE__, "4VLongitudinalStringDecay::SetStrangenessSuppression after FragmentString() not allowed");
	} else {
		StrangeSuppress = aValue;
	}
}

//----------------------------------------------------------------------------------------------------------

void PBARVLongitudinalStringDecay::SetDiquarkSuppression(G4double aValue)
{
	if ( PastInitPhase ) {
		throw G4HadronicException(__FILE__, __LINE__, "4VLongitudinalStringDecay::SetDiquarkSuppression after FragmentString() not allowed");
	} else {
		DiquarkSuppress = aValue;
	}
}

//----------------------------------------------------------------------------------------

void PBARVLongitudinalStringDecay::SetDiquarkBreakProbability(G4double aValue)
{
	if ( PastInitPhase ) {
		throw G4HadronicException(__FILE__, __LINE__, "4VLongitudinalStringDecay::SetDiquarkBreakProbability after FragmentString() not allowed");
	} else {
		DiquarkBreakProb = aValue;
	}
}

//----------------------------------------------------------------------------------------------------------

void PBARVLongitudinalStringDecay::SetVectorMesonProbability(G4double aValue)
{
	if ( PastInitPhase ) {
		throw G4HadronicException(__FILE__, __LINE__, "PBARVLongitudinalStringDecay::SetVectorMesonProbability after FragmentString() not allowed");
	} else {
		pspin_meson = aValue;
		delete hadronizer;
		hadronizer = new G4HadronBuilder(pspin_meson,pspin_barion,
		   		scalarMesonMix,vectorMesonMix);
	}
}

//----------------------------------------------------------------------------------------------------------

void PBARVLongitudinalStringDecay::SetSpinThreeHalfBarionProbability(G4double aValue)
{
	if ( PastInitPhase ) {
		throw G4HadronicException(__FILE__, __LINE__, "PBARVLongitudinalStringDecay::SetSpinThreeHalfBarionProbability after FragmentString() not allowed");
	} else {
		pspin_barion = aValue;
		delete hadronizer;
		hadronizer = new G4HadronBuilder(pspin_meson,pspin_barion,
		   		scalarMesonMix,vectorMesonMix);
	}
}

//----------------------------------------------------------------------------------------------------------

void PBARVLongitudinalStringDecay::SetScalarMesonMixings(std::vector<G4double> aVector)
{
	if ( PastInitPhase ) {
		throw G4HadronicException(__FILE__, __LINE__, "PBARVLongitudinalStringDecay::SetScalarMesonMixings after FragmentString() not allowed");
	} else {
	  if ( aVector.size() < 6 ) 
	      throw G4HadronicException(__FILE__, __LINE__, "PBARVLongitudinalStringDecay::SetScalarMesonMixings( argument Vector too small");
	  scalarMesonMix[0] = aVector[0];
	  scalarMesonMix[1] = aVector[1];
	  scalarMesonMix[2] = aVector[2];
	  scalarMesonMix[3] = aVector[3];
	  scalarMesonMix[4] = aVector[4];
	  scalarMesonMix[5] = aVector[5];
	  delete hadronizer;
	  hadronizer = new G4HadronBuilder(pspin_meson,pspin_barion,
		   		scalarMesonMix,vectorMesonMix);
	}
}

//----------------------------------------------------------------------------------------------------------

void PBARVLongitudinalStringDecay::SetVectorMesonMixings(std::vector<G4double> aVector)
{
	if ( PastInitPhase ) {
		throw G4HadronicException(__FILE__, __LINE__, "PBARVLongitudinalStringDecay::SetVectorMesonMixings after FragmentString() not allowed");
	} else {
	  if ( aVector.size() < 6 ) 
	      throw G4HadronicException(__FILE__, __LINE__, "PBARVLongitudinalStringDecay::SetVectorMesonMixings( argument Vector too small");
	  vectorMesonMix[0] = aVector[0];
	  vectorMesonMix[1] = aVector[1];
	  vectorMesonMix[2] = aVector[2];
	  vectorMesonMix[3] = aVector[3];
	  vectorMesonMix[4] = aVector[4];
	  vectorMesonMix[5] = aVector[5];
	  delete hadronizer;
	  hadronizer = new G4HadronBuilder(pspin_meson,pspin_barion,
		   		scalarMesonMix,vectorMesonMix);
  
	}
}

//-------------------------------------------------------------------------------------------
void PBARVLongitudinalStringDecay::SetStringTensionParameter(G4double aValue)// Uzhi 20 June 08
{
          Kappa = aValue * GeV/fermi;
}	
//**************************************************************************************

#endif
