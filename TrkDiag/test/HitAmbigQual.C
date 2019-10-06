#define HitAmbigQual_cxx
#include "HitAmbigQual.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void HitAmbigQual::Loop()
{
//   In a ROOT session, you can do:
//      root> .L HitAmbigQual.C
//      root> HitAmbigQual t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     if(dem__status > 0){
       // loop over hits, find ones in the smame plane
       //
       Float_t spaspanel(0.0), spasplane(0.0), spasstation(0.0);
       Float_t spasdpanel(0.0), spasdplane(0.0), spasdstation(0.0);
       Float_t spapanel(0.0), spaplane(0.0), spastation(0.0);
       Float_t spadpanel(0.0), spadplane(0.0), spadstation(0.0);

      Float_t dm = dem__mom - demmcent_mom;

       for(Int_t ihit = 0;ihit < demtsh_; ++ihit){
	 if(demtsh__active[ihit] && demtsh__ambig[ihit]!= 0 ){
	   for(Int_t jhit = ihit+1;jhit < demtsh_; ++jhit){
	     if(demtsh__active[jhit] && demtsh__ambig[jhit]!= 0 &&
		 demtsh__plane[jhit]/2 ==  demtsh__plane[ihit]/2){
	       float df = demtsh__trklen[jhit]-demtsh__trklen[ihit];
	       int da = demtsh__ambig[jhit]*demtsh__ambig[ihit]; 
	       float pd = demtsh__rdrift[jhit]*demtsh__rdrift[ihit]*da;
	       _pastation->Fill(pd);
	       _pdstation->Fill(df);
	       spastation += abs(pd);
	       spadstation += abs(pd)/df;
	       if(da < 0){
		 spasstation += pd;
		 spasdstation += pd/df;
	       }
	       if(demtsh__plane[jhit] == demtsh__plane[ihit] ){
		 _paplane->Fill(pd);
		 _pdplane->Fill(df);
		 spaplane += abs(pd);
		 spadplane += abs(pd)/df;
		 if(da < 0){
		   spasplane += pd;
		   spasdplane += pd/df;
		 }
		 if(demtsh__panel[jhit] == demtsh__panel[ihit] ){
		   _papanel->Fill(pd);
		   _pdpanel->Fill(df);
		   spapanel += abs(pd);
		   spadpanel += abs(pd)/df;
		   if(da < 0){
		     spaspanel += pd;
		     spasdpanel += pd/df;
		   }
		 }
	       }
	     }
	   }
	 }
       }
       _spaplane->Fill(spaplane);
       _spapanel->Fill(spapanel);
       _spastation->Fill(spastation);
       _spadplane->Fill(spadplane);
       _spadpanel->Fill(spadpanel);
       _spadstation->Fill(spadstation);
       _spasplane->Fill(-spasplane);
       _spaspanel->Fill(-spaspanel);
       _spasstation->Fill(-spasstation);
       _spasdplane->Fill(-spasdplane);
       _spasdpanel->Fill(-spasdpanel);
       _spasdstation->Fill(-spasdstation);

       _spasfplane->Fill(-spasplane/spaplane);
       _spasfpanel->Fill(-spaspanel/spapanel);
       _spasfstation->Fill(-spasstation/spastation);
       _spasdfplane->Fill(-spasdplane/spadplane);
       _spasdfpanel->Fill(-spasdpanel/spadpanel);
       _spasdfstation->Fill(-spasdstation/spadstation);
           
       _dmplane->Fill(-spasdplane/spadplane,dm);
       _dmpanel->Fill(-spasdpanel/spadpanel,dm);
       _dmstation->Fill(-spasdstation/spadstation,dm);

       if(dem__trkqual > 0.4){
	 _dmtplane->Fill(-spasdplane/spadplane,dm);
	 _dmtpanel->Fill(-spasdpanel/spadpanel,dm);
	 _dmtstation->Fill(-spasdstation/spadstation,dm);
       }


     }
   }
}
