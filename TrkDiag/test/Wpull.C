#define Wpull_cxx
#include "Wpull.h"

void Wpull::Loop()
{
//   In a ROOT session, you can do:
//      root> .L Wpull.C
//      root> Wpull t
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

     double maxwp(4.0);

     if (dem__status > 0 && demmc_proc == 56 && dem__mom > 100.0) {
       double wpavg;
       unsigned nactive(0), nwpactive(0);
       for (int i = 0; i < demtshmc_ ; ++i){
	 if(demtsh__active[i] && demtshmc__rel[i] == 0){
	   ++nactive;
	   double wp = fabs(demtsh__wdist[i]-demtsh__hlen[i])/demtsh__werr[i];
	   wpavg += wp;
	   if(wp > maxwp)++nwpactive;
	 }
       }
       wpavg /= nactive;
       double wpfrac = nwpactive/(float)nactive;
       double dm = dem__mom - demmcent_mom;
       _wpavg->Fill(wpavg);
       _wpfrac->Fill(wpfrac);
       _wpavgmom->Fill(dm,wpavg);
       _wpfracmom->Fill(dm,wpfrac);
       _wpavgtq->Fill(dem__trkqual,wpavg);
       _wpfractq->Fill(dem__trkqual,wpfrac);
     }
   }

}
