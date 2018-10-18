// Plot the clustering parameters (see TrkReco/src/TNTCluster for definition) for clusters 
// - produced by a single particle and entirely contained in one cluster (extact)
// - produced by a single particle and split over several clusters (split)
// - produced by a several particles (mult)
//
// In you analyzer, you must you must run
//     BD: @local::BD
//    TCD : @local::TCD
//and set 
// physics.producers.FlagBkgHits.SaveBkgClusters :true
// physics.analyzers.TCD.diagLevel :3
//


#include <TFile.h>
#include <TH1.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <vector>
#include <iostream>


void cluAnaParameter()
{
   
   //   int entry = 0;
   
   gStyle->SetOptStat(0);
   
   TH1F *_valT[3],*_valR[3],*_valP[3],*_valTR[3],*_valR2[3];
   for (int i=0;i<3;++i)_valT[i]  = new TH1F(Form("valT%i",i),"valT",100,1,51);
   for (int i=0;i<3;++i)_valR[i]  = new TH1F(Form("valR%i",i),"valR",100,1,101);
   for (int i=0;i<3;++i)_valP[i]  = new TH1F(Form("valP%i",i),"valP",100,0.01,1);
   for (int i=0;i<3;++i)_valTR[i] = new TH1F(Form("valTR%i",i),"valTR",100,1,21);
   for (int i=0;i<3;++i)_valR2[i] = new TH1F(Form("valR2%i",i),"valR2",100,1,201);

   for (int i=0;i<3;++i){   
     _valT[i]->SetLineColor(i+1);
     _valTR[i]->SetLineColor(i+1);
     _valR[i]->SetLineColor(i+1);
     _valR2[i]->SetLineColor(i+1);
     _valP[i]->SetLineColor(i+1);
   }
   
   TTreeReader fReader;  
   TTreeReaderValue<Int_t> iev      = {fReader, "iev"};
   TTreeReaderValue<Int_t> cluIdx   = {fReader, "cluIdx"};
   TTreeReaderValue<Int_t> nchits   = {fReader, "nchits"};
   TTreeReaderValue<Int_t> nshits   = {fReader, "nshits"};
   TTreeReaderValue<Float_t> mvaout = {fReader, "mvaout"};
   TTreeReaderValue<Float_t> pmom   = {fReader, "pmom"};
   TTreeReaderValue<Int_t> ppid     = {fReader, "ppid"};
   TTreeReaderValue<Int_t> ppdg     = {fReader, "ppdg"};
   TTreeReaderValue<Int_t> pgen     = {fReader, "pgen"};
   TTreeReaderValue<Int_t> pproc    = {fReader, "pproc"};
   TTreeReaderValue<Int_t> nprimary = {fReader, "nprimary"};
   TTreeReaderValue<Int_t> nconv    = {fReader, "nconv"};
   TTreeReaderValue<Int_t> ndelta   = {fReader, "ndelta"};
   TTreeReaderValue<Int_t> ncompt   = {fReader, "ncompt"};
   TTreeReaderValue<Int_t> ngconv   = {fReader, "ngconv"};
   TTreeReaderValue<Int_t> nebkg    = {fReader, "nebkg"};
   TTreeReaderValue<Int_t> nprot    = {fReader, "nprot"};
   TTreeReaderValue<Int_t> ncontrib = {fReader, "ncontrib"};
   TTreeReaderArray<Int_t> icontrib = {fReader, "icontrib"};
   TTreeReaderValue<Int_t> nindex   = {fReader, "nindex"};
   TTreeReaderArray<Int_t> hitidx   = {fReader, "hitidx"};

   TTreeReader fReader2;  
   TTreeReaderValue<Int_t>   iev2      = {fReader2, "iev"};
   TTreeReaderValue<Int_t>   nhits     = {fReader2, "nhits"};
   TTreeReaderArray<Int_t>   hitNcombo = {fReader2, "hitNcombo"};
   TTreeReaderArray<Float_t> hitRad    = {fReader2, "hitRad"};
   TTreeReaderArray<Float_t> hitPhi    = {fReader2, "hitPhi"};
   TTreeReaderArray<Float_t> hitTime   = {fReader2, "hitTime"};
   TTreeReaderArray<Int_t>   hitPdg    = {fReader2, "hitPdg"};
   TTreeReaderArray<Int_t>   hitCrCode = {fReader2, "hitCrCode"};
   TTreeReaderArray<Int_t>   hitGen    = {fReader2, "hitGen"};

   TTreeReader fReader3;  
   TTreeReaderValue<Int_t>   ncluIter  = {fReader3, "ncluIter"};
   TTreeReaderArray<Int_t>   cluId     = {fReader3, "cluId"};
   TTreeReaderArray<Int_t>   cluNpass  = {fReader3, "cluNpass"};
   TTreeReaderArray<Int_t>   cnhi      = {fReader3, "cnhi"};
   TTreeReaderArray<Float_t> cluRad    = {fReader3, "cluRad"};
   TTreeReaderArray<Float_t> cluPhi    = {fReader3, "cluPhi"};
   TTreeReaderArray<Float_t> cluTime   = {fReader3, "cluTime"};
   TTreeReaderArray<Float_t> cluR2diff = {fReader3, "cluR2diff"};
   TTreeReaderArray<Float_t> cluRdiff  = {fReader3, "cluRdiff"};
   TTreeReaderArray<Float_t> cluPdiff  = {fReader3, "cluPdiff"};
   TTreeReaderArray<Float_t> cluTdiff  = {fReader3, "cluTdiff"};
   TTreeReaderArray<Float_t> cluTRdiff = {fReader3, "cluTRdiff"};   
   TTreeReaderValue<Int_t>   nhitClu   = {fReader3, "nhitClu"};
   TTreeReaderArray<Int_t>   hcIdxClu  = {fReader3, "hcIdxClu"};
   TTreeReaderArray<Int_t>   hcIdxHit  = {fReader3, "hcIdxHit"};
   TTreeReaderArray<Int_t>   hcNpass   = {fReader3, "hcNpass"};   
   TTreeReaderValue<Int_t>   niter     = {fReader3, "niter"};
   TTreeReaderArray<Int_t>   nclu      = {fReader3, "nclu"};
   TTreeReaderArray<Int_t>   nChanged  = {fReader3, "nChanged"};
   TTreeReaderArray<Float_t> odist     = {fReader3, "odist"};
   TTreeReaderArray<Float_t> tdist     = {fReader3, "tdist"};
   
 
   TFile f("trkDiag.root");
   
   TTree* tree = (TTree*) f.Get("BD/bkgdiag");   
   TTree* tree2 = (TTree*) f.Get("BD/bkgdiag2");   
   TTree* tree3 = (TTree*) f.Get("FlagBkgHits/idiag");   
   fReader.SetTree(tree);      
   fReader2.SetTree(tree2);      
   fReader3.SetTree(tree3);      
   
   
   
   for (int entry = 1; entry < 30;++entry){
   
       //set fReader2 to correct entry
       int in2(0);
       for (in2=0;in2<tree2->GetEntriesFast();++in2){
           fReader2.SetEntry(in2);       
           if (*iev2 == entry) break;
       }

       if (*iev2 != entry) continue;

       //synchronize cluster tree and fill map
       fReader3.SetEntry(in2);       

       map<int,float> contmapT,contmapTR,contmapR,contmapP,contmapR2;     
       for (int ic=0;ic<*ncluIter;++ic){
          contmapT[cluId[ic]] = cluTdiff[ic];
          contmapR[cluId[ic]] = cluRdiff[ic];
          contmapP[cluId[ic]] = cluPdiff[ic];
          contmapR2[cluId[ic]] = cluR2diff[ic];
          contmapTR[cluId[ic]] = cluTRdiff[ic];
       }


       //now fill the contribution list
       map<int,vector<int>> contmap;
       for (int in=0;in<tree->GetEntriesFast();++in){
          fReader.SetEntry(in);
          if (*iev != entry) continue;
          for (size_t ib=0;ib <*ncontrib;++ib) contmap[icontrib[ib]].push_back(in);
       }     
   
       for (int in=0;in<tree->GetEntriesFast();++in){
         fReader.SetEntry(in);
         if (*iev != entry) continue;
         
         //at least 3 hits or this won't say much
         if (*nindex<3) continue; 
         //if (*mvaout>0.5) continue; 
         
         if (*ncontrib==1 && contmap[icontrib[0]].size()==1) _valT[0]->Fill(contmapT[*cluIdx]);        
         if (*ncontrib==1 && contmap[icontrib[0]].size()>1)  _valT[1]->Fill(contmapT[*cluIdx]);     
         if (*ncontrib>1)                                    _valT[2]->Fill(contmapT[*cluIdx]);
         if (*ncontrib==1 && contmap[icontrib[0]].size()==1) _valTR[0]->Fill(contmapTR[*cluIdx]);        
         if (*ncontrib==1 && contmap[icontrib[0]].size()>1)  _valTR[1]->Fill(contmapTR[*cluIdx]);     
         if (*ncontrib>1)                                    _valTR[2]->Fill(contmapTR[*cluIdx]);
         if (*ncontrib==1 && contmap[icontrib[0]].size()==1) _valR[0]->Fill(contmapR[*cluIdx]);        
         if (*ncontrib==1 && contmap[icontrib[0]].size()>1)  _valR[1]->Fill(contmapR[*cluIdx]);     
         if (*ncontrib>1)                                    _valR[2]->Fill(contmapR[*cluIdx]);
         if (*ncontrib==1 && contmap[icontrib[0]].size()==1) _valP[0]->Fill(contmapP[*cluIdx]);        
         if (*ncontrib==1 && contmap[icontrib[0]].size()>1)  _valP[1]->Fill(contmapP[*cluIdx]);     
         if (*ncontrib>1)                                    _valP[2]->Fill(contmapP[*cluIdx]);
         if (*ncontrib==1 && contmap[icontrib[0]].size()==1) _valR2[0]->Fill(contmapR2[*cluIdx]);        
         if (*ncontrib==1 && contmap[icontrib[0]].size()>1)  _valR2[1]->Fill(contmapR2[*cluIdx]);     
         if (*ncontrib>1)                                    _valR2[2]->Fill(contmapR2[*cluIdx]);
       }
    }
    
    
    


   TLegend leg(0.7,0.7,0.9,0.9);
   leg.AddEntry(_valT[0],"exact","l");
   leg.AddEntry(_valT[1],"split","l");
   leg.AddEntry(_valT[2],"mult","l");
  
   
   TCanvas *c1 = new TCanvas("c1","c1",1000,1000);   
   c1->Divide(2,2);
   
   c1->cd(1);   
   _valT[1]->Draw();
   _valT[0]->Draw("same");
   _valT[2]->Draw("same");
   leg.Draw();
   c1->cd(2);   
   _valR[1]->Draw();
   _valR[0]->Draw("same");
   _valR[2]->Draw("same");
   c1->cd(3);   
   _valP[1]->Draw();
   _valP[0]->Draw("same");
   _valP[2]->Draw("same");
   c1->SaveAs("cluParam1.pdf");
   
   TCanvas *c2 = new TCanvas("c2","c2",1000,1000);   
   c2->Divide(2,2);
   c2->cd(1);   
   _valTR[1]->Draw();
   _valTR[0]->Draw("same");
   _valTR[2]->Draw("same");
   leg.Draw();
   c2->cd(2);   
   _valR2[1]->Draw();
   _valR2[0]->Draw("same");
   _valR2[2]->Draw("same");
   c2->SaveAs("cluParam2.pdf");

   f.Close();
}
