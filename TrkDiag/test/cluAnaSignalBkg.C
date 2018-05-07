// Plot the r,phi,time coordinates of clusters with signal CE in them (black = bkg, red = CE). The bottom right plane 
// plots the rescontructed clusters with signal hits failing the MVA CUT
//
// SET plotAllSig = true to display all signal clusters in the bottom right plot
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

TH2F *_hview1[100],*_hview2[100],*_hview3[100],*_hview4[100];

bool plotAllSig = true;

void docluAnaSignalBkg(int entry)
{
   
   cout<<"Inspect entry "<<entry<<endl;

   for (int i=0;i<10;++i) _hview1[i]->Reset(); 
   for (int i=0;i<10;++i) _hview2[i]->Reset(); 
   for (int i=0;i<10;++i) _hview3[i]->Reset(); 
   for (int i=0;i<100;++i) _hview4[i]->Reset(); 
   
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
  

   TFile f("trkDiag.root");
   
   TTree* tree = (TTree*) f.Get("BD/bkgdiag");   
   TTree* tree2 = (TTree*) f.Get("BD/bkgdiag2");   
   fReader.SetTree(tree);      
   fReader2.SetTree(tree2);      
   
   //set fReader2 to correct entry
   int in2(0);
   for (in2=0;in2<tree2->GetEntriesFast();++in2){
       fReader2.SetEntry(in2);       
       if (*iev2 == entry) break;
   }
   
   if (*iev2 != entry) return;




   int cluToWatch(0);        
   for (int in=0;in<tree->GetEntriesFast();++in){
      fReader.SetEntry(in);
      if (*iev != entry) continue;
      
      int ic = *cluIdx;
      if (*nconv > 0 && *mvaout > 0.5) ++cluToWatch;
      if (*nconv > 0 && *mvaout > 0.5) cout<<"Signal problem cluIdx="<<ic<<"   nShits="<<*nshits<<"   nChits="<<*nchits<<"   nConv="<<*nconv<<"   col="<<cluToWatch%10+1<<endl;
      
      for (int id=0;id<*nindex;++id)
      {
         int ih = hitidx[id];
         int ipass = (hitGen[ih]==2) ? 1 : 0;
         _hview1[ipass]->Fill(hitRad[ih],hitPhi[ih]);
         _hview2[ipass]->Fill(hitRad[ih],hitTime[ih]);
         _hview3[ipass]->Fill(hitPhi[ih],hitTime[ih]); 
         if (plotAllSig && *nconv > 0) _hview4[cluToWatch%10+1]->Fill(hitRad[ih],hitPhi[ih]);       
         if (!plotAllSig && *nconv > 0 && *mvaout > 0.5) _hview4[cluToWatch%10+1]->Fill(hitRad[ih],hitPhi[ih]);       
      }
   }
   
    

   
   TCanvas *c1 = new TCanvas("c1","c1",1000,1000);   
   c1->Divide(2,2);
      
   c1->cd(1);   
   _hview1[0]->Draw("");
   _hview1[1]->Draw("same");
   c1->cd(2);   
   _hview2[0]->Draw("");
   _hview2[1]->Draw("same");
   c1->cd(3);   
   _hview3[0]->Draw("");
   _hview3[1]->Draw("same");
   c1->cd(4); 
   _hview4[0]->Draw();
   for (int i=1;i<20;++i) _hview4[i]->Draw("same");
   c1->SaveAs(Form("cluSignalFail_%i.pdf",entry));

   delete c1;
   f.Close();
}


void cluAnaSignalBkg(int entry =-1){

   gStyle->SetOptStat(0);
   gErrorIgnoreLevel = kWarning;
   
   
   for (int i=0;i<10;++i)  _hview1[i]  = new TH2F(Form("hview1_%i",i),      " ",           200, 350, 750, 100, -3.15, 3.15);
   for (int i=0;i<10;++i)  _hview2[i]  = new TH2F(Form("hview2_%i",i),      " ",           200, 350, 750, 200, 400,2000);
   for (int i=0;i<10;++i)  _hview3[i]  = new TH2F(Form("hview3_%i",i),      " ",           100, -3.15, 3.15, 200, 400,2000);
   for (int i=0;i<100;++i) _hview4[i]  = new TH2F(Form("hview4_%i",i),      " ",           200, 350, 750, 100, -3.15, 3.15);

   for (int i=0;i<10;++i){     
     _hview1[i]->GetXaxis()->SetTitle("radius");
     _hview1[i]->GetYaxis()->SetTitle("phi");
     _hview2[i]->GetXaxis()->SetTitle("radius");
     _hview2[i]->GetYaxis()->SetTitle("time");
     _hview3[i]->GetXaxis()->SetTitle("phi");
     _hview3[i]->GetYaxis()->SetTitle("time");
     _hview1[i]->SetMarkerSize(0.6);   
     _hview1[i]->SetMarkerStyle(21);   
     _hview2[i]->SetMarkerSize(0.6);   
     _hview2[i]->SetMarkerStyle(21);   
     _hview3[i]->SetMarkerStyle(21);   
     _hview3[i]->SetMarkerSize(0.6);   
          
     int icol = i%10+1;
     if (icol==10) icol = 92;     
     _hview1[i]->SetMarkerColor(icol);   
     _hview2[i]->SetMarkerColor(icol);   
     _hview3[i]->SetMarkerColor(icol);   
   }
   _hview1[1]->SetMarkerSize(0.8);   
   _hview2[1]->SetMarkerSize(0.8);   
   _hview3[1]->SetMarkerSize(0.8);   
   
   for (int i=0;i<100;++i){
     _hview4[i]->SetMarkerStyle(21);   
     _hview4[i]->SetMarkerSize(0.6);   
     _hview4[i]->GetXaxis()->SetTitle("radius");
     _hview4[i]->GetYaxis()->SetTitle("phi");
     int icol = i%10+1;
     if (icol==10) icol = 92;
      _hview4[i]->SetMarkerColor(icol);   
  
   }     
   
   
   if (entry==-1) for (int i=1;i<500;++i) docluAnaSignalBkg(i);
   else docluAnaSignalBkg(entry);

   
}
