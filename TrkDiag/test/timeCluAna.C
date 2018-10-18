// Plot the r,phi,time coordinates of hits in timeClusters (black = bkg, red = CE) for timeClusters having at least one CE hit
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



void dotimeCluAna(int entry)
{
   
   //cout<<"Inspect entry "<<entry<<endl;

   for (int i=0;i<10;++i) _hview1[i]->Reset(); 
   for (int i=0;i<10;++i) _hview2[i]->Reset(); 
   for (int i=0;i<10;++i) _hview3[i]->Reset(); 
   for (int i=0;i<10;++i) _hview4[i]->Reset(); 
  
   TTreeReader fReader;  
   TTreeReaderValue<Int_t>   iev        = {fReader, "iev"};
   TTreeReaderValue<Int_t>   tcnHits    = {fReader, "tcnHits"};
   TTreeReaderArray<Int_t>   tcHitIdx   = {fReader, "tcHitIdx"};
   TTreeReaderValue<Int_t>   tcover     = {fReader, "tcover"};
   TTreeReaderValue<Int_t>   tcnce      = {fReader, "tcnce"};
   TTreeReaderValue<Float_t> tctime     = {fReader, "tctime"};
   TTreeReaderValue<Float_t> tcphi      = {fReader, "tcphi"};
   TTreeReaderValue<Float_t> tcminhtime = {fReader, "tcminhtime"};
   TTreeReaderValue<Float_t> tcmaxhtime = {fReader, "tcmaxhtime"};

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
   
   TTree* tree = (TTree*) f.Get("TCD/tcdiag2");   
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
   
   
   int nCEtot(0);
   for (int in=0;in<tree->GetEntriesFast();++in)
   {
      fReader.SetEntry(in);
      if (*iev != entry) continue;      
      for (int id=0;id<*tcnHits;++id)
      {
         int ih = tcHitIdx[id];
         if (hitGen[ih]==2) ++nCEtot;
      }
   }


   for (int in=0;in<tree->GetEntriesFast();++in)
   {
      fReader.SetEntry(in);
      if (*iev != entry) continue;
      
      int nCe(0);
      for (int id=0;id<*tcnHits;++id) if (hitGen[tcHitIdx[id]]==2) ++nCe;
      if (nCe>0 ) cout<<"Entry "<<*iev<<"   CE hits "<<nCe<<" / "<<nCEtot<<"   "<<*tcnHits<<endl;
                  
      if (float(nCe)/float(nCEtot) < 0.2) continue; 
      
      for (int id=0;id<*tcnHits;++id)
      {
         int ih = tcHitIdx[id];
         int ipass = (hitGen[ih]==2) ? 1 : 0;
         _hview1[ipass]->Fill(hitRad[ih],hitPhi[ih]);
         _hview2[ipass]->Fill(hitRad[ih],hitTime[ih]);
         _hview3[ipass]->Fill(hitPhi[ih],hitTime[ih]); 
         _hview4[ipass]->Fill(hitRad[ih]*cos(hitPhi[ih]),hitRad[ih]*sin(hitPhi[ih]));
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
   _hview4[0]->Draw("");
   _hview4[1]->Draw("same");
   
   c1->SaveAs(Form("tcCluAna_%i.pdf",entry));

   delete c1;
   f.Close();
}


void timeCluAna(int entry =-1){

   gStyle->SetOptStat(0);
   gErrorIgnoreLevel = kWarning;
   
   
   for (int i=0;i<10;++i)  _hview1[i]  = new TH2F(Form("hview1_%i",i),      " ",  200, 400, 720, 100, -3.15, 3.15);
   for (int i=0;i<10;++i)  _hview2[i]  = new TH2F(Form("hview2_%i",i),      " ",  200, 400, 720, 200, 400,2000);
   for (int i=0;i<10;++i)  _hview3[i]  = new TH2F(Form("hview3_%i",i),      " ",  100, -3.15, 3.15, 200, 400,2000);
   for (int i=0;i<10;++i)  _hview4[i]  = new TH2F(Form("hview4_%i",i),      " ",  200, -700,700,200,-700,700);

   for (int i=0;i<10;++i){     
     _hview1[i]->GetXaxis()->SetTitle("radius");
     _hview1[i]->GetYaxis()->SetTitle("phi");
     _hview2[i]->GetXaxis()->SetTitle("radius");
     _hview2[i]->GetYaxis()->SetTitle("time");
     _hview3[i]->GetXaxis()->SetTitle("phi");
     _hview3[i]->GetYaxis()->SetTitle("time");
     _hview4[i]->GetXaxis()->SetTitle("x");
     _hview4[i]->GetYaxis()->SetTitle("y");
     _hview1[i]->SetMarkerSize(0.6);   
     _hview1[i]->SetMarkerStyle(21);   
     _hview2[i]->SetMarkerSize(0.6);   
     _hview2[i]->SetMarkerStyle(21);   
     _hview3[i]->SetMarkerStyle(21);   
     _hview3[i]->SetMarkerSize(0.6);   
     _hview4[i]->SetMarkerStyle(21);   
     _hview4[i]->SetMarkerSize(0.6);   
          
     int icol = i%10+1;
     if (icol==10) icol = 92;     
     _hview1[i]->SetMarkerColor(icol);   
     _hview2[i]->SetMarkerColor(icol);   
     _hview3[i]->SetMarkerColor(icol);   
     _hview4[i]->SetMarkerColor(icol);   
   }
   

   
   if (entry==-1) for (int i=1;i<700;++i) dotimeCluAna(i);
   else dotimeCluAna(entry);

   
}
