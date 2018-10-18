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

void cluAnaIdea()
{
      
   TH1F* _nSig = new TH1F("nSig","nSig",20,0,20);
   TH1F* _nBkg = new TH1F("nBkg","nBkg",20,0,20);
   TH1F* _nRecS = new TH1F("nRecS","nSigH",20,0,20);
   TH1F* _nRecB = new TH1F("nRecB","nBkgH",20,0,20);
   
   for (int i=0;i<10;++i)  _hview1[i]  = new TH2F(Form("hview1_%i",i),      " ",           200, 400, 720, 100, -3.15, 3.15);
   for (int i=0;i<10;++i)  _hview2[i]  = new TH2F(Form("hview2_%i",i),      " ",           200, 400, 720, 200, 400,2000);
   for (int i=0;i<10;++i)  _hview3[i]  = new TH2F(Form("hview3_%i",i),      " ",           100, -3.15, 3.15, 200, 400,2000);

   for (int i=0;i<10;++i){     
     _hview1[i]->GetXaxis()->SetTitle("radius");
     _hview1[i]->GetYaxis()->SetTitle("phi");
     _hview2[i]->GetXaxis()->SetTitle("radius");
     _hview2[i]->GetYaxis()->SetTitle("time");
     _hview3[i]->GetXaxis()->SetTitle("phi");
     _hview3[i]->GetYaxis()->SetTitle("time");
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
   TTreeReaderValue<Float_t> HitRho          = {fReader, "HitRho"};
   TTreeReaderValue<Float_t> HitRhoSpread    = {fReader, "HitRhoSpread"};
   TTreeReaderValue<Float_t> ClusterRho      = {fReader, "ClusterRho"};
   TTreeReaderValue<Float_t> TimeSpread      = {fReader, "TimeSpread"};
   TTreeReaderValue<Float_t> ZMin            = {fReader, "ZMin"};
   TTreeReaderValue<Float_t> ZMax            = {fReader, "ZMax"};
   TTreeReaderValue<Float_t> ZGap            = {fReader, "ZGap"};
   TTreeReaderValue<Float_t> NPlanes         = {fReader, "NPlanes"};
   TTreeReaderValue<Float_t> NExpectedPlanes = {fReader, "NExpectedPlanes"};
   TTreeReaderValue<Float_t> PlaneFraction   = {fReader, "PlaneFraction"};
   TTreeReaderValue<Float_t> NPlaneHits      = {fReader, "NPlaneHits"};
   TTreeReaderValue<Float_t> NHits           = {fReader, "NHits"};

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
   


   for (int entry = 1; entry < tree2->GetEntriesFast();++entry){

      //set fReader2 to correct entry
      int in2(0);
      for (in2=0;in2<tree2->GetEntriesFast();++in2){
          fReader2.SetEntry(in2);       
          if (*iev2 == entry) break;
      }

      if (*iev2 != entry) continue;

      vector<double> badCluTime,badCluPhi;  
      for (int in=0;in<tree->GetEntriesFast();++in){
         fReader.SetEntry(in);
         if (*iev != entry) continue;
         if (*nshits<3 || *NPlanes<2) continue; 
         bool sigSel = (*NHits <25) && (*ZGap < 50 || *ZGap >800) && (*PlaneFraction< 0.41 || *PlaneFraction >0.95);

         if (sigSel && *mvaout >0.001) badCluTime.push_back(hitTime[hitidx[0]]);
         if (sigSel && *mvaout >0.001) badCluPhi.push_back(hitPhi[hitidx[0]]);
         //if (*mvaout < 0.5  && *mvaout >0.001 && *NHits < 15)  badCluTime.insert(hitTime[hitidx[0]]);
      }
      int nRecS(0),nRecB(0);
      for (int in=0;in<tree->GetEntriesFast();++in){
         fReader.SetEntry(in);
         if (*iev != entry) continue;
         if (*nshits<3 || *NPlanes<2) continue; 

         bool sigSel = (*NHits <25) && (*ZGap < 50 || *ZGap >800) && (*PlaneFraction< 0.41 || *PlaneFraction >0.95);
         if (sigSel || *NHits > 15) continue;


         int nCluOk(0);
         for (size_t i=0;i<badCluTime.size();++i){ 
           double dPhi = badCluPhi[i]-hitPhi[hitidx[0]];
           while (dPhi < -3.14159)dPhi +=6.283185;
           while (dPhi > 3.14159) dPhi -=6.283185;           
           if ((std::abs(badCluTime[i]-hitTime[hitidx[0]]) < 30) && abs(dPhi) < 1.0) ++nCluOk;// +=badCluN[i];
         }
         if (nCluOk < 2 && *nconv>0) cout<<*nshits<<" "<<*nconv<<" "<<hitTime[hitidx[0]]<<" "<<hitPhi[hitidx[0]]<<"    "<<entry<<"   "<<nCluOk<<endl;   

         if (*nconv>0) _nSig->Fill(nCluOk);
         else          _nBkg->Fill(nCluOk);
         
         if (*nconv>0  && nCluOk>2) ++nRecS;
         if (*nconv==0 && nCluOk>2) ++nRecB;         
      }
      
      if (nRecS)_nRecS->Fill(nRecS);
      _nRecB->Fill(nRecB);
   }
   



   TCanvas *c1 = new TCanvas("c1","c1",1000,1000);   
   c1->Divide(2,2);
      
   c1->cd(1);   
   _nSig->Draw();
   c1->cd(2);   
   _nBkg->Draw();
   c1->cd(3);   
   _nRecS->Draw();
   c1->cd(4);   
   _nRecB->Draw();
   
   c1->SaveAs("clu.pdf");

   delete c1;
   f.Close();
}
