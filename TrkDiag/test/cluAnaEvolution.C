//
// Plot the r,phi,time coordinates of hits in clusters for each cluster and iteration of the cluster algorithm
//
// In you analyzer, you must you must run
//     BD: @local::BD
//    TCD : @local::TCD
//and set 
// physics.producers.FlagBkgHits.SaveBkgClusters :true
// physics.analyzers.TCD.diagLevel :3
//


#include <TFile.h>
#include <TH1F.h>
#include <TH3F.h>
#include <TH2F.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <vector>
#include <iostream>
#include <TStyle.h>
#include <TCanvas.h>


void cluAnaEvolution(int entry = 0, bool do3d=0)
{
   //   int entry = 0;
   
   
   int icol[17]={2,3,4,5,6,7,8,9,13,20,28,32,46,51,67,88,94};
   
   gStyle->SetOptStat(0);
   
   const int NMAX(6000);
   TH2F *_hview1[10000],*_hview2[10000],*_hview3[10000],*_null;   
   for (int i=0;i<NMAX;++i) _hview1[i]  = new TH2F(Form("hview1_%i",i),      " ",           100, -3.15, 3.15,200,400,720);
   for (int i=0;i<NMAX;++i) _hview2[i]  = new TH2F(Form("hview2_%i",i),      " ",           200,400,720, 200, 400,2000);
   for (int i=0;i<NMAX;++i) _hview3[i]  = new TH2F(Form("hview3_%i",i),      " ",           100, -3.15, 3.15, 200, 400,2000);
   _null  = new TH2F("null",      " ",           10,0,1,10,0,1);

   TH3F *_hview4[10000];
   if (do3d)for (int i=0;i<100;++i) _hview4[i] = new TH3F(Form("hview4_%i",i),      " ",           200,400,720,100, -3.15, 3.15, 200, 400,2000);



   

   TTreeReader fReader;  
   TTreeReaderValue<Int_t> nhits = {fReader, "nhits"};
   TTreeReaderArray<Float_t> hitRad = {fReader, "hitRad"};
   TTreeReaderArray<Float_t> hitPhi = {fReader, "hitPhi"};
   TTreeReaderArray<Float_t> hitTime = {fReader, "hitTime"};
   TTreeReaderArray<Int_t> hitNcombo = {fReader, "hitNcombo"};
   TTreeReaderValue<Int_t> ncluIter = {fReader, "ncluIter"};
   TTreeReaderArray<Int_t> cluId = {fReader, "cluId"};
   TTreeReaderArray<Int_t> cluNpass = {fReader, "cluNpass"};
   TTreeReaderArray<Int_t> cnhi = {fReader, "cnhi"};
   TTreeReaderArray<Float_t> cluRad = {fReader, "cluRad"};
   TTreeReaderArray<Float_t> cluPhi = {fReader, "cluPhi"};
   TTreeReaderArray<Float_t> cluTime = {fReader, "cluTime"};
   TTreeReaderArray<Float_t> cluR2diff = {fReader, "cluR2diff"};
   TTreeReaderArray<Float_t> cluRdiff = {fReader, "cluRdiff"};
   TTreeReaderArray<Float_t> cluPdiff = {fReader, "cluPdiff"};
   TTreeReaderArray<Float_t> cluTdiff = {fReader, "cluTdiff"};
   TTreeReaderArray<Float_t> cluTRdiff = {fReader, "cluTRdiff"};
   TTreeReaderValue<Int_t> nhitClu = {fReader, "nhitClu"};
   TTreeReaderArray<Int_t> hcIdxClu = {fReader, "hcIdxClu"};
   TTreeReaderArray<Int_t> hcIdxHit = {fReader, "hcIdxHit"};
   TTreeReaderArray<Int_t> hcNpass = {fReader, "hcNpass"};
   TTreeReaderValue<Int_t> niter = {fReader, "niter"};
   TTreeReaderArray<Int_t> nclu = {fReader, "nclu"};
   TTreeReaderArray<Int_t> nChanged = {fReader, "nChanged"};
   TTreeReaderArray<Float_t> odist = {fReader, "odist"};
   TTreeReaderArray<Float_t> tdist = {fReader, "tdist"};

   

   TFile f("trkDiag.root");
   TTree* tree = (TTree*) f.Get("FlagBkgHits/idiag");   
   fReader.SetTree(tree);      
   fReader.SetEntry(entry);

   
   for (int i=0;i<NMAX;++i){     
     _hview1[i]->GetXaxis()->SetTitle("phi");
     _hview1[i]->GetYaxis()->SetTitle("radius");
     _hview2[i]->GetXaxis()->SetTitle("radius");
     _hview2[i]->GetYaxis()->SetTitle("time");
     _hview3[i]->GetXaxis()->SetTitle("phi");
     _hview3[i]->GetYaxis()->SetTitle("time");
     
     _hview1[i]->SetMarkerColor(icol[i%17]);   
     _hview2[i]->SetMarkerColor(icol[i%17]);   
     _hview3[i]->SetMarkerColor(icol[i%17]);   
     _hview1[i]->SetMarkerSize(0.6);   
     _hview2[i]->SetMarkerSize(0.6);   
     _hview3[i]->SetMarkerSize(0.6);   
     _hview1[i]->SetMarkerStyle(21);   
     _hview2[i]->SetMarkerStyle(21);   
     _hview3[i]->SetMarkerStyle(21);   
   }
   
   for (int i=0;i<100;++i){     
    if (do3d) _hview4[i]->SetMarkerColor(icol[i%17]);   
    if (do3d) _hview4[i]->SetMarkerSize(0.6);   
    if (do3d) _hview4[i]->SetMarkerStyle(21);   
   }
   
   
   std::cout<<"Iteration , number of cluster changed"<<std::endl;
   for (int i=0;i<*niter;++i) cout<<i<<" "<<nChanged[i]<<std::endl;      
    
   for (int in=0;in<*nhitClu;++in){
      int ih = hcIdxHit[in];
      int ipass = 1000*hcNpass[in]+hcIdxClu[in];
      if (ipass>NMAX-1) break;
      //if (hitNcombo[ih]<2) continue;
      
      _hview1[ipass]->Fill(hitPhi[ih],hitRad[ih]);
      _hview2[ipass]->Fill(hitRad[ih],hitTime[ih]);
      _hview3[ipass]->Fill(hitPhi[ih],hitTime[ih]);
      if (do3d && ipass<100) _hview4[ipass]->Fill(hitRad[ih],hitPhi[ih],hitTime[ih]);      
   }
   
   
   
   TCanvas *c1 = new TCanvas("c1","c1",1000,1000);   
   c1->Divide(2,2);
      
   for (int j=0;j<4;++j){
     c1->cd(j+1);   
     _hview1[1000*j]->Draw("");
     for (int i=1000*j+1;i<1000*j+1000;++i) _hview1[i]->Draw("same");
   }
   c1->SaveAs(Form("rphi_%i.pdf",entry));

   for (int j=0;j<4;++j){
     c1->cd(j+1);   
     _hview2[1000*j]->Draw("");
     for (int i=1000*j+1;i<1000*j+1000;++i) _hview2[i]->Draw("same");
   }
   c1->SaveAs(Form("rt_%i.pdf",entry));

   for (int j=0;j<4;++j){
     c1->cd(j+1);   
     _hview3[1000*j]->Draw("");
     for (int i=100*j+1;i<1000*j+1000;++i) _hview3[i]->Draw("same");
   }
   c1->SaveAs(Form("phit_%i.pdf",entry));

   for (int j=0;j<3;++j){
      c1->cd(1);   
      _hview1[1000*j]->Draw("");
      for (int i=1;i<1000*j+1000;++i) _hview1[i]->Draw("same");
      c1->cd(2);   
      _hview2[1000*j]->Draw("");
      for (int i=1;i<1000*j+1000;++i) _hview2[i]->Draw("same");
      c1->cd(3);   
      _hview3[1000*j]->Draw("");
      for (int i=1;i<1000*j+1000;++i) _hview3[i]->Draw("same");
      c1->cd(4);   
      _null->Draw("");
      c1->SaveAs(Form("All_%i_%i.pdf",j,entry));
   }

   if (do3d)
   {   
     TCanvas *c2 = new TCanvas("c2","c2",1000,1000);   
     _hview4[0]->Draw();
     for (int i=1;i<100;++i) _hview4[i]->Draw("same");
     c2->SaveAs(Form("init3d_%i.pdf",entry));
   }

   f.Close();
}
