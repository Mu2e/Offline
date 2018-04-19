// Plot the r,phi,time coordinates of signal hits, separating the pass/fail in mva if the signal hit content > 20%
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


void cluAnaContent(int entry = 0)
{
   
   //   int entry = 0;
   
   gStyle->SetOptStat(0);
   
   TH2F *_hview1[100],*_hview2[100],*_hview3[100],*_hview4[100];   
   for (int i=0;i<100;++i) _hview1[i]  = new TH2F(Form("hview1_%i",i),      " ",           200, 400, 720, 100, -3.15, 3.15);
   for (int i=0;i<100;++i) _hview2[i]  = new TH2F(Form("hview2_%i",i),      " ",           200, 400, 720, 200, 400,2000);
   for (int i=0;i<100;++i) _hview3[i]  = new TH2F(Form("hview3_%i",i),      " ",           100, -3.15, 3.15, 200, 400,2000);
   for (int i=0;i<100;++i) _hview4[i]  = new TH2F(Form("hview4_%i",i),      " ",           100, -3.15, 3.15, 200, 400,2000);

   TH1F *_valT[3],*_valR[3],*_valP[3];
   for (int i=0;i<3;++i)_valT[i] = new TH1F(Form("valT%i",i),"valT",100,0,50);
   for (int i=0;i<3;++i)_valR[i] = new TH1F(Form("valR%i",i),"valR",100,0,100);
   for (int i=0;i<3;++i)_valP[i] = new TH1F(Form("valP%i",i),"valP",100,0,1);


   

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
   
   TTreeReader fReader2;  
   TTreeReaderValue<Int_t> iev      = {fReader2, "iev"};
   TTreeReaderValue<Int_t> cluIdx   = {fReader2, "cluIdx"};
   TTreeReaderValue<Int_t> nchits   = {fReader2, "nchits"};
   TTreeReaderValue<Int_t> nshits   = {fReader2, "nshits"};
   TTreeReaderValue<Float_t> mvaout = {fReader2, "mvaout"};
   TTreeReaderValue<Float_t> pmom   = {fReader2, "pmom"};
   TTreeReaderValue<Int_t> ppid     = {fReader2, "ppid"};
   TTreeReaderValue<Int_t> ppdg     = {fReader2, "ppdg"};
   TTreeReaderValue<Int_t> pgen     = {fReader2, "pgen"};
   TTreeReaderValue<Int_t> pproc    = {fReader2, "pproc"};
   TTreeReaderValue<Int_t> nprimary = {fReader2, "nprimary"};
   TTreeReaderValue<Int_t> nconv    = {fReader2, "nconv"};
   TTreeReaderValue<Int_t> ndelta   = {fReader2, "ndelta"};
   TTreeReaderValue<Int_t> ncompt   = {fReader2, "ncompt"};
   TTreeReaderValue<Int_t> ngconv   = {fReader2, "ngconv"};
   TTreeReaderValue<Int_t> nebkg    = {fReader2, "nebkg"};
   TTreeReaderValue<Int_t> nprot    = {fReader2, "nprot"};
   TTreeReaderValue<Int_t> ncontrib = {fReader2, "ncontrib"};
   TTreeReaderArray<Int_t> icontrib = {fReader2, "icontrib"};
   TTreeReaderValue<Int_t> nindex   = {fReader2, "nindex"};
   TTreeReaderArray<Int_t> hitidx   = {fReader2, "hitidx"};

   TTreeReader fReader3;  
   TTreeReaderValue<Int_t> nhits3    = {fReader3, "nhits"};
   TTreeReaderArray<Int_t> hitPdg    = {fReader3, "hitPdg"};
   TTreeReaderArray<Int_t> hitCrCode = {fReader3, "hitCrCode"};
   TTreeReaderArray<Int_t> hitGen    = {fReader3, "hitGen"};
  

   TFile f("trkDiag.root");
   TTree* tree = (TTree*) f.Get("FlagBkgHits/idiag");   
   fReader.SetTree(tree);      
   fReader.SetEntry(entry);
   
   TTree* tree2 = (TTree*) f.Get("BD/bkgdiag");   
   fReader2.SetTree(tree2);      

   TTree* tree3 = (TTree*) f.Get("BD/bkgdiag2");   
   fReader3.SetTree(tree3);      
   fReader3.SetEntry(entry);
   
   
   for (int i=0;i<100;++i){     
     _hview1[i]->GetXaxis()->SetTitle("radius");
     _hview1[i]->GetYaxis()->SetTitle("phi");
     _hview2[i]->GetXaxis()->SetTitle("radius");
     _hview2[i]->GetYaxis()->SetTitle("time");
     _hview3[i]->GetXaxis()->SetTitle("phi");
     _hview3[i]->GetYaxis()->SetTitle("time");
     _hview4[i]->GetXaxis()->SetTitle("phi");
     _hview4[i]->GetYaxis()->SetTitle("time");
     int icol = i%10+1;
     if (icol==10) icol = 92;
     
     _hview1[i]->SetMarkerColor(icol);   
     _hview2[i]->SetMarkerColor(icol);   
     _hview3[i]->SetMarkerColor(icol);   
     _hview4[i]->SetMarkerColor(icol);   
     _hview1[i]->SetMarkerSize(0.6);   
     _hview2[i]->SetMarkerSize(0.6);   
     _hview3[i]->SetMarkerSize(0.6);   
     _hview4[i]->SetMarkerSize(0.6);   
     _hview1[i]->SetMarkerStyle(21);   
     _hview2[i]->SetMarkerStyle(21);   
     _hview3[i]->SetMarkerStyle(21);   
     _hview4[i]->SetMarkerStyle(21);   
   }
   
   _valT[0]->SetLineColor(1);
   _valT[1]->SetLineColor(2);
   _valT[2]->SetLineColor(4);
   _valR[0]->SetLineColor(1);
   _valR[1]->SetLineColor(2);
   _valR[2]->SetLineColor(4);
   _valP[0]->SetLineColor(1);
   _valP[1]->SetLineColor(2);
   _valP[2]->SetLineColor(4);
   

    //map all the SimParticles to cluIdx
    map<int,vector<int>> contmap;
    for (int in=0;in<tree2->GetEntriesFast();++in){
      fReader2.SetEntry(in);
      if (*iev != (entry+1)) continue;
      
      for (size_t ib=0;ib <*ncontrib;++ib) 
        contmap[icontrib[ib]].push_back(in);
    }
       
   //map cluster to value to plot   
   map<int,float> contmapT,contmapR,contmapP;     
   for (int ic=0;ic<*ncluIter;++ic){
      contmapT[cluId[ic]] = cluTdiff[ic];
      contmapR[cluId[ic]] = cluRdiff[ic];
      contmapP[cluId[ic]] = cluPdiff[ic];
   }
    

    //check if a cluIdx has all hits from a SimParticles, or if another one shares hits
    
    for (int in=0;in<tree2->GetEntriesFast();++in){
      fReader2.SetEntry(in);
      if (*iev != (entry+1)) continue;               
      if (*nindex<2) continue; 
      if (*ncontrib==1 && contmap[icontrib[0]].size()==1) _valT[0]->Fill(contmapT[*cluIdx]);        
      if (*ncontrib==1 && contmap[icontrib[0]].size()>1)  _valT[1]->Fill(contmapT[*cluIdx]);     
      if (*ncontrib>1)                                    _valT[2]->Fill(contmapT[*cluIdx]);
      if (*ncontrib==1 && contmap[icontrib[0]].size()==1) _valR[0]->Fill(contmapR[*cluIdx]);        
      if (*ncontrib==1 && contmap[icontrib[0]].size()>1)  _valR[1]->Fill(contmapR[*cluIdx]);     
      if (*ncontrib>1)                                    _valR[2]->Fill(contmapR[*cluIdx]);
      if (*ncontrib==1 && contmap[icontrib[0]].size()==1) _valP[0]->Fill(contmapP[*cluIdx]);        
      if (*ncontrib==1 && contmap[icontrib[0]].size()>1)  _valP[1]->Fill(contmapP[*cluIdx]);     
      if (*ncontrib>1)                                    _valP[2]->Fill(contmapP[*cluIdx]);
    }
      
      
      
    for (auto& kv : contmap)
    {               

      if (kv.second.size()>1) continue;

      int icol = kv.first%10+1;        
      for (int in : kv.second) {

        fReader2.SetEntry(in);
        for (int id=0;id<*nindex;++id){
           int ih = hitidx[id];
           _hview1[icol]->Fill(hitRad[ih],hitPhi[ih]);
           _hview2[icol]->Fill(hitRad[ih],hitTime[ih]);
           _hview3[icol]->Fill(hitPhi[ih],hitTime[ih]);
        }
      }
    }      

    
    
    


  
  
   
   TCanvas *c1 = new TCanvas("c1","c1",1000,1000);   
   c1->Divide(2,2);
   
   /*   
    c1->cd(1);   
   _hview1[0]->Draw("");
   for(int i=1;i<10;++i) _hview1[i]->Draw("same");
    c1->cd(2);   
   _hview2[0]->Draw("");
   for(int i=1;i<10;++i) _hview2[i]->Draw("same");
    c1->cd(3);   
   _hview3[0]->Draw("");
   for(int i=1;i<10;++i) _hview3[i]->Draw("same");
   */
    c1->cd(1);   
    _valT[1]->Draw();
    _valT[0]->Draw("same");
    _valT[2]->Draw("same");
    c1->cd(2);   
    _valR[1]->Draw();
    _valR[0]->Draw("same");
    _valR[2]->Draw("same");
    c1->cd(3);   
    _valP[1]->Draw();
    _valP[0]->Draw("same");
    _valP[2]->Draw("same");
    c1->SaveAs(Form("AllContent_%i.pdf",entry));

   f.Close();
}
/*
   TTreeReader fReader2;  
   TTreeReaderValue<Int_t> iev = {fReader2, "iev"};
   TTreeReaderValue<Double_t> dx = {fReader2, "cpos.dx"};
   TTreeReaderValue<Double_t> dy = {fReader2, "cpos.dy"};
   TTreeReaderValue<Double_t> dz = {fReader2, "cpos.dz"};
   TTreeReaderValue<Float_t> ctime2 = {fReader2, "ctime"};
   TTreeReaderValue<Char_t> isbkg = {fReader2, "isbkg"};
   TTreeReaderValue<Char_t> isref = {fReader2, "isref"};
   TTreeReaderValue<Char_t> isolated = {fReader2, "isolated"};
   TTreeReaderValue<Char_t> stereo = {fReader2, "stereo"};
   TTreeReaderValue<Float_t> mindt = {fReader2, "mindt"};
   TTreeReaderValue<Float_t> mindrho = {fReader2, "mindrho"};
   TTreeReaderValue<Int_t> nchits = {fReader2, "nchits"};
   TTreeReaderValue<Int_t> nshits = {fReader2, "nshits"};
   TTreeReaderValue<Int_t> nactive = {fReader2, "nactive"};
   TTreeReaderValue<Int_t> nstereo = {fReader2, "nstereo"};
   TTreeReaderValue<Int_t> nsactive = {fReader2, "nsactive"};
   TTreeReaderValue<Int_t> nbkg = {fReader2, "nbkg"};
   TTreeReaderValue<Float_t> HitRho = {fReader2, "HitRho"};
   TTreeReaderValue<Float_t> HitRhoSpread = {fReader2, "HitRhoSpread"};
   TTreeReaderValue<Float_t> ClusterRho = {fReader2, "ClusterRho"};
   TTreeReaderValue<Float_t> TimeSpread = {fReader2, "TimeSpread"};
   TTreeReaderValue<Float_t> ZMin = {fReader2, "ZMin"};
   TTreeReaderValue<Float_t> ZMax = {fReader2, "ZMax"};
   TTreeReaderValue<Float_t> ZGap = {fReader2, "ZGap"};
   TTreeReaderValue<Float_t> NPlanes = {fReader2, "NPlanes"};
   TTreeReaderValue<Float_t> NExpectedPlanes = {fReader2, "NExpectedPlanes"};
   TTreeReaderValue<Float_t> PlaneFraction = {fReader2, "PlaneFraction"};
   TTreeReaderValue<Float_t> NPlaneHits = {fReader2, "NPlaneHits"};
   TTreeReaderValue<Float_t> NHits = {fReader2, "NHits"};
   TTreeReaderValue<Float_t> StereoFraction = {fReader2, "StereoFraction"};
   TTreeReaderValue<Float_t> mvaout = {fReader2, "mvaout"};
   TTreeReaderValue<Int_t> mvastat = {fReader2, "mvastat"};
   TTreeReaderValue<Float_t> pmom = {fReader2, "pmom"};
   TTreeReaderValue<Int_t> ppid = {fReader2, "ppid"};
   TTreeReaderValue<Int_t> ppdg = {fReader2, "ppdg"};
   TTreeReaderValue<Int_t> pgen = {fReader2, "pgen"};
   TTreeReaderValue<Int_t> pproc = {fReader2, "pproc"};
   TTreeReaderValue<Int_t> nprimary = {fReader2, "nprimary"};
   TTreeReaderValue<Int_t> nconv = {fReader2, "nconv"};
   TTreeReaderValue<Int_t> ndelta = {fReader2, "ndelta"};
   TTreeReaderValue<Int_t> ncompt = {fReader2, "ncompt"};
   TTreeReaderValue<Int_t> ngconv = {fReader2, "ngconv"};
   TTreeReaderValue<Int_t> nebkg = {fReader2, "nebkg"};
   TTreeReaderValue<Int_t> nprot = {fReader2, "nprot"};
*/
