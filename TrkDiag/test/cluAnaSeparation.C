// Plot the r,phi,time coordinates of hits in clusters that have one [article contributing to it (black)
// or more than one (red) 
// Bottom-right plot the reconstructed clusters that have more than one particle contributing to it. 
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


void cluAnaSeparation(int entry = 0)
{
   
   //   int entry = 0;
   
   gStyle->SetOptStat(0);
   
   TH2F *_hview1[100],*_hview2[100],*_hview3[100],*_hview4[100];   
   for (int i=0;i<100;++i) _hview1[i]  = new TH2F(Form("hview1_%i",i),      " ",           200, 400, 720, 100, -3.15, 3.15);
   for (int i=0;i<100;++i) _hview2[i]  = new TH2F(Form("hview2_%i",i),      " ",           200, 400, 720, 200, 400,2000);
   for (int i=0;i<100;++i) _hview3[i]  = new TH2F(Form("hview3_%i",i),      " ",           100, -3.15, 3.15, 200, 400,2000);
   for (int i=0;i<100;++i) _hview4[i]  = new TH2F(Form("hview4_%i",i),      " ",           100, -3.15, 3.15, 200, 400,2000);

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
   cout<<"Analyze evt "<<*iev2<<endl;
   
   
   
   std::set<int> ibToWatch;    
   for (int in=0;in<tree->GetEntriesFast();++in){
      fReader.SetEntry(in);
      if (*iev != entry) continue;
      
      if (*ncontrib>1) for (int ib=0;ib<*ncontrib;++ib) ibToWatch.insert(icontrib[ib]);
            
      for (int id=0;id<*nindex;++id){
         int ih = hitidx[id];
         int ipass = (*ncontrib>1) ? 1 : 0;
         _hview1[ipass]->Fill(hitRad[ih],hitPhi[ih]);
         _hview2[ipass]->Fill(hitRad[ih],hitTime[ih]);
         _hview3[ipass]->Fill(hitPhi[ih],hitTime[ih]);        
      }
   }
   
  
  
   //cout<<"Now the clusters containing those hits"<<endl;
   for (int in=0;in<tree->GetEntriesFast();++in){
      fReader.SetEntry(in);
      if (*iev != entry) continue;
     
      bool reject(true);
      for (int ib=0;ib<*ncontrib;++ib) if (ibToWatch.find(icontrib[ib]) != ibToWatch.end()) reject=false;      
      if (reject) continue;

      for (int id=0;id<*nindex;++id){
        int ih = hitidx[id];
        _hview4[*cluIdx%10+1]->Fill(hitPhi[ih],hitTime[ih]);
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
   for (int i=1;i<22;++i)_hview4[i]->Draw("same");  
   c1->SaveAs(Form("cluSeparation_%i.pdf",entry));


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
