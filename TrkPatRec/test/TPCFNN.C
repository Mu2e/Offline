// Macro to produce data to train NN for TimeAndPhiClusterFinder module
// must run th emodule with diag=1 to produce the input file.

#include <TFile.h>
#include <TH1.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TStyle.h>
#include <vector>
#include <iostream>
#include <cmath>


int  NHitCut_(15);




void TPCFNN(){

   TFile saveData("data.root","RECREATE");
   TTree dTuple("Data","Data");
   float var2,var3,var4,var5,var6,var7,var8,var9;
   int var1;
   TBranch *b_var1 = dTuple.Branch("var1",   &var1,   "var1/I");
   TBranch *b_var2 = dTuple.Branch("var2",   &var2,   "var2/F");
   TBranch *b_var3 = dTuple.Branch("var3",   &var3,   "var3/F");
   TBranch *b_var4 = dTuple.Branch("var4",   &var4,   "var4/F");
   TBranch *b_var5 = dTuple.Branch("var5",   &var5,   "var5/F");
   TBranch *b_var6 = dTuple.Branch("var6",   &var6,   "var6/F");
   TBranch *b_var7 = dTuple.Branch("var7",   &var7,   "var7/F");
   TBranch *b_var8 = dTuple.Branch("var8",   &var8,   "var8/F");
   TBranch *b_var9 = dTuple.Branch("var9",   &var9,   "var9/F");



   // TimeClusterDiag reader
   TTreeReader fReader;
   TTreeReaderValue<Int_t>   iev     = {fReader, "iev"};
   TTreeReaderValue<Int_t>   Nch     = {fReader, "Nch"};
   TTreeReaderArray<Int_t>   chSel   = {fReader, "chSel"};
   TTreeReaderArray<Int_t>   chNhit  = {fReader, "chNhit"};
   TTreeReaderArray<Float_t> chTime  = {fReader, "chTime"};
   TTreeReaderArray<Float_t> chRad   = {fReader, "chRad"};
   TTreeReaderArray<Float_t> chPhi   = {fReader, "chPhi"};
   TTreeReaderArray<Float_t> chX     = {fReader, "chX"};
   TTreeReaderArray<Float_t> chY     = {fReader, "chY"};
   TTreeReaderArray<Float_t> chZ     = {fReader, "chZ"};
   TTreeReaderArray<Int_t>   chPdg   = {fReader, "chPdg"};
   TTreeReaderArray<Int_t>   chCrCode= {fReader, "chCrCode"};
   TTreeReaderArray<Int_t>   chSimId = {fReader, "chSimId"};
   TTreeReaderValue<Int_t>   nhit1   = {fReader, "nhit1"};
   TTreeReaderArray<Int_t>   hitIdx1 = {fReader, "hitIdx1"};
   TTreeReaderArray<Int_t>   nclu1   = {fReader, "nclu1"};
   TTreeReaderValue<Int_t>   nhit2   = {fReader, "nhit2"};
   TTreeReaderArray<Int_t>   hitIdx2 = {fReader, "hitIdx2"};
   TTreeReaderArray<Int_t>   nclu2   = {fReader, "nclu2"};



   TFile f("trkDiag1BBBest.root");
   if (!f.IsOpen()) return   ;
   TTree* tree = (TTree*) f.Get("TimeClusterFinderDe/tpcdiag");
   fReader.SetTree(tree);


   for (int entry=0;entry<tree->GetEntries();++entry){
       fReader.SetEntry(entry);

       // FILL THE DATA AND REQUIRE AT LEAST ONE CONVERSION TRACK
       map<int, vector<int>> cluHitMap;
       for (int i=0;i<*nhit1;++i) cluHitMap[nclu1[i]].emplace_back(hitIdx1[i]);

       int convCluId(-1), nMatch(0),nHits(0);
       for (auto& kv : cluHitMap) {
          int nCE(0),nh(0);
          for (auto hit : kv.second) {nh += chNhit[hit];if (chCrCode[hit]==167) nCE += chNhit[hit];}
          if (nCE>nMatch) {nMatch = nCE;convCluId = kv.first;nHits=nh;}
       }

       if (convCluId<0 || nHits < NHitCut_) continue;

       //FILL DATA FOR NN
       auto& hits = cluHitMap[convCluId];
       float dr(0),dt(0),dx(0),dy(0),dz(0);
       for (auto hit : hits) {dr += chRad[hit];dx += chRad[hit]*cos(chPhi[hit]);dy += chRad[hit]*sin(chPhi[hit]);dz += chZ[hit];dt += chTime[hit];}
       dr /= hits.size();
       dt /= hits.size();
       dx /= hits.size();
       dy /= hits.size();
       dz /= hits.size();
       float dp = atan2(dy,dx);

       for (const auto& hit : hits){
          float dphi = chPhi[hit]- dp;
          if (dphi > 3.1415)  dphi -= 6.2830;
          if (dphi < -3.1415) dphi += 6.2830;
          var1 = chCrCode[hit]==167 ? 1 : 0;
          var2 = chRad[hit]-dr;
          var3 = dphi;
          var4 = chTime[hit]-dt;
          var5 = chZ[hit]-dz;
          var6 = chNhit[hit];


          dTuple.Fill();
       }
   }

   saveData.Write();
   saveData.Close();
}
