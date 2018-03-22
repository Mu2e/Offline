//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Feb  8 15:32:58 2013 by ROOT version 5.30/02
// from TTree Calo/Calo
// found on file: test.root
//////////////////////////////////////////////////////////

#ifndef AnalysisBaseMatch_h
#define AnalysisBaseMatch_h

#include <vector>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH2.h>
#include <TH2Poly.h>
#include <TGraph.h>
#include <TVector3.h>
#include "Src/Disk.hh"


class AnalysisBaseMatch : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Int_t           evt;
   Int_t           nMatch;
   Int_t           mTrkId[1024];   //[nMatch]
   Int_t           mCluId[1024];   //[nMatch]
   Float_t         mChi2[1024];   //[nMatch]
   Float_t         mChi2Pos[1024];   //[nMatch]
   Float_t         mChi2Time[1024];   //[nMatch]
   Int_t           nTrk;
   Float_t         trkX[1024];   //[nTrk]
   Float_t         trkY[1024];   //[nTrk]
   Float_t         trkZ[1024];   //[nTrk]
   Float_t         trkFFX[1024];   //[nTrk]
   Float_t         trkFFY[1024];   //[nTrk]
   Float_t         trkFFZ[1024];   //[nTrk]
   Float_t         trkCposX[1024];   //[nTrk]
   Float_t         trkCposY[1024];   //[nTrk]
   Float_t         trkt[1024];   //[nTrk]
   Float_t         trke[1024];   //[nTrk]
   Float_t         trkpX[1024];   //[nTrk]
   Float_t         trkpY[1024];   //[nTrk]
   Float_t         trkpZ[1024];   //[nTrk]
   Int_t           trknHit[1024];   //[nTrk]
   Int_t           trkStat[1024];   //[nTrk]
   Float_t         trkprob[1024];   //[nTrk]
   Float_t         trkd0[1024];   //[nTrk]
   Float_t         trkz0[1024];   //[nTrk]
   Float_t         trkphi0[1024];   //[nTrk]
   Float_t         trkomega[1024];   //[nTrk]
   Float_t         trkcdip[1024];   //[nTrk]
   Float_t         trkdlen[1024];   //[nTrk]
   Int_t           trkCluIdx[1024];   //[nTrk]
   Int_t           nCluster;
   Float_t         cluEnergy[1024];   //[nCluster]
   Float_t         cluTime[1024];   //[nCluster]
   Float_t         cluCogX[1024];   //[nCluster]
   Float_t         cluCogY[1024];   //[nCluster]
   Float_t         cluCogZ[1024];   //[nCluster]
   Float_t         cluE1[1024];   //[nCluster]
   Float_t         cluE9[1024];   //[nCluster]
   Float_t         cluE25[1024];   //[nCluster]
   Float_t         clu2Mom[1024];   //[nCluster]
   Float_t         cluAngle[1024];   //[nCluster]
   Int_t           cluConv[1024];   //[nCluster]
   Int_t           cluSimIdx[1024];   //[nCluster]
   Int_t           cluSimLen[1024];   //[nCluster]
   vector<vector<int> > *cluList;
   Int_t           nCluSim;
   Int_t           clusimId[1024];   //[nCluSim]
   Int_t           clusimPdgId[1024];   //[nCluSim]
   Int_t           clusimGenIdx[1024];   //[nCluSim]
   Int_t           clusimCrCode[1024];   //[nCluSim]
   Float_t         clusimMom[1024];   //[nCluSim]
   Float_t         clusimPosX[1024];   //[nCluSim]
   Float_t         clusimPosY[1024];   //[nCluSim]
   Float_t         clusimPosZ[1024];   //[nCluSim]
   Float_t         clusimTime[1024];   //[nCluSim]
   Float_t         clusimEdep[1024];   //[nCluSim]
   Int_t           nCry;
   Int_t           cryId[101];   //[nCry]
   Int_t           crySecId[101];   //[nCry]
   Float_t         cryPosX[101];   //[nCry]
   Float_t         cryPosY[101];   //[nCry]
   Float_t         cryPosZ[101];   //[nCry]
   Float_t         cryEdep[101];   //[nCry]
   Float_t         cryTime[101];   //[nCry]
   Int_t           nVd;
   Int_t           vdId[1024];   //[nVd]
   Int_t           vdPdgId[1024];   //[nVd]
   Float_t         vdMom[1024];   //[nVd]
   Float_t         vdMomX[1024];   //[nVd]
   Float_t         vdMomY[1024];   //[nVd]
   Float_t         vdMomZ[1024];   //[nVd]
   Float_t         vdPosX[1024];   //[nVd]
   Float_t         vdPosY[1024];   //[nVd]
   Float_t         vdPosZ[1024];   //[nVd]
   Float_t         vdTime[1024];   //[nVd]
   Int_t           vdGenIdx[1024];   //[nVd]

   // List of branches
   TBranch        *b_evt;   //!
   TBranch        *b_nMatch;   //!
   TBranch        *b_mTrkId;   //!
   TBranch        *b_mCluId;   //!
   TBranch        *b_mChi2;   //!
   TBranch        *b_mChi2Pos;   //!
   TBranch        *b_mChi2Time;   //!
   TBranch        *b_nTrk;   //!
   TBranch        *b_trkX;   //!
   TBranch        *b_trkY;   //!
   TBranch        *b_trkZ;   //!
   TBranch        *b_trkFFX;   //!
   TBranch        *b_trkFFY;   //!
   TBranch        *b_trkFFZ;   //!
   TBranch        *b_trkCposX;   //!
   TBranch        *b_trkCposY;   //!
   TBranch        *b_trkt;   //!
   TBranch        *b_trke;   //!
   TBranch        *b_trkpX;   //!
   TBranch        *b_trkpY;   //!
   TBranch        *b_trkpZ;   //!
   TBranch        *b_trknHit;   //!
   TBranch        *b_trkStat;   //!
   TBranch        *b_trkprob;   //!
   TBranch        *b_trkd0;   //!
   TBranch        *b_trkz0;   //!
   TBranch        *b_trkphi0;   //!
   TBranch        *b_trkomega;   //!
   TBranch        *b_trkcdip;   //!
   TBranch        *b_trkdlen;   //!
   TBranch        *b_trkCluIdx;   //!
   TBranch        *b_nCluster;   //!
   TBranch        *b_cluEnergy;   //!
   TBranch        *b_cluTime;   //!
   TBranch        *b_cluCogX;   //!
   TBranch        *b_cluCogY;   //!
   TBranch        *b_cluCogZ;   //!
   TBranch        *b_cluE1;   //!
   TBranch        *b_cluE9;   //!
   TBranch        *b_cluE25;   //!
   TBranch        *b_clu2Mom;   //!
   TBranch        *b_cluAngle;   //!
   TBranch        *b_cluConv;   //!
   TBranch        *b_cluSimIdx;   //!
   TBranch        *b_cluSimLen;   //!
   TBranch        *b_cluList;   //!
   TBranch        *b_nCluSim;   //!
   TBranch        *b_clusimId;   //!
   TBranch        *b_clusimPdgId;   //!
   TBranch        *b_clusimGenIdx;   //!
   TBranch        *b_clusimCrCode;   //!
   TBranch        *b_clusimMom;   //!
   TBranch        *b_clusimPosX;   //!
   TBranch        *b_clusimPosY;   //!
   TBranch        *b_clusimPosZ;   //!
   TBranch        *b_clusimTime;   //!
   TBranch        *b_clusimEdep;   //!
   TBranch        *b_nCry;   //!
   TBranch        *b_cryId;   //!
   TBranch        *b_crySecId;   //!
   TBranch        *b_cryPosX;   //!
   TBranch        *b_cryPosY;   //!
   TBranch        *b_cryPosZ;   //!
   TBranch        *b_cryEdep;   //!
   TBranch        *b_cryTime;   //!
   TBranch        *b_nVd;   //!
   TBranch        *b_vdId;   //!
   TBranch        *b_vdPdgId;   //!
   TBranch        *b_vdMom;   //!
   TBranch        *b_vdMomX;   //!
   TBranch        *b_vdMomY;   //!
   TBranch        *b_vdMomZ;   //!
   TBranch        *b_vdPosX;   //!
   TBranch        *b_vdPosY;   //!
   TBranch        *b_vdPosZ;   //!
   TBranch        *b_vdTime;   //!
   TBranch        *b_vdGenIdx;   //!

   AnalysisBaseMatch(TTree * /*tree*/ =0) { }
   virtual ~AnalysisBaseMatch() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree*) {};
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate() {};
   virtual void    Terminate();


   TH1F*    Book1DF(TString name, TString title, Int_t nbins, Double_t low, Double_t high, TString xaxis="", TString units="");
   TH1I*    Book1DI(TString name, TString title, Int_t nbins, Double_t low, Double_t high, TString xaxis="", TString units="");
   TH2F*    Book2DF(TString name, TString title, Int_t nbinsx, Double_t xlow, Double_t xhigh, Int_t nbinsy, Double_t ylow, Double_t yhigh, TString xaxis="", TString yaxis="");
   void     SetStyle(TH1 *histo,TString xaxis, TString units);
   void     SetStyleGr(TGraph *gr);
   double   cluDistance(int ic1, int ic2); 
   double   targetExtend(int ic, int it);
   TVector3 calcCog(int ic);
   TVector3 getPosTrk(int it, double dlen);
   
   int   vdFinder();

   TFile *fHistFile;
   TString histname;
   
   Disk* disk;

   ClassDef(AnalysisBaseMatch,0);
};

#endif

#ifdef AnalysisBaseMatch_cxx
void AnalysisBaseMatch::Init(TTree *tree)
{

   cluList = 0;
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("nMatch", &nMatch, &b_nMatch);
   fChain->SetBranchAddress("mTrkId", mTrkId, &b_mTrkId);
   fChain->SetBranchAddress("mCluId", mCluId, &b_mCluId);
   fChain->SetBranchAddress("mChi2", mChi2, &b_mChi2);
   fChain->SetBranchAddress("mChi2Pos", mChi2Pos, &b_mChi2Pos);
   fChain->SetBranchAddress("mChi2Time", mChi2Time, &b_mChi2Time);
   fChain->SetBranchAddress("nTrk", &nTrk, &b_nTrk);
   fChain->SetBranchAddress("trkX", trkX, &b_trkX);
   fChain->SetBranchAddress("trkY", trkY, &b_trkY);
   fChain->SetBranchAddress("trkZ", trkZ, &b_trkZ);
   fChain->SetBranchAddress("trkFFX", trkFFX, &b_trkFFX);
   fChain->SetBranchAddress("trkFFY", trkFFY, &b_trkFFY);
   fChain->SetBranchAddress("trkFFZ", trkFFZ, &b_trkFFZ);
   fChain->SetBranchAddress("tkrCposX", trkCposX, &b_trkCposX);
   fChain->SetBranchAddress("tkrCposY", trkCposY, &b_trkCposY);
   fChain->SetBranchAddress("trkt", trkt, &b_trkt);
   fChain->SetBranchAddress("trke", trke, &b_trke);
   fChain->SetBranchAddress("trkpX", trkpX, &b_trkpX);
   fChain->SetBranchAddress("trkpY", trkpY, &b_trkpY);
   fChain->SetBranchAddress("trkpZ", trkpZ, &b_trkpZ);
   fChain->SetBranchAddress("trknHit", trknHit, &b_trknHit);
   fChain->SetBranchAddress("trkStat", trkStat, &b_trkStat);
   fChain->SetBranchAddress("trkprob", trkprob, &b_trkprob);
   fChain->SetBranchAddress("trkd0", trkd0, &b_trkd0);
   fChain->SetBranchAddress("trkz0", trkz0, &b_trkz0);
   fChain->SetBranchAddress("trkphi0", trkphi0, &b_trkphi0);
   fChain->SetBranchAddress("trkomega", trkomega, &b_trkomega);
   fChain->SetBranchAddress("trkcdip", trkcdip, &b_trkcdip);
   fChain->SetBranchAddress("trkdlen", trkdlen, &b_trkdlen);
   fChain->SetBranchAddress("trkCluIdx", trkCluIdx, &b_trkCluIdx);
   fChain->SetBranchAddress("nCluster", &nCluster, &b_nCluster);
   fChain->SetBranchAddress("cluEnergy", cluEnergy, &b_cluEnergy);
   fChain->SetBranchAddress("cluTime", cluTime, &b_cluTime);
   fChain->SetBranchAddress("cluCogX", cluCogX, &b_cluCogX);
   fChain->SetBranchAddress("cluCogY", cluCogY, &b_cluCogY);
   fChain->SetBranchAddress("cluCogZ", cluCogZ, &b_cluCogZ);
   fChain->SetBranchAddress("cluE1", cluE1, &b_cluE1);
   fChain->SetBranchAddress("cluE9", cluE9, &b_cluE9);
   fChain->SetBranchAddress("cluE25", cluE25, &b_cluE25);
   fChain->SetBranchAddress("clu2Mom", clu2Mom, &b_clu2Mom);
   fChain->SetBranchAddress("cluAngle", cluAngle, &b_cluAngle);
   fChain->SetBranchAddress("cluConv", cluConv, &b_cluConv);
   fChain->SetBranchAddress("cluSimIdx", cluSimIdx, &b_cluSimIdx);
   fChain->SetBranchAddress("cluSimLen", cluSimLen, &b_cluSimLen);
   fChain->SetBranchAddress("cluList", &cluList, &b_cluList);
   fChain->SetBranchAddress("nCluSim", &nCluSim, &b_nCluSim);
   fChain->SetBranchAddress("clusimId", clusimId, &b_clusimId);
   fChain->SetBranchAddress("clusimPdgId", clusimPdgId, &b_clusimPdgId);
   fChain->SetBranchAddress("clusimGenIdx", clusimGenIdx, &b_clusimGenIdx);
   fChain->SetBranchAddress("clusimCrCode", clusimCrCode, &b_clusimCrCode);
   fChain->SetBranchAddress("clusimMom", clusimMom, &b_clusimMom);
   fChain->SetBranchAddress("clusimPosX", clusimPosX, &b_clusimPosX);
   fChain->SetBranchAddress("clusimPosY", clusimPosY, &b_clusimPosY);
   fChain->SetBranchAddress("clusimPosZ", clusimPosZ, &b_clusimPosZ);
   fChain->SetBranchAddress("clusimTime", clusimTime, &b_clusimTime);
   fChain->SetBranchAddress("clusimEdep", clusimEdep, &b_clusimEdep);
   fChain->SetBranchAddress("nCry", &nCry, &b_nCry);
   fChain->SetBranchAddress("cryId", cryId, &b_cryId);
   fChain->SetBranchAddress("crySecId", crySecId, &b_crySecId);
   fChain->SetBranchAddress("cryPosX", cryPosX, &b_cryPosX);
   fChain->SetBranchAddress("cryPosY", cryPosY, &b_cryPosY);
   fChain->SetBranchAddress("cryPosZ", cryPosZ, &b_cryPosZ);
   fChain->SetBranchAddress("cryEdep", cryEdep, &b_cryEdep);
   fChain->SetBranchAddress("cryTime", cryTime, &b_cryTime);
   fChain->SetBranchAddress("nVd", &nVd, &b_nVd);
   fChain->SetBranchAddress("vdId", vdId, &b_vdId);
   fChain->SetBranchAddress("vdPdgId", vdPdgId, &b_vdPdgId);
   fChain->SetBranchAddress("vdMom", vdMom, &b_vdMom);
   fChain->SetBranchAddress("vdMomX", vdMomX, &b_vdMomX);
   fChain->SetBranchAddress("vdMomY", vdMomY, &b_vdMomY);
   fChain->SetBranchAddress("vdMomZ", vdMomZ, &b_vdMomZ);
   fChain->SetBranchAddress("vdPosX", vdPosX, &b_vdPosX);
   fChain->SetBranchAddress("vdPosY", vdPosY, &b_vdPosY);
   fChain->SetBranchAddress("vdPosZ", vdPosZ, &b_vdPosZ);
   fChain->SetBranchAddress("vdTime", vdTime, &b_vdTime);
   fChain->SetBranchAddress("vdGenIdx", vdGenIdx, &b_vdGenIdx);
}

Bool_t AnalysisBaseMatch::Notify()
{
   return kTRUE;
}

#endif // #ifdef AnalysisBaseMatch_cxx
