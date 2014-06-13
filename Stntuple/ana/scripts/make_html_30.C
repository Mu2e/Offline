
#if !defined (__CINT__) || defined (__MAKECINT__)
#include "THtml.h"
#include "TSystem.h"
#endif
//-----------------------------------------------------------------------------
// not functional until ROOT 3.00
//-----------------------------------------------------------------------------
class MyHtml: public THtml {

public:

  MyHtml() {};
  ~MyHtml() {};

  void CreateIndex(char **names, int n) {
    THtml::CreateIndex((const char**) names,n);
  }
};



  char* cl[] = { 
					// detectors
    "TBaseMuonStub",
    "TCalTower", 
    "TCesCluster",
    "TClcChannel", "TClcLayer", "TClcModule", 
    "TCmudHit", 
    "TCmueHit",
    "TCmxdHit", 
    "TCprCluster",
    "TCmpTdcWord",
    "TCmpTdcHeader",
    "TCspTdcWord",
    "TCspTdcHeader",
    "TPreFred", 
    "TSumet", 
    "TCaltrg", 
    "TMutrg", 
    "TBsctrg",
    "TTrktrg", 
    "TMulti",
					// ****** utility classes
    "TBitset",
    "TStnArrayI",
    "TMatrix33", "TMatrix55",
					// raw data blocks
    "TCalDataBlock", 
    "TCesDataBlock", 
    "TClcDataBlock", 
    "TCmuDataBlock", 
    "TCmpDataBlock", 
    "TCmxDataBlock",
    "TCprDataBlock",
    "TDcasDataBlock",
    "TFwdDetDataBlock",
    "TGenParticle",
    "TGenpBlock",

    "TObspBlock",
    "TObspParticle",
    "TObsvVertex",
					// 
    "TPesDataBlock",
    "TPesCorrectedDataBlock",
    "TPesCluster",
    "TPi0", 

    "TSiIsectBlock",
    "TStnAna", 
    "TStnBeamPos",
    "TStnClusterBlock", 
    "TStnDataBlock", 
    "TStnDBManager", 
    "TStnDileptonBlock",
    "TStnElectron", 
    "TStnElectronBlock", 
    "TStnHeaderBlock",
    "TStnJet", 
    "TStnJetBlock", 
    "TStnL1Info", 
    "TStnL2Info", 
    "TStnL3Info", 
    "TStnLepton", 
    "TStnLinkBlock",
    "TStnMassBlock",
    "TStnMetBlock",
    "TStnModule", 
    "TStnMuon", 
    "TStnNode", 
    "TStnMuonBlock",
    "TStnPhoton", 
    "TStnPhotonBlock", 
    "TStnSecVtxTag",
    "TStnSecVtxTagBlock",
    "TStnSiDigiCode",
    "TStnSiGeantIsect",
    "TStnSiHit",
    "TStnSiIsect",
    "TStnTag", 
    "TStnTagBlock",
    "TStnTau", 
    "TStnTauBlock", 
    "TStnTopSummaryBlock", 
    "TStnTrack", 
    "TStnTrackBlock", 
    "TStnTrackLinkBlock", 
    "TStnTrigger",
    "TStnTriggerTable",
    "TStnTriggerBlock",
    "TStnVertex", 
    "TStnVertexBlock", 
    "TStnEvent", 
					// reconstructed data blocks
					// algorithms
    "TStntuple",
					// modules
    "FillStntupleModule",
    "InitStntupleModule", 
    "StntupleFilterModule",
    "StntupleMakerModule",
    "StntupleModule", 
					// Run1 code
    //    "TStnRun1Event", 
    //    "TStnRun1InputModule", 

    "TSvtDataBlock",
    "TSvtTrack",
    "TSvXtTrack",

    "TSvxDataBlock",
    "TStnSiHit",
    "TStnSiIsect",

    "TTdcHeader",
    "TTdcModule",
    "TTdcWord",
    "TTfrd", 
    "TTl1d",   
    "TTl2d",   
    "TTl2dCluster",   
    "TTl2dClusterIso",   
    "TTl2dL1Decision",   
    "TTl2dL2Decision",   
    "TTl2dXtrpTrack",   
    "TTl2dSvtTrack",   
    "TTl3d",
    "TTsid",

    "TXftBlock",
    "TXftHit",
    "TXftPixel",
    "TXftTrack",


    0};

//_____________________________________________________________________________
void make_html(const char* html_dir = "~/www/Stntuple/html",
	       const char* opt      = "") 
{
  // opt = "index" : index only

  MyHtml html;
  html.SetSourceDir("Stntuple/data:include/Stntuple/data:Stntuple/obj:include/Stntuple/obj:Stntuple/mod:include/Stntuple/mod:Stntuple/alg:include/Stntuple/alg:~cdfsoft/dist/releases/development/include");
  html.SetOutputDir(html_dir);

  char cmd[200];

  int ncl = 0;
  for (int i=0; cl[i] != 0; i++) {
    if (strcmp(opt,"index") != 0) {
      sprintf(cmd,"rm -f %s/%s.hh %s/%s.h",html_dir,cl[i],html_dir,cl[i]);
      gSystem->Exec(cmd);
      sprintf(cmd,"rm -f %s/%s.html",html_dir,cl[i]);
      gSystem->Exec(cmd);
      sprintf(cmd,"rm -f %s/%s.ps",html_dir,cl[i]);
      gSystem->Exec(cmd);
      sprintf(cmd,"rm -f %s/src/%s.cxx.html",html_dir,cl[i]);
      gSystem->Exec(cmd);
      html.MakeClass(cl[i]);
    }
    ncl++;
  }

  html.CreateIndex(cl,ncl);

  sprintf(cmd,"mv %s/ClassIndex.html %s/ClassIndex_Stntuple.html",
 	  html_dir,html_dir);
  printf("%s...\n",cmd);
  gSystem->Exec(cmd);
  
} 

//_____________________________________________________________________________
int a(const char** x) {
  // a little utility routine which allows to calculate the number of
  // classes on the fly
  int i=0;
  while (x[i++] != 0) ;
  return i;
}

