///////////////////////////////////////////////////////////////////////////////
// the 1st histogram is drawn with "ep" (data)
// the secong one    - as a histogram   (MC)
// assume that we have module_name/Hist/hist_name
// routines: get_hist, gh1, gh2
///////////////////////////////////////////////////////////////////////////////

#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TKey.h"
#include "TFile.h"
#include "TFolder.h"
#include "TROOT.h"

TH1F* get_hist(TString FName, TString ModuleName, TString HistName)
{

  TFile   *file(NULL);
  TFolder *fol;
  TH1F    *h; 

  char folder_name[200];

  if (FName != "root") {
    file = gROOT->GetFile(FName.Data());
    if (!file) file = new TFile(FName.Data());
    if (! file->IsOpen()) {
      printf(" get_hist: can't open file %s \n",FName.Data());
      return 0;
    }

    TString s = FName;
    s.ReplaceAll("/","_");
    sprintf(folder_name,"%s_Ana",s.Data());
  }
  else {
    sprintf(folder_name,"Ana");
  }
//-----------------------------------------------------------------------------
// check whether Ana is a folder or a directory, assume "Ana" is always there
//-----------------------------------------------------------------------------
  TKey* key = (TKey*) file->GetListOfKeys()->FindObject("Ana");
  if (key) {
    if (strcmp(key->GetClassName(),"TFolder") == 0) {
//-----------------------------------------------------------------------------
// file has Ana folder in it, preserve backward compatibility
//-----------------------------------------------------------------------------
      fol = (TFolder*) gROOT->GetRootFolder()->FindObject(folder_name);

      if (! fol) {
	fol = (TFolder*) file->Get("Ana");
	fol->SetName(folder_name);
	gROOT->GetRootFolder()->Add(fol);
      }

      TString hname = ModuleName+"/Hist/"+HistName;
      h  = (TH1F*) fol->FindObject(hname.Data());
    }
    else {
//-----------------------------------------------------------------------------
// "Ana" is a directory
//-----------------------------------------------------------------------------
      gROOT->cd();
      TString hname = "Ana/"+ModuleName+"/Hist/"+HistName;
      h = (TH1F*) file->Get(hname.Data());
    }
  }
  else {
//-----------------------------------------------------------------------------
// no "Ana"
//-----------------------------------------------------------------------------
    gROOT->cd();
    TString hname = ModuleName+"/Hist/"+HistName;
    h = (TH1F*) file->Get(hname.Data());
  }
  return h;
}

//_____________________________________________________________________________
TH1F* gh1(TString FName, TString ModuleName, TString HistName) {
  return get_hist(FName,ModuleName,HistName);
}


//_____________________________________________________________________________
TH2F* gh2(TString FName, TString ModuleName, TString HistName) {
  return (TH2F*) get_hist(FName,ModuleName,HistName);
}

