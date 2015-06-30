///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include <map>
#include "TClass.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TFolder.h"
#include "TKey.h"
#include "TObjArray.h"
#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TBrowser.h"
#include "TCanvas.h"
#include "Stntuple/val/THistComp.hh"
#include "Stntuple/val/TGoodFolder.hh"
#include "Stntuple/val/TBadFolder.hh"
#include "Stntuple/val/stntuple_val_functions.hh"

//-----------------------------------------------------------------------------
TDirectory* get_file(const char* Filename) {
  TFile*     file;
  //  TFolder*   folder = NULL;

  file = (TFile*) gROOT->GetListOfFiles()->FindObject(Filename);
  if (! file) {
    file = TFile::Open(Filename);
    if (! file) {
      printf("can\'t open %s\n",Filename);
      return NULL;
    }
  }
  return file;
}

//-----------------------------------------------------------------------------
TFolder* read_folder(const char* Filename, const char* FolderName) {
  // find a folder in a file, if folder doesn't exist - create it

  TFile*     file;
  TFolder*   folder = NULL;
  bool opened = false;

  file = (TFile*) gROOT->GetListOfFiles()->FindObject(Filename);
  if (! file) {
    file = TFile::Open(Filename);
    if (! file) {
      printf("can\'t open %s\n",Filename);
      return NULL;
    }
    opened = true;
  }
//-----------------------------------------------------------------------------
//  file is open, look for a folder
//-----------------------------------------------------------------------------
  char folder_name[200];

  TString    s;
  s = Filename;
  s.ReplaceAll("/","_");

  sprintf(folder_name,"%s_%s",s.Data(),FolderName);

  folder = (TFolder*) gROOT->GetRootFolder()->FindObject(folder_name);

  if (! folder) {
    folder = (TFolder*) file->Get(FolderName);
    if (!folder) {
      printf("can\'t find folder %s in %s\n",FolderName,Filename);
      if(opened) file->Close();
      return NULL;
    }
  }

  folder->SetName(folder_name);
  gROOT->GetRootFolder()->Add(folder);

  return folder;
}

//_____________________________________________________________________________
int histograms_identical(TH1* Hist1, TH1* Hist2) {
  // returns 1 if 2 histograms are identical, 0 otherwise

  // make sure the pointers are not NULL:
  if (Hist1==NULL || Hist2==NULL) { return 0; }
  // make sure the number of bins is the same:
  if (Hist1->GetNbinsX() != Hist2->GetNbinsX()) { return 0; }
  // do a quick check by looking at the number of entries:
  if (Hist1->GetEntries() != Hist2->GetEntries()) { return 0; }
  // Bins and number of entries are the same, so loop over the 
  // histograms and make sure each bin matches:
  for (int ii=0; ii<=Hist1->GetNbinsX(); ii++) {
    if (Hist1->GetBinContent(ii) != Hist2->GetBinContent(ii)) {
      return 0;
    }
  }
  return 1;
}

//_____________________________________________________________________________
double histograms_fraction(TH1* Hist1, TH1* Hist2) {
  // finds the largest integral descrepancy
  // divided by the total event count

  // make sure the pointers are not NULL:
  if (Hist1==NULL || Hist2==NULL) { return 0.0; }
  // make sure the number of bins is the same:
  if (Hist1->GetNbinsX() != Hist2->GetNbinsX()) { return 0.0; }
  // do a quick check by looking at the number of entries:
  if (Hist1->GetEntries()==0 || Hist2->GetEntries()==0 ) { return 0.0; }
  // Bins and number of entries are the same, so loop over the 
  // histograms and make sure each bin matches:
  double sum1=0;
  double sum2=0;
  double maxdiff=0;
  if(Hist1->InheritsFrom("TH2")) {
    for (int ii=0; ii<=Hist1->GetNbinsX()+1; ii++) {
      for (int jj=0; jj<=Hist1->GetNbinsY()+1; jj++) {
	sum1 += Hist1->GetBinContent(ii,jj);
	sum2 += Hist2->GetBinContent(ii,jj);
	double sumdiff = TMath::Abs(sum1-sum2);
	if(sumdiff>maxdiff) maxdiff = sumdiff;
      }
    }
  } else if(Hist1->InheritsFrom("TH1")) {
    for (int ii=0; ii<=Hist1->GetNbinsX()+1; ii++) {
      sum1 += Hist1->GetBinContent(ii);
      sum2 += Hist2->GetBinContent(ii);
      double sumdiff = TMath::Abs(sum1-sum2);
      if(sumdiff>maxdiff) maxdiff = sumdiff;
    }
  } else {
    printf("Error: Fraciton comparison not implemented for %s\n",
	   Hist1->ClassName());
  }
  //printf("%s %s\n",Hist1->GetName(),Hist1->ClassName());
  //printf("%f %f %f %f\n",sum1,sum2,Hist1->GetEntries(),Hist2->GetEntries());
  return 1.0-maxdiff/TMath::Max(sum1,sum2);

}

//_____________________________________________________________________________
double compare_histograms(TH1* Hist1, TH1* Hist2, Double_t MinProb, Int_t flag) {
  double res, failed, ok;

  ok     = 1;

  if (MinProb > 0) failed = 0;
  else             failed = MinProb-1;
  if(MinProb==-2.0) failed= 0.999999;


  if (Hist2 == NULL) return failed;

  if (Hist1->GetEntries() == 0) {
    if (Hist2->GetEntries() == 0) {
      res = ok;
    }
    else {
      res = failed;
    }
  }
  else if (Hist2->GetEntries() == 0) {
    res = failed;
  }
  else {
    if (Hist1->Integral() == 0) {
      if (Hist2->Integral() == 0) {
	res = ok;
      }
      else {
	res = failed;
      }
    }
    else if (Hist2->Integral() == 0) {
      res = failed;
    }
    else {
      if(flag==1) {
	res = histograms_fraction(Hist1,Hist2);
	//printf("fraction res = %f \n",res);
      }else {
	res = Hist1->KolmogorovTest(Hist2,"UO");
      }
      if (MinProb < 0 && res >=1.0 ) {
	if (histograms_identical(Hist1,Hist2)) res = ok;
	else                                   res = failed;
      }
    }
  }

  return res;
}

//_____________________________________________________________________________
int compare_arrays(TObjArray*  Arr1, 
		   TObjArray*  Arr2, 
		   Double_t    MinProb,
		   TObjArray*  Results,
		   Int_t       flag ) 
{
  // arrays are not allowed to have subfolders

  TObjArrayIter it1(Arr1);
  TObjArrayIter it2(Arr2);

  TObject    *o1, *o2;
  TH1        *h1, *h2; 
  TObjArray  *a1, *a2;

  THistComp* hc;

  double   prob;

  while ((o1 = it1.Next())) {
    o2 = Arr2->FindObject(o1->GetName());

    if (o1->InheritsFrom("TH1")) {
      h1   = (TH1*) o1;
      h2   = (TH1*) o2;
      prob = compare_histograms(h1,h2,MinProb,flag);

      if (prob < MinProb) {
	printf("%-30s %-20s histogram, KS(prob) = %10.5f \n",
	       o1->GetName(),o1->ClassName(),prob);
      }

      if (prob < MinProb) {
	hc = new TBadHistComp(h1,h2,prob);
      }
      else {
	hc = new TGoodHistComp(h1,h2,prob);
      }

      Results->Add(hc); 
    }
    else if (strcmp(o1->ClassName(),"TObjArray") == 0) {
//-----------------------------------------------------------------------------
// compare lists of histograms
//-----------------------------------------------------------------------------
      a1 = (TObjArray*) o1;
      a2 = (TObjArray*) o2;
      compare_arrays(a1,a2,MinProb,Results,flag);
    }
    else {
      printf("%-30s %-20s *************** EMOE *****************\n",
	     o1->GetName(),o1->ClassName());
    }
  }
  return 0;
}

//_____________________________________________________________________________
int compare_folders(TFolder*   Fol1, 
		    TFolder*   Fol2, 
		    Double_t   MinProb, 
		    TObjArray* Results,
		    Int_t      flag) {

  Fol1->Print();
  Fol1->GetListOfFolders()->Print();

  TIter it1(Fol1->GetListOfFolders());
  TObject  *o1, *o2;
  TH1      *h1, *h2; 
  TFolder  *f1, *f2;
  TObjArray *a1, *a2;
  THistComp* hc;

  double   prob;

//-----------------------------------------------------------------------------
//  crate array with the results
//-----------------------------------------------------------------------------
  TObjArray* array = new TObjArray();
  array->SetName(Fol1->GetName());
  Results->Add(array);
//-----------------------------------------------------------------------------
//  loop over the objects stored in this folder
//-----------------------------------------------------------------------------
  while ((o1 = it1.Next())) {
    o2 = Fol2->FindObject(o1->GetName());

    if (o1->InheritsFrom("TH1")) {
      h1 = (TH1*) o1;
      h2 = (TH1*) o2;
      prob = compare_histograms(h1,h2,MinProb,flag);

      if (prob < MinProb) {
	printf("%-30s %-20s histogram, KS(prob) = %10.5f \n",
	       o1->GetName(),o1->ClassName(),prob);
      }

      if (prob < MinProb) {
	printf("%s\n",Fol1->GetName());
	hc = new TBadHistComp(h1,h2,prob);
      }
      else {
	hc = new TGoodHistComp(h1,h2,prob);
      }

      array->Add(hc);
    }
    else if (strcmp(o1->ClassName(),"TFolder") == 0) {
      printf("%-30s %-20s folder \n",o1->GetName(),o1->ClassName());
      f1 = (TFolder*) o1;
      f2 = (TFolder*) o2;
//-----------------------------------------------------------------------------
// compare histogram contents of 2 folders
//-----------------------------------------------------------------------------
      if (f2 ) {
	compare_folders(f1,f2,MinProb,array);
      }
      else {
	printf("folder %s doesnt exist in the 2nd file\n",
	       o1->GetName());
      }
    }
    else if (strcmp(o1->ClassName(),"TObjArray") == 0) {
//-----------------------------------------------------------------------------
// compare lists of histograms
//-----------------------------------------------------------------------------
      a1 = (TObjArray*) o1;
      a2 = (TObjArray*) o2;
      compare_arrays(a1,a2,MinProb,array,flag);
    }
    else {
      printf("%-30s %-20s *************** EMOE *****************\n",
	     o1->GetName(),o1->ClassName());
    }
  }
  return 0;
}

//_____________________________________________________________________________
int compare_directories(TDirectory*   Dir1, 
			TDirectory*   Dir2, 
			Double_t      MinProb, 
			TObjArray*    Results,
			Int_t         flag) {

  TIter it1(Dir1->GetListOfKeys());
  TKey        *k1;
  TH1         *h1, *h2; 
  TDirectory  *d1, *d2;
  TObjArray   *a1, *a2;
  TObject     *o1, *o2;
  THistComp*   hc;

  double   prob;

//-----------------------------------------------------------------------------
//  crate array with the results
//-----------------------------------------------------------------------------
  TObjArray* array = new TObjArray();
  array->SetName(Dir1->GetName());
  Results->Add(array);
//-----------------------------------------------------------------------------
//  loop over the objects stored in this folder
//-----------------------------------------------------------------------------
  while ((k1 = (TKey*) it1.Next())) {

    o1 = Dir1->Get(k1->GetName());
    o2 = Dir2->Get(k1->GetName());

    if (o2 == 0) {
      printf("%-40s directory doesn't exist in the 2nd file\n",
	     o1->GetName());
    }
    else if (o1->InheritsFrom("TH1")) {
      h1 = (TH1*) o1;
      h2 = (TH1*) o2;
      prob = compare_histograms(h1,h2,MinProb,flag);

      if (prob < MinProb || (MinProb<0.0 && prob<1.0) ) {
	printf("%-25s %-15s histogram, KS(prob) = %12.7f \n",
	       o1->GetName(),o1->ClassName(),prob);
	hc = new TBadHistComp(h1,h2,prob);
	TString temp(Dir1->GetPath());
	temp.Remove(temp.Index(":"));
	temp.Append(" vs ");
	temp.Append(Dir2->GetPath());
	temp.Append("/");
	temp.Append(h1->GetName());
	hc->SetHistory(temp);
      }
      else {
	hc = new TGoodHistComp(h1,h2,prob);
      }

      array->Add(hc);
    }
    else if (strcmp(o1->ClassName(),"TDirectory") == 0 ||
	     strcmp(o1->ClassName(),"TDirectoryFile") == 0  ) {
      printf("%-25s %-15s directory \n",o1->GetName(),o1->ClassName());
      d1 = (TDirectory*) o1;
      d2 = (TDirectory*) o2;
//-----------------------------------------------------------------------------
// compare histogram contents of 2 folders
//-----------------------------------------------------------------------------
      compare_directories(d1,d2,MinProb,array,flag);
    }
    else if (strcmp(o1->ClassName(),"TObjArray") == 0) {
//-----------------------------------------------------------------------------
// compare lists of histograms
//-----------------------------------------------------------------------------
      a1 = (TObjArray*) o1;
      a2 = (TObjArray*) o2;
      compare_arrays(a1,a2,MinProb,array,flag);
    }
    else {
      printf("%-30s %-20s *************** EMOE *****************\n",
	     o1->GetName(),o1->ClassName());
    }
  }
  return 0;
}


//_____________________________________________________________________________
Double_t  calc_result(TObjArray* Array, Double_t MinProb, TFolder* Results) {

  TObjArrayIter it(Array);
  TObject*      o;
  TFolder*      res_folder;
//-----------------------------------------------------------------------------
// first pass
//-----------------------------------------------------------------------------
  double res = MinProb;
  double prob (-1.);

  while ((o = it.Next())) {

    if (o->InheritsFrom("THistComp")) {
      THistComp* hc = (THistComp*) o;
      prob = hc->GetKsProb();
      if ( ((prob < MinProb) || (MinProb<0.0 && prob<1.0) ) && Results) {
	Results->Add(hc);
      }
    }
    else if (strcmp(o->ClassName(),"TObjArray") == 0) {
 
      prob = calc_result((TObjArray*)o,MinProb,NULL);

      if (Results) {
	if ( (prob < MinProb) || (MinProb<0.0 && prob<1.0) ) {
	  res_folder = new TBadFolder();

	  res_folder->SetName(o->GetName());

	  TFolder* x = res_folder->AddFolder("x","x");
	  res_folder->GetListOfFolders()->Clear();
	  delete x;

	  Results->Add(res_folder);
//-----------------------------------------------------------------------------
// and do the second pass storing results
//-----------------------------------------------------------------------------
	  prob = calc_result((TObjArray*)o,MinProb,res_folder);
	}
      }
    }
//-----------------------------------------------------------------------------
// when Results=0, only need to figure the icon...
//-----------------------------------------------------------------------------
    if (prob < res) {
      res = prob;
    }
  }

  return res;
}

//_____________________________________________________________________________
void compare_stn_hist(const char* Filename1, 
		      const char* Filename2, 
		      Double_t    MinProb,
		      Int_t       flag) 
{

  TFolder*     fol1;
  TFolder*     fol2;
  TFolder*     res_folder;
  //  double       res;

  TObjArray*   results = new TObjArray(10);
  results->SetName("HistComparison");
//-----------------------------------------------------------------------------
// make sure both files are opened
//-----------------------------------------------------------------------------
  fol1 = read_folder(Filename1,"Ana");
  fol2 = read_folder(Filename2,"Ana");

  compare_folders(fol1,fol2,MinProb,results,flag);
//-----------------------------------------------------------------------------
// presentation part: display the results
//-----------------------------------------------------------------------------
  calc_result(results,MinProb,NULL);

  res_folder = gROOT->GetRootFolder()->AddFolder("STNTUPLE_RESULTS",
						 "STNTUPLE_RESULTS");

  calc_result(results,MinProb,res_folder);

  TBrowser* b = new TBrowser;

  if (b == NULL) printf("compare_stn_hist ERROR: coudn't open the TBrowser window\n");
} 


//_____________________________________________________________________________
void compare_prod_hist(const char* Filename1, 
		       const char* Filename2, 
		       Double_t    MinProb,
		       Int_t       flag ) 
{

  TDirectory*     dir1;
  TDirectory*     dir2;
  TFolder*        res_folder;
  //  double          res;

  TObjArray*   results = new TObjArray(10);
  results->SetName("HistComparison");
//-----------------------------------------------------------------------------
// make sure both files are opened
//-----------------------------------------------------------------------------
  dir1 = get_file(Filename1);
  dir2 = get_file(Filename2);

  compare_directories(dir1,dir2,MinProb,results,flag);
//-----------------------------------------------------------------------------
// presentation part: display the results
//-----------------------------------------------------------------------------
  calc_result(results,MinProb,NULL);

  res_folder = gROOT->GetRootFolder()->AddFolder("STNTUPLE_RESULTS",
						 "STNTUPLE_RESULTS");

  calc_result(results,MinProb,res_folder);

  TBrowser* b = new TBrowser;

  if (b == NULL) printf("compare_stn_hist ERROR: coudn't open the TBrowser window\n");
} 


//_____________________________________________________________________________
int compare_files(const char* Filename1, 
		  const char* Filename2, 
		  Double_t    MinProb  ,
		  Int_t       flag     ) {

  TFolder*        fol1;
  TFolder*        fol2;
  TDirectory*     dir1;
  TDirectory*     dir2;
  TFolder*        res_folder;
  int             rc(0);

  TObjArray*   results = new TObjArray(10);
  results->SetName("HistComparison");

  // try folders first
  fol1 = read_folder(Filename1,"Ana");
  fol2 = read_folder(Filename2,"Ana");
  if(fol1 && fol2 && strstr(fol1->ClassName(),"TFolder")!=0) {
    compare_folders(fol1,fol2,MinProb,results,flag);
  } else {
    dir1 = get_file(Filename1);
    dir2 = get_file(Filename2);
    compare_directories(dir1,dir2,MinProb,results,flag);
  }

//-----------------------------------------------------------------------------
// presentation part: display the results
//-----------------------------------------------------------------------------
  calc_result(results,MinProb,NULL);

  res_folder = gROOT->GetRootFolder()->AddFolder("STNTUPLE_RESULTS",
						 "STNTUPLE_RESULTS");

  calc_result(results,MinProb,res_folder);

  TBrowser* b = new TBrowser;

  if (b == NULL) {
    printf("compare_stn_hist ERROR: coudn't open the TBrowser window\n");
    rc = -1;
  }

  return rc;
} 




//_____________________________________________________________________________
int add_arrays(TObjArray*  Arr1, TObjArray*  Arr2) {
  // add contents of Arr2 to Arr1, arrays are not allowed to have subfolders

  TObjArrayIter it1(Arr1);
  TObjArrayIter it2(Arr2);

  TObject    *o1, *o2;
  TH1        *h1, *h2; 
  TObjArray  *a1, *a2;

  while ((o1 = it1.Next())) {
    o2 = Arr2->FindObject(o1->GetName());

    if (o1->InheritsFrom("TH1")) {
      h1   = (TH1*) o1;
      h2   = (TH1*) o2;

      h1->Add(h1,h2,1,1);
    }
    else if (strcmp(o1->ClassName(),"TObjArray") == 0) {
//-----------------------------------------------------------------------------
// lists
//-----------------------------------------------------------------------------
      a1 = (TObjArray*) o1;
      a2 = (TObjArray*) o2;
      add_arrays(a1,a2);
    }
    else {
      printf("%-30s %-20s ** add_arrays: unknown class, do not add**\n",
	     o1->GetName(),o1->ClassName());
    }
  }
  return 0;
}


//_____________________________________________________________________________
 int add_folders(TFolder* Fol1, TFolder* Fol2) {

  // for each histogram from Fol1 add corresponding histogram from Fol2 to it

  TIter it1(Fol1->GetListOfFolders());

  TObject  *o1, *o2;
  TH1      *h1, *h2; 
  TFolder  *f1, *f2;
  TObjArray *a1, *a2;
//-----------------------------------------------------------------------------
//  loop over the objects stored in Fol1 and for each object find its vis-a-vis
//  in Fol2
//-----------------------------------------------------------------------------
  while ((o1 = it1.Next())) {
    o2 = Fol2->FindObject(o1->GetName());

    if (o1->InheritsFrom("TH1")) {
      h1 = (TH1*) o1;
      h2 = (TH1*) o2;

      h1->Add(h1,h2,1,1);
    }
    else if (strcmp(o1->ClassName(),"TFolder") == 0) {
      printf("%-30s %-20s folder \n",o1->GetName(),o1->ClassName());
      f1 = (TFolder*) o1;
      f2 = (TFolder*) o2;
//-----------------------------------------------------------------------------
// add histogram contents of 2 folders
//-----------------------------------------------------------------------------
      add_folders(f1,f2);
    }
    else if (strcmp(o1->ClassName(),"TObjArray") == 0) {
//-----------------------------------------------------------------------------
// add lists (assuming they are lists of histograms)
//-----------------------------------------------------------------------------
      a1 = (TObjArray*) o1;
      a2 = (TObjArray*) o2;
      add_arrays(a1,a2);
    }
    else {
//-----------------------------------------------------------------------------
// ** unknown class **
//-----------------------------------------------------------------------------
      printf("%-30s %-20s ** unknown class, do not add\n",
	     o1->GetName(),o1->ClassName());
    }
  }
  return 0;
}



//_____________________________________________________________________________
void merge_stn_hist(const char* List, const char* OutputFile)
{
  // given List of STNTUPLE histogram files (i.e. "a/*.root") merges them and 
  // writes output into OutputFile

  FILE* file = gSystem->OpenPipe(Form("ls %s",List),"r");

  TFile      *output_file;
  TFolder    *fol1(NULL), *fol2(NULL);

  int first = 1;

  char fn[200];

  while ( fscanf(file,"%s",fn) > 0) {
    printf("merge_stn_hist: read  %s \n",fn); 
    fol2 = read_folder(fn,"Ana");
    if (first) {
					// make output file out of the 1st one 
      first = 0;
      fol1  = fol2;
    }
    else {
					// for each histogram add its contents
					// to the same histogram in fol1
      add_folders(fol1,fol2);
    }
  }
  
  printf("merge_stn_hist: writing %s\n",OutputFile); 
  output_file = new TFile(OutputFile,"update");
  fol1->SetName("Ana");
  fol1->Write();
  output_file->Close();
} 


//_____________________________________________________________________________
int write_directories(TDirectory* Dir1, TDirectory* Dir2) {
  // traverse subdirectories in Dir1 and for each subdirectory create 
  // the same in Dir2

  TObject     *o1 /*, *o2*/;
  //  TH1         *h1, *h2; 
  TDirectory  *d1, *d2;
  //  TObjArray   *a1, *a2;
  //  TKey        *key;
  const char  *name, *class_name; 
//-----------------------------------------------------------------------------
//  loop over the objects stored in Dir1 and for each object find its vis-a-vis
//  in Dir2
//-----------------------------------------------------------------------------
  TIter it(Dir1->GetList());

  Dir2->cd();

  while ((o1 = it.Next())) {
    name = o1->GetName();
    class_name = o1->ClassName();

    if (strcmp(class_name,"TDirectory") == 0) {

      printf("%-30s new directory \n",name);

      d1 = (TDirectory*) o1;
      d2 = Dir2->mkdir(name);
      write_directories(d1,d2);
      Dir2->cd();
    }
    else {
//-----------------------------------------------------------------------------
//  write new object in Dir2
//-----------------------------------------------------------------------------
      o1->Write();
    }
  }
  return 0;
}

//_____________________________________________________________________________
int create_directories(TDirectory* Dir1, TDirectory* Dir2) {
  // traverse subdirectories in Dir1 and for each subdirectory create 
  // the same in Dir2

  TObject     *o1, *o2;
  //  TH1         *h1, *h2; 
  TDirectory  *d1, *d2;
  //  TObjArray   *a1, *a2;
  TKey        *key;
  const char  *name, *class_name; 
//-----------------------------------------------------------------------------
//  loop over the objects stored in Dir1 and for each object find its vis-a-vis
//  in Dir2
//-----------------------------------------------------------------------------
  TIter it1(Dir1->GetListOfKeys());

  Dir2->cd();
  while ((key = (TKey*) it1.Next())) {
    name = key->GetName();
    class_name = key->GetClassName();
    if (strcmp(class_name,"TDirectory") == 0) {


      printf("%-30s new directory \n",name);

      d1 = (TDirectory*) Dir1->Get(name);
      d2 = Dir2->mkdir(name);
      create_directories(d1,d2);
      Dir2->cd();
    }
    else {
//-----------------------------------------------------------------------------
//  create new object in Dir2
//-----------------------------------------------------------------------------
      o1 = Dir1->Get(name);
      TClass* cl = gROOT->GetClass(class_name);
      o2 = (TObject*) cl->New();
      o1->Copy(*o2);
    }
  }
  return 0;
}


//_____________________________________________________________________________
int add_directories(TDirectory* Dir1, TDirectory* Dir2) {

  // Dir1 is in memory; add corresponding histogram from Dir2 to it
  // Dir2 will be the output

  TObject     *o1;
  TH1         *h1, *h2; 
  TDirectory  *d1, *d2;
  TObjArray   *a1, *a2;
  //  TKey        *key;
  const char  *name, *class_name;
//-----------------------------------------------------------------------------
//  loop over the objects stored in Fol1 and for each object find 
//  its vis-a-vis in Fol2
//-----------------------------------------------------------------------------
  Dir1->cd();
  TIter it1(Dir1->GetList());

  while ((o1 = it1.Next())) {

    name       = o1->GetName  ();
    class_name = o1->ClassName();

    if (o1->InheritsFrom("TH1")) {
      h1 = (TH1*) o1;
      h2 = (TH1*) Dir2->Get(o1->GetName());

      h1->Add(h1,h2,1,1);
    }
    else if (strcmp(class_name,"TDirectory") == 0) {
      printf("%-30s %-20s directory \n",name,class_name);
      d1 = (TDirectory*) o1;
      d2 = (TDirectory*) Dir2->Get(name);
//-----------------------------------------------------------------------------
// add histogram contents of 2 folders
//-----------------------------------------------------------------------------
      add_directories(d1,d2);
    }
    else if (strcmp(class_name,"TObjArray") == 0) {
//-----------------------------------------------------------------------------
// add lists (assuming they are lists of histograms)
//-----------------------------------------------------------------------------
      a1 = (TObjArray*) o1;
      a2 = (TObjArray*) Dir2->Get(name);
      add_arrays(a1,a2);
    }
    else {
//-----------------------------------------------------------------------------
// ** unknown class **
//-----------------------------------------------------------------------------
      printf("%-30s %-20s ** unknown class, do not add\n",name,class_name);
    }
  }
  return 0;
}

//_____________________________________________________________________________
void merge_prod_hist(const char* List, const char* OutputFile)
{
  // given List of STNTUPLE histogram files (i.e. "a/*.root") merges them 
  // and writes output into OutputFile

  FILE* file = gSystem->OpenPipe(Form("ls %s",List),"r");

  TFile       *output_file, *input_file;
  //  TDirectory  *dir1, *dir2;
  //  TDirectory  *output_dir;

  int first = 1;

  char fn[200];

  while ( fscanf(file,"%s",fn) > 0) {
//-----------------------------------------------------------------------------
// for each histogram add its contents to the same histogram in output_dir
//-----------------------------------------------------------------------------
    printf("merge_prod_hist: read  %s \n",fn);

    input_file = TFile::Open(fn);

    if (first) {
      first = 0;
//-----------------------------------------------------------------------------
// recreate directory structure
//-----------------------------------------------------------------------------
      create_directories(input_file,gROOT);
    }
    else {
      add_directories(gROOT,input_file);
    }
    input_file->Close();
  }
//-----------------------------------------------------------------------------
// now write the output
//-----------------------------------------------------------------------------
  output_file = TFile::Open(OutputFile,"recreate");

  write_directories(gROOT,output_file);

  output_file->Write();
  output_file->Close();
} 



//_____________________________________________________________________________
int write_web_page(const char* wfile, int flag) {

  TString wdir(wfile);
  int nfile = strlen(wfile);
  int nbase = strlen(gSystem->BaseName(wfile));
  wdir.Remove(nfile-nbase,nbase);

  TFolder* fol;

  fol = (TFolder*) gROOT->GetRootFolder()->FindObject("STNTUPLE_RESULTS");

  if(fol==NULL) {
    printf("Could not find directory STNTUPLE_RESULTS\n");
    return -1;
  }

  TObjArray arr;
  arr.Add(fol);

  TCanvas * ccc = new TCanvas();
  ccc->cd();

  TObject* o;
  int ind = 0;
  std::multimap<float,THistComp*> hists;

  while (arr.GetEntries()>0) {
    o = arr.Last();
    //printf("%s %s\n",o->ClassName(),o->GetName());
    if(o->InheritsFrom("TFolder")){
      arr.RemoveLast();
      TIter it(((TFolder*)o)->GetListOfFolders());
      TObject* o1;
      while ((o1=it.Next())) {
        arr.AddLast(o1);
      }
    } else if(!o->InheritsFrom("THistComp")) {
      printf("%4d  Skipping %s %s \n" ,ind,
	     o->ClassName(),o->GetName());
      ind++;
      arr.RemoveLast();
    } else {
      THistComp* c = (THistComp*)o;
      if(c->GetHist1()->InheritsFrom("TH2")) {
	printf("%4d  Skipping %s %s \n" ,ind,
	       o->ClassName(),o->GetName());
      } else {
	printf("%4d  Writing %s %s \n" ,ind,
	       o->ClassName(),o->GetName());
	if(flag==0) {
	  hists.insert(std::pair<float,THistComp*>(float(ind),c));
	} else {
	  hists.insert(std::pair<float,THistComp*>(c->GetKsProb(),c));
	}
      }
      ind++;
      arr.RemoveLast();
    }
  }


  FILE* pfile = fopen(wfile,"w");
  if(!pfile) {
    printf("Could not open web page\n");
    return 1;
  }
  fprintf(pfile,"<html>");
  fprintf(pfile,"<body>");

  if(flag==0) {
    fprintf(pfile,"<H2>Input Order, not sorted</H2>\n");
  } else {
    fprintf(pfile,"<H2>Sorted, lowest prob first</H2>\n");
  }

  // write them out, order determined by flag
  std::map<float,THistComp*>::iterator it;
  it = hists.begin();
  fprintf(pfile,"<TABLE>\n<TR><TD>Plot</TD><TD>Prob</TD></TR>\n");
  while(it!=hists.end()) {
    THistComp* c = it->second;
    TString str(c->GetHistory());
    str.ReplaceAll(" ","_");
    str.ReplaceAll("\t","_");
    fprintf(pfile,"<TR><TD><a href=\"#%s\">%s</a></TD><TD>%10.8f</TD></TR>\n",
	    str.Data(),c->GetHistory().Data(),c->GetKsProb());
    it++;
  }
  fprintf(pfile,"</TABLE>\n");

  it = hists.begin();
  while(it!=hists.end()) {
    fprintf(pfile,"<BR><BR><BR><HR>\n");
    THistComp* c = it->second;
    TString str(c->GetHistory());
    str.ReplaceAll(" ","_");
    str.ReplaceAll("\t","_");
    fprintf(pfile,"<a name=\"%s\"></a>\n",str.Data());

    str = c->GetHistory();
    str.Remove(0,str.Index("/")+1);
    str.ReplaceAll(" ","_");
    str.ReplaceAll("\t","_");
    str.ReplaceAll("/","_");

    str.Append(".gif");
    fprintf(pfile,"<H2>%s</H2>\n",c->GetHistory().Data());
    fprintf(pfile,"<img src=\"%s\"></img>\n",str.Data());
    str.Prepend(wdir.Data());
    c->Draw("e0");
    ccc->SaveAs(str.Data());
    it++;
  }
		 
  fprintf(pfile,"</body>");
  fprintf(pfile,"</html>");
  fclose(pfile);
  ccc->Close();
  delete ccc;

  return 0;
}

