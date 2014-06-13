#ifndef stntuple_val_functions_hh
#define stntuple_val_functions_hh

#include "Rtypes.h"
#include "TString.h"

class TH1;
class TH1F;
class TH2F;
class TCanvas;
class TFolder;

void merge_stn_hist(const char* InputList, 
		    const char* OutputFile);

void merge_prod_hist(const char* InputList, 
		     const char* OutputFile);

void compare_stn_hist(const char* Filename1, 
		      const char* Filename2,
		      Double_t    MinProb,
		      Int_t       flag=0);

void compare_prod_hist(const char* Filename1, 
		       const char* Filename2,
		       Double_t    MinProb,
		       Int_t       flag=0);

int  compare_files(const char* Filename1, 
		   const char* Filename2,
		   Double_t    MinProb,
		   Int_t       flag=0);

int    compare_folders   (TFolder*   Fol1, 
			  TFolder*   Fol2, 
			  Double_t   MinProb,
			  TObjArray* Result,
			  Int_t      flag=0);

int    compare_directories(TDirectory*   Dir1,
			  TDirectory*   Dir2, 
			  Double_t   MinProb,
			  TObjArray* Result,
			  Int_t      flag=0);

int    compare_arrays    (TObjArray* Arr1, 
			  TObjArray* Arr2, 
			  Double_t   MinProb,
			  TObjArray* Result,
			  Int_t      flag=0);

double compare_histograms(TH1* H1,TH1* H2, Double_t   MinProb, Int_t flag=0);


int histograms_identical(TH1 * Hist1, TH1 * Hist2);

double histograms_fraction(TH1 * Hist1, TH1 * Hist2);

int folder_to_dir(const char* InputFile, 
		  const char* OutputFile);

TCanvas* new_slide(const char* name, 
		   const char* slide_title, 
		   int nx=2, int ny=3,
		   int XSize = 600, int YSize = 800);

int write_web_page(const char* wfile, int flag=0);

void overlay_2h(TH1*        h1, 
		TH1*        h2, 
		TString     Opt        = "e,p",
		int         MinXBin    = 0,
		int         MaxXBin    = 0,
		const float MarkerSize = 0.5,
		int         Scale      = 1);

void overlay_2h_1f(TString     FName1, 
		   TString     ModuleName,
		   TString     Hist1Name,
		   TString     Hist2Name,
		   TString     Opt        = "ep",
		   const int   MinXBin    = 0,
		         int   MaxXBin    = 0,
		   const float MarkerSize = 0.5,
		   Int_t       Scale      = 1);

void overlay_2h_2f(TString     FName1, 
		   TString     FName2, 
		   TString     ModuleName,
		   TString     Hist1Name,
		   TString     Hist2Name  = "",
		   TString     Opt        = "e,p",
		   int         MinXBin    = 0,
		   int         MaxXBin    = 0,
		   const float MarkerSize = 0.5,
		   int         Scale      = 1);

TH1F* get_hist(TString FName, TString ModuleName, TString HistName);
TH1F* gh1     (TString FName, TString ModuleName, TString HistName);
TH2F* gh2     (TString FName, TString ModuleName, TString HistName);

TH1F* create_eff_hist(const char* FName,
		      const char* Module,
		      const char* DenomHist, 
		      const char* NumerHist,
		      const char* EffHist     = "eff",
		      const char* EffTitle    = "",
		      Int_t       MarkerStyle = 20, 
		      Float_t     MarkerSize  = 0.5);

void plot_hist(const char* fn, 
	       const char* Module, 
	       const char* Hist, 
	       const char* Opt   = "",
	       const int   MarkerStyle   = -1,
	       const int   MarkerSize    = -1,
	       const int   PrintFileName = 1) ;

#endif

