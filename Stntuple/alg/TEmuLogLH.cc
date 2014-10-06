//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
#include "Stntuple/alg/TEmuLogLH.hh"
#include "TMath.h"
// #include "murat/plot/smooth_new.hh"

ClassImp(TEmuLogLH)
//-----------------------------------------------------------------------------
// default constructor
//-----------------------------------------------------------------------------
TEmuLogLH::TEmuLogLH() {
  int const nsl = 7;
  float path[nsl+1] = { 0., 50., 100., 150., 200., 250., 300., 1.e12};

  fEleEpVsPath = NULL;
  fMuoEpVsPath = NULL;

  fEleDtHist = NULL;
  fEleXsHist = NULL;
  fEleDeHist = NULL;

  fMuoDtHist = NULL;
  fMuoXsHist = NULL;
  fMuoDeHist = NULL;

  fNEpSlices = nsl;
  for (int i=0; i<nsl+1; i++) {
    fPath[i] = path[i];
  }

  for (int i=0; i<nsl; i++) {
    fEleEpHist[i] = NULL;
    fMuoEpHist[i] = NULL;
  }
}

//-----------------------------------------------------------------------------
// real costructor
//-----------------------------------------------------------------------------
TEmuLogLH::TEmuLogLH(const char* FnEle,  const char* FnMuo) {
  //  TH1  *h1;

  printf(" >>> ERROR : TEmuLogLH::TEmuLogLH(const char* FnEle,  const char* FnMuo) called by mistake\n");
  
  // h1 = get_mu2e_histogram(FnEle,"TCalm002/Hist/trk_1/dt");
  // fEleDtHist = (TH1*) h1->Clone("ele_dt");

  // // rebin if necessary
  // fEleDtHist->Rebin(4);
  // fEleDtFunc = new smooth_new(fEleDtHist);

  // h1 = get_mu2e_histogram(FnEle,"TCalm002/Hist/trk_1/ep");
  // fEleEpHist = (TH1*) h1->Clone("ele_ep");

  // // rebin if necessary
  // fEleEpHist->Rebin(1);
  // fEleEpFunc = new smooth_new(fEleEpHist);

//-----------------------------------------------------------------------------
// muon part
//-----------------------------------------------------------------------------
  // h1 = get_mu2e_histogram(FnMuo,"TCalm002/Hist/trk_1/dt");
  // fMuoDtHist = (TH1*) h1->Clone("muo_dt");

  // // rebin if necessary
  // fMuoDtHist->Rebin(4);
  // fMuoDtFunc = new smooth_new(fMuoDtHist);

  // h1 = get_mu2e_histogram(FnMuo,"TCalm002/Hist/trk_1/ep");
  // fMuoEpHist = (TH1*) h1->Clone("muo_ep");

  // // rebin if necessary
  // fMuoEpHist->Rebin(1);
  // fMuoEpFunc = new smooth_new(fMuoEpHist);
}


//-----------------------------------------------------------------------------
TEmuLogLH::~TEmuLogLH() {
  if (fEleDtHist) delete fEleDtHist;
  if (fEleXsHist) delete fEleXsHist;
  if (fEleDeHist) delete fEleDeHist;

  if (fMuoDtHist) delete fMuoDtHist;
  if (fMuoXsHist) delete fMuoXsHist;
  if (fMuoDeHist) delete fMuoDeHist;

  for (int i=0; i<fNEpSlices; i++) {
    if (fEleEpHist[i]) delete fEleEpHist[i];
    if (fMuoEpHist[i]) delete fMuoEpHist[i];
  }
}

//-----------------------------------------------------------------------------
void  TEmuLogLH::SetEleXsHist(TH1* Hist) { 
  if (fEleXsHist) delete fEleXsHist;
  fEleXsHist         = (TH1F*) Hist->Clone(Form("%s_EleXsHist_clone",Hist->GetName())); 
  fEleXsHist->Scale(1./Hist->Integral());
}

// //-----------------------------------------------------------------------------
// void  TEmuLogLH::SetEleDeHist(TH1* Hist) { 
//   if (fEleDeHist) delete fEleDeHist;
//   fEleDeHist         = (TH1*) Hist->Clone(Form("%s_EleDeHist_clone",Hist->GetName())); 
//   fEleDeHist->Scale(1./Hist->Integral());
// }

//-----------------------------------------------------------------------------
void  TEmuLogLH::SetMuoXsHist(TH1* Hist) { 
  if (fMuoXsHist) delete fMuoXsHist;
  fMuoXsHist         = (TH1F*) Hist->Clone(Form("%s_MuoXsHist_clone",Hist->GetName())); 
  fMuoXsHist->Scale(1./Hist->Integral());
}


// //-----------------------------------------------------------------------------
// void  TEmuLogLH::SetMuoDeHist(TH1* Hist) { 
//   if (fMuoDeHist) delete fMuoDeHist;
//   fMuoDeHist         = (TH1*) Hist->Clone(Form("%s_MuoDeHist_clone",Hist->GetName())); 
//   fMuoDeHist->Scale(1./Hist->Integral());
// }


//-----------------------------------------------------------------------------
double TEmuLogLH::LogLHDt(PidData_t* Data, int PdgCode) {

  double log_lh, p1, dt;
  int    i1;

  TH1    *h1(0);   

//-----------------------------------------------------------------------------
// assume constant bin
//-----------------------------------------------------------------------------
  if      (abs(PdgCode) == 11) h1 = fEleDtHist;
  else if (abs(PdgCode) == 13) h1 = fMuoDtHist;

  dt   = Data->fDt;
  i1 = (dt-h1->GetXaxis()->GetXmin())/h1->GetBinWidth(1)+1;
  p1 = h1->GetBinContent(i1);
  if (p1 < 1.e-15) p1 = 1.e-15;

  log_lh = TMath::Log(p1);
  
  return log_lh;
}


//-----------------------------------------------------------------------------
double TEmuLogLH::LogLHEp(PidData_t* Data, int PdgCode) {

  double log_lh, ep, path, p1;
  int    isl(-1), i1;

  TH1    *h1(0);   
//-----------------------------------------------------------------------------
// first, figure out the EP slice number
//-----------------------------------------------------------------------------
  path = Data->fPath;
  ep   = Data->fEp;

  for (int i=0; i<fNEpSlices; i++) {
    if (path < fPath[i+1]) {
      isl = i;
      break;
    }
  }
//-----------------------------------------------------------------------------
// choose the histogram
//-----------------------------------------------------------------------------
  if      (abs(PdgCode) == 11) h1 = fEleEpHist[isl];
  else if (abs(PdgCode) == 13) h1 = fMuoEpHist[isl];

  i1 = (ep-h1->GetXaxis()->GetXmin())/h1->GetBinWidth(1)+1;
  p1 = h1->GetBinContent(i1);
  if (p1 < 1.e-15) p1 = 1.e-15;

  log_lh = TMath::Log(p1);
  
  return log_lh;
}


//-----------------------------------------------------------------------------
double TEmuLogLH::LogLHXs(double Xs, int PdgCode) {

  double log_lh, p1;
  int    i1;

  TH1    *h1(0);   

  if      (abs(PdgCode) == 11) h1 = fEleXsHist;
  else if (abs(PdgCode) == 13) h1 = fMuoXsHist;

  i1 = (Xs-h1->GetXaxis()->GetXmin())/h1->GetBinWidth(1)+1;
  p1 = h1->GetBinContent(i1);
  if (p1 < 1.e-15) p1 = 1.e-15;

  log_lh = TMath::Log(p1);
  
  return log_lh;
}


//-----------------------------------------------------------------------------
double TEmuLogLH::LogLHCal(PidData_t* Data, int PdgCode) {

  double llh;
  llh = LogLHDt(Data,PdgCode)+LogLHEp(Data,PdgCode);

  return llh;
}

// //-----------------------------------------------------------------------------
// double TEmuLogLH::LogLHTrk(TrkData_t* Data, int PdgCode) {

//   double llh, dt, ep, p1, p2;
//   int    i1, i2;

//   TH1    *h1(0), *h2(0);   

//   xs = Data->fXs;
//   de = Data->fDe;

//   llh = LogLHXs(xs,PdgCode)+LogLHDe(de,PdgCode);

//   return llh;
// }

//-----------------------------------------------------------------------------
// 11 and 13 - electron and muon PDG codes
//-----------------------------------------------------------------------------
double TEmuLogLH::LogLHRDt(PidData_t* Data) {
  double llhr;
  llhr = LogLHDt(Data,11)-LogLHDt(Data,13);
  return llhr;
}

//-----------------------------------------------------------------------------
double TEmuLogLH::LogLHREp(PidData_t* Data) {
  double llhr;
  llhr = LogLHEp(Data,11)-LogLHEp(Data,13);
  return llhr;
}

//-----------------------------------------------------------------------------
double TEmuLogLH::LogLHRXs(double Xs) {
  double llhr;
  llhr = LogLHXs(Xs,11)-LogLHXs(Xs,13);
  return llhr;
}

// //-----------------------------------------------------------------------------
// double TEmuLogLH::LogLHRDe(double De) {
//   double llhr;
//   llhr = LogLHDe(De,11)-LogLHDe(De,13);
//   return llhr;
// }


//-----------------------------------------------------------------------------
// 11 and 13 - electron and muon PDG codes
//-----------------------------------------------------------------------------
double TEmuLogLH::LogLHRCal(PidData_t* Data) {
  double llhr = LogLHRDt(Data)+LogLHREp(Data);
  return llhr;
}

// double TEmuLogLH::LogLHRTrk(TrkData_t* Data) {
//   double llhr = LogLHRXs(Data->fXs)+LogLHRDe(Data->fDe);
//   return llhr;
// }


//-----------------------------------------------------------------------------
// so far it is trivial
//-----------------------------------------------------------------------------
void TEmuLogLH::InitEleDtHist(const char* Fn) {
  double hint;
  ReadHistogram1D(Fn,&fEleDtHist);

  hint = fEleDtHist->Integral();
  fEleDtHist->Scale(1./hint);
}

//-----------------------------------------------------------------------------
void TEmuLogLH::InitEleDtHist(const TH1F* Hist) {
  double hint;

  if (fEleDtHist) delete fEleDtHist;

  fEleDtHist = (TH1F*) Hist->Clone();

  hint = fEleDtHist->Integral();
  fEleDtHist->Scale(1./hint);
}

//-----------------------------------------------------------------------------
void TEmuLogLH::InitMuoDtHist(const char* Fn) {
  double hint;
  ReadHistogram1D(Fn,&fMuoDtHist);

  hint = fMuoDtHist->Integral();
  fMuoDtHist->Scale(1./hint);
}

//-----------------------------------------------------------------------------
void TEmuLogLH::InitMuoDtHist(const TH1F* Hist) {
  double hint;

  if (fMuoDtHist) delete fMuoDtHist;

  fMuoDtHist = (TH1F*) Hist->Clone();

  hint = fMuoDtHist->Integral();
  fMuoDtHist->Scale(1./hint);
}

//-----------------------------------------------------------------------------
void TEmuLogLH::InitEleXsHist(const char* Fn) {
  double hint;
  ReadHistogram1D(Fn,&fEleXsHist);

  hint = fEleXsHist->Integral();
  fEleXsHist->Scale(1./hint);
}

//-----------------------------------------------------------------------------
void TEmuLogLH::InitMuoXsHist(const char* Fn) {
  double hint;
  ReadHistogram1D(Fn,&fMuoXsHist);

  hint = fMuoXsHist->Integral();
  fMuoXsHist->Scale(1./hint);
}


//-----------------------------------------------------------------------------
int TEmuLogLH::ReadHistogram1D(const char* Fn, TH1F** Hist) {
  FILE  *f;
  int    done = 0, nbx, loc(0), ix, line(0);
  char   c[1000], title[200], name[200];
  float  val, xmin, xmax;

  f = fopen(Fn,"r");
  if (f == 0) {
    Error("TEmuLogLH::ReadHistogram",Form("missing file %s\n",Fn));
    return -2;
  }

  if ((*Hist) != NULL) delete (*Hist);

  while ( ((c[0]=getc(f)) != EOF) && !done) {

					// check if it is a comment line
    if (c[0] != '#') {
      ungetc(c[0],f);

      if (line == 0) {
	fscanf(f,"title: %s" ,title);
	line++;
      }
      else if (line == 1) {
	fscanf(f,"name: %s"  ,name);
	line++;
      }
      else if (line ==2) {
	fscanf(f,"nbx,xmin,xmax: %i %f %f"  ,&nbx,&xmin,&xmax);
	*Hist = new TH1F(name,title,nbx,xmin,xmax);
	line++;
      }
      else {
	for (int i=0; i<10; i++) {
	  fscanf(f,"%f" ,&val);
	  ix = loc + 1;
	  (*Hist)->SetBinContent(ix,val);
	  loc++;
	}
	line++;
      }
    }
					// skip the rest of the line
    fgets(c,100,f);
    
  }

  fclose(f);
  return 0;
}

//-----------------------------------------------------------------------------
int TEmuLogLH::ReadHistogram2D(const char* Fn, TH2F** Hist) {
  FILE  *f;
  int    done = 0, nbx, nby, loc(0), ix, iy, line(0);
  char   c[1000], title[200], name[200];
  float  val, xmin, xmax, ymin, ymax;

  f = fopen(Fn,"r");
  if (f == 0) {
    Error("TEmuLogLH::ReadHistogram",Form("missing file %s\n",Fn));
    return -2;
  }

  if ((*Hist) != NULL) delete (*Hist);

  while ( ((c[0]=getc(f)) != EOF) && !done) {

					// check if it is a comment line
    if (c[0] != '#') {
      ungetc(c[0],f);
      if (line == 0) {
	fscanf(f,"title: %s" ,title);
	line++;
      }
      else if (line == 1) {
	fscanf(f,"name: %s"  ,name);
	line++;
      }
      else if (line ==2) {
	fscanf(f,"nbx,xmin,xmax,nby,ymin,ymax: %i %f %f %i %f %f"  ,&nbx,&xmin,&xmax,&nby,&ymin,&ymax);
	*Hist = new TH2F(name,title,nbx,xmin,xmax,nby,ymin,ymax);
	line++;
      }
      else {
					// read channel number
	for (int i=0; i<10; i++) {
	  fscanf(f,"%f" ,&val);
	  iy = loc / nbx + 1;
	  ix = loc % nbx + 1;
	  (*Hist)->SetBinContent(ix,iy,val);
	  loc++;
	}
	line++;
      }
    }
					// skip the rest of the line
    fgets(c,100,f);
  }

  fclose(f);
  return 0;
}

//-----------------------------------------------------------------------------
// assume path is from 0 to 50..
// slices step by 5 cm
//-----------------------------------------------------------------------------
void TEmuLogLH::InitEleEpHist(const char* Fn) {
  int        imin[10], imax[10], nx;
  double     hint, x1;
  TH1D       *hpx;

  ReadHistogram2D(Fn,&fEleEpVsPath);

  hpx = fEleEpVsPath->ProjectionX("hpx_EleEpVsPath");

  nx                 = hpx->GetNbinsX();

  imin[0]            = 1;		// first bin
  imax[fNEpSlices-1] = nx;

  int isl=0;
  for (int ix=1; ix<=nx; ix++) {
    x1 = hpx->GetBinLowEdge(ix);
    if (x1 >= fPath[isl+1]) {
      imax[isl  ] = ix-1;
      imin[isl+1] = ix;
      isl++;
    }
  }
					// create slices
  for (int i=0; i<fNEpSlices; i++) {
    if (fEleEpHist[i] != NULL) delete fEleEpHist[i];
    fEleEpHist[i] = (TH1F*) fEleEpVsPath->ProjectionY(Form("EleEpHist_slice_%02i",i),
						      imin[i],imax[i]); 
    hint = fEleEpHist[i]->Integral();
    fEleEpHist[i]->Scale(1./hint);
  }

  delete hpx;
}

//-----------------------------------------------------------------------------
// assume path is from 0 to 50..
// slices step by 5 cm
//-----------------------------------------------------------------------------
void TEmuLogLH::InitEleEpHist(const TH2F* Hist) {
  int        imin[10], imax[10], nx;
  double     hint, x1;

  if (fEleEpVsPath) delete fEleEpVsPath;
  fEleEpVsPath = (TH2F*) Hist->Clone();

  nx                 = fEleEpVsPath->GetNbinsX();

  imin[0]            = 1;		// first bin
  imax[fNEpSlices-1] = nx;

  int isl=0;
  for (int ix=1; ix<=nx; ix++) {
    x1 = fEleEpVsPath->GetBinLowEdge(ix);
    if (x1 >= fPath[isl+1]) {
      imax[isl  ] = ix-1;
      imin[isl+1] = ix;
      isl++;
    }
  }
					// create slices
  for (int i=0; i<fNEpSlices; i++) {
    if (fEleEpHist[i] != NULL) delete fEleEpHist[i];
    fEleEpHist[i] = (TH1F*) fEleEpVsPath->ProjectionY(Form("EleEpHist_slice_%02i",i),
						      imin[i],imax[i]); 
    hint = fEleEpHist[i]->Integral();
    fEleEpHist[i]->Scale(1./hint);
  }
}

//-----------------------------------------------------------------------------
void TEmuLogLH::InitMuoEpHist(const char* Fn) {

  ReadHistogram2D(Fn,&fMuoEpVsPath);

  int        imin[10], imax[10], nx;
  double     hint, x1;
  TH1D       *hpx;

  hpx                = fMuoEpVsPath->ProjectionX("hpx_MuoEpVsPath");
  nx                 = hpx->GetNbinsX();
  imin[0]            = 1;		// first bin
  imax[fNEpSlices-1] = nx;

  int isl=0;
  for (int ix=1; ix<=nx; ix++) {
    x1 = hpx->GetBinLowEdge(ix);
    if (x1 >= fPath[isl+1]) {
      imax[isl  ] = ix-1;
      imin[isl+1] = ix;
      isl++;
    }
  }
					// create slices
  for (int i=0; i<fNEpSlices; i++) {
    if (fMuoEpHist[i]) delete fMuoEpHist[i];
    fMuoEpHist[i] = (TH1F*) fMuoEpVsPath->ProjectionY(Form("MuoEpHist_slice_%02i",i),
						      imin[i],imax[i]); 
    hint = fMuoEpHist[i]->Integral();
    fMuoEpHist[i]->Scale(1./hint);
  }

  delete hpx;
}

//-----------------------------------------------------------------------------
void TEmuLogLH::InitMuoEpHist(const TH2F* Hist) {

  int        imin[10], imax[10], nx;
  double     hint, x1;
  TH1D       *hpx;

  if (fMuoEpVsPath) delete fMuoEpVsPath;
  fMuoEpVsPath = (TH2F*) Hist->Clone();

  hpx                = fMuoEpVsPath->ProjectionX("hpx_MuoEpVsPath");
  nx                 = hpx->GetNbinsX();

  imin[0]            = 1;		// first bin
  imax[fNEpSlices-1] = nx;

  int isl=0;
  for (int ix=1; ix<=nx; ix++) {
    x1 = hpx->GetBinLowEdge(ix);
    if (x1 >= fPath[isl+1]) {
      imax[isl  ] = ix-1;
      imin[isl+1] = ix;
      isl++;
    }
  }
					// create slices
  for (int i=0; i<fNEpSlices; i++) {
    if (fMuoEpHist[i]) delete fMuoEpHist[i];
    fMuoEpHist[i] = (TH1F*) fMuoEpVsPath->ProjectionY(Form("MuoEpHist_slice_%02i",i),
						      imin[i],imax[i]); 
    hint = fMuoEpHist[i]->Integral();
    fMuoEpHist[i]->Scale(1./hint);
  }

  delete hpx;
}

//-----------------------------------------------------------------------------
// default initialization: use electron and muon templates from e00s1412 and m00s1412
// TrackAna trk_13
//-----------------------------------------------------------------------------
int TEmuLogLH::Init_v4_2_4() {

  char f_ele_ep_vs_path[256], f_ele_dt[256], f_ele_xs[256];
  char f_muo_ep_vs_path[256], f_muo_dt[256], f_muo_xs[256];

  const char* dir = getenv("MU2E_BASE_RELEASE");

  sprintf(f_ele_ep_vs_path,"%s/ConditionsService/data/pid_ele_ep_vs_path_v4_2_4.tab",dir);
  sprintf(f_ele_dt        ,"%s/ConditionsService/data/pid_ele_dt_v4_2_4.tab",dir);
  sprintf(f_ele_xs        ,"%s/ConditionsService/data/pid_ele_xs_v4_2_4.tab",dir);
  sprintf(f_muo_ep_vs_path,"%s/ConditionsService/data/pid_muo_ep_vs_path_v4_2_4.tab",dir);
  sprintf(f_muo_dt        ,"%s/ConditionsService/data/pid_muo_dt_v4_2_4.tab",dir);
  sprintf(f_muo_xs        ,"%s/ConditionsService/data/pid_muo_xs_v4_2_4.tab",dir);

  InitEleEpHist(f_ele_ep_vs_path);
  InitEleDtHist(f_ele_dt);
  InitEleXsHist(f_ele_xs);

  InitMuoEpHist(f_muo_ep_vs_path);
  InitMuoDtHist(f_muo_dt);
  InitMuoXsHist(f_muo_xs);

  return 0;
}

//-----------------------------------------------------------------------------
// default initialization: use electron and muon templates from e00s1412 and m00s1412
// TrackAna trk_13
//-----------------------------------------------------------------------------
int TEmuLogLH::Init(const char* Version) {
  
  TString ver;
  
  ver = Version;
  ver.ToLower();
  
  if (ver == "v4_2_4") Init_v4_2_4();
  else {
    printf(" >>> ERROR in TEmuLogLH::Init: unknown version : %s, BAILING OUT\n",Version); 
  }
  
  return 0;
}
