//_____________________________________________________________________________
void draw (TH1F* hist, TPad* pad, int subpad) {
  pad->cd(subpad);
  hist->Draw();
}


//_____________________________________________________________________________
void draw1d (TH1F* hist, int b1, int b2) {
  // draw a given bin subrange of a histogram 

  TH1F* dummy  = new TH1F(*hist);
  dummy->GetXaxis()->SetRange(b1,b2);
  dummy->Draw();
}

//_____________________________________________________________________________
void flpr(const char* Name, const char* Printer = 0) {
  // print a given file using `flpr' on a given printer

  char  command[200];
  if (Printer == 0) sprintf (command,"flpr %s&",Name);
  else              sprintf (command,"flpr -q %s %s&",Printer,Name);
  system(command);
}

//_____________________________________________________________________________
void gv(char* name) {
  // run `ghostview'

  char  command[200];
  sprintf (command,"ghostview %s&",name);
  system(command);
}

//_____________________________________________________________________________
void print1d () {
  TList *l = gDirectory->GetList();
  TIter next(l);
  TObject *obj;
  cout << "----------------------------------- list of existing histograms" << endl;
  while((obj = (TObject*) next())) {
    if (obj->InheritsFrom("TH1")) {
      printf("%20s  <%s>\n",obj->GetName(),obj->GetTitle());
    }
  }
}

//-----------------------------------------------------------------------------
// reload a shared library 'SharedLib'
//-----------------------------------------------------------------------------
void reload(const char* SharedLib) {
  gSystem->Unload(SharedLib);
  gSystem->Load  (SharedLib);
}

//_____________________________________________________________________________
void save1d (const char * file) {
  // save all the histograms from the current directory into a `file'

  TDirectory* old_dir = gDirectory;

  TFile *hfile = new TFile(file,"RECREATE","");
  TList *l = old_dir->GetList();
  TIter next(l);
  TObject *obj;

  while((obj = (TObject*) next())) {
    if (obj->InheritsFrom("TH1")) obj->Write();
  }
  hfile->Write();
  delete hfile;
  old_dir->cd();

}


//_____________________________________________________________________________
void set_style(const char* Style = "") {
  gROOT->SetStyle("Plain");
}

//_____________________________________________________________________________
void stage(int run_number) {
  // call `stage' command to see all the files of a given run
  char command[100];
  sprintf(command,"stage list -f | grep %x",run_number);
  gSystem->Exec(command);
}


//_____________________________________________________________________________
void tb(int sx=350, int sy=600) {
  // run `ghostview'
  static int nb = 0;
  nb++;
  new TBrowser(Form("tb_%i",nb),Form("tb_%i",nb),sx,sy);
}

//_____________________________________________________________________________
void tohex(int i) {
  // print integer in hex format: useful for run number conversions
  printf("%x\n",i);
}

//_____________________________________________________________________________
void toi(int x) {
  // print HEX number in integer format, useful for run number conversions
  printf("%i\n",x);
}

//_____________________________________________________________________________
int zone(int  nx, int ny) {
  printf("blin\n");
  gPad->Divide(nx,ny);
  gPad->cd(1);
}
