#ifndef StnAna_global_vars_h
#define StnAna_global_vars_h

///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
  struct  StnAnaGlobals_t {
    TStnAna*    x              = NULL;
    TString     Task           = "undefined";
    TString     JobName        = "";
    TString     JobNumber      =  "0000";
    Int_t       JobType        =  0;  // 0:Run, 1:EventList
    Int_t       DoMc           = -1;  // to be redefined in stnana
    TString     CalibPass      = "";
    TString     GoodRunList    = "";
    Int_t       MinRun         =  1;
    Int_t       MaxRun         = 999999999;
    TEventList* EventList      =  0;
    Int_t*      RunEventList   =  0;
    TString     L3TrigPath     = "";
    Int_t       DoLittle       =  0;
    TString     LittleFileName = "";
    TString     HistFileName   = "";
    TString     OutputFileName = "";
    Int_t       Debug          = 0 ;
    Int_t       IDMode         = 2 ;
    Int_t       NEvents        = 0 ;
    TObjArray*  ListOfTasks    = NULL;


    TChain*           chain         = NULL;
    TStnCatalog*      catalog       = NULL;
    TStnDataset*      dataset       = NULL;
    TStnGoodRunList*  gGoodRunList  = NULL;

  };

extern StnAnaGlobals_t         g;

class def_name {
 public:
  def_name() {} 
  def_name(const char* x) { 

    if (g.ListOfTasks == 0) g.ListOfTasks = new TObjArray();

    g.ListOfTasks->Add(new TObjString(x)); 
  }
};

#endif
