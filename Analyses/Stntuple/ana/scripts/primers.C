///////////////////////////////////////////////////////////////////////////////
void primer_005(const char* Book    = "stntuple/dev_242",
                const char* Dataset = "sqcd00",
                const char* Fileset = "sqcd00.0000") {

  TStnAna*     x;
  TStnCatalog* ctl = new TStnCatalog();
  TStnDataset* ds  = new TStnDataset();

  ctl->InitDataset(ds,Book,Dataset,Fileset,"");
  ds->Print();
  x = new TStnAna(ds);
                               // ... run analysis loop on 100 events
  int nevents = 100;

  x->PrintStat(100);
}
