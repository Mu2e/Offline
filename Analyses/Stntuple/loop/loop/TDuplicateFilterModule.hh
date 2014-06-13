#ifndef TDuplicateFilterModule_hh
#define TDuplicateFilterModule_hh

#include "TObjArray.h"
#include "Stntuple/loop/TStnModule.hh"

class TStnRunRecord;

class TDuplicateFilterModule: public TStnModule {
  // everything is public, it is communism
public:

  TObjArray*      fListOfRunRecords;
  TStnRunRecord*  fCurrentRunRecord;
  int             fNDuplicateEvents;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TDuplicateFilterModule(const char* name  = "DuplicateFilter", 
			 const char* title = "DuplicateFilter");
  ~TDuplicateFilterModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// setters ... SetOutputDir leaks memory. I know.
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// other methods
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// overloaded methods of TStnModule
//-----------------------------------------------------------------------------
  int       BeginJob       ();
  int       BeginRun       ();
  int       Event          (int ientry);
  int       EndJob         ();

  ClassDef(TDuplicateFilterModule,0)
};
#endif
