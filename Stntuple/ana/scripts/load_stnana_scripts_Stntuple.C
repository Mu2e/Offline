///////////////////////////////////////////////////////////////////////////////
// classes in Stntuple/ana murat/ana have the same names, don't load 
// libmurat_ana.so and libStntuple_ana.so simultaneously
///////////////////////////////////////////////////////////////////////////////
#include "TInterpreter.h"
#include "modules.hh"

// TG3Pythia6*              py      = NULL;
//-----------------------------------------------------------------------------
int load_stnana_scripts_Stntuple() {
  char        macro[200];

  const char* script[] = { 
    "catalog.C",
    "debug.C",
    "lumi.C",
    "validation.C",
    0 
  };

  const char* work_dir = gSystem->Getenv("MU2E_TEST_RELEASE");

  TInterpreter* cint = gROOT->GetInterpreter();
  
  for (int i=0; script[i] != 0; i++) {
    sprintf(macro,"%s/Stntuple/ana/scripts/%s",work_dir,script[i]);
    if (! cint->IsLoaded(macro)) {
      cint->LoadMacro(macro);
    }
  }
  
  return 0;
}
