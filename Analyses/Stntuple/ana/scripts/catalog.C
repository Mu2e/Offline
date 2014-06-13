///////////////////////////////////////////////////////////////////////////////
// new tasks
///////////////////////////////////////////////////////////////////////////////
#include "Stntuple/scripts/global_vars.h"
#include "Stntuple/ana/scripts/modules.hh"

def_name catalog_001("catalog");
///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// catalog, use gSystem->Setenv("DFC_PrintLevel","2") to change the 
// print level
//-----------------------------------------------------------------------------
int catalog(int Level=2, int ReportFreq=1000) {
  stntuple::m_dfc = (TDFCModule*) g.x->AddModule("TDFCModule",0,"dfc","dfc");
  stntuple::m_dfc->SetPrintLevel(Level);
  g.x->SetNEventsToReport(ReportFreq);
}
