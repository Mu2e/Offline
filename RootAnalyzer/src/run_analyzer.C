//
// cint (not c++) Root "script" to run Analyzer.C
//
// $Id: run_analyzer.C,v 1.2 2010/06/25 22:08:34 genser Exp $
// $Author: genser $
// $Date: 2010/06/25 22:08:34 $
//
// Original author KLG
//

void run_analyzer(char const * fname="data_03.root",
                  ULong64_t nevents=300000000, 
                  ULong64_t maxPrintEvents=2, 
                  char const * g4ModuleLabel="g4run",
                  Double_t minEnergy=0.001,
                  char const * cformat="png")
{

  // Do not do it in compiled scripts and named macros
  //  gROOT->Reset();

  // Get rid of grey background on print out.
  gROOT->SetStyle("Plain");

  gSystem->Load("libCintex"); 
  Cintex::Enable();
  gSystem->Load("libCore.so");
  gSystem->Load("libMatrix.so");
  gSystem->Load("libTree.so");
  gSystem->Load("libDataFormatsCommon.so");
  gSystem->Load("libDataFormatsCommon_map_plugin.so");
  gSystem->Load("libDataFormatsCommon_dict_plugin.so");
  gSystem->Load("libDataFormatsStdDictionaries_map_plugin.so");
  gSystem->Load("libDataFormatsStdDictionaries_dict_plugin.so");
  gSystem->Load("libDataFormatsProvenance.so");
  gSystem->Load("libDataFormatsProvenance_map_plugin.so");
  gSystem->Load("libDataFormatsProvenance_dict_plugin.so");
  gSystem->Load("libFWCoreFramework.so");
  gSystem->Load("libToyDP_dict_plugin.so");
  gSystem->Load("libToyDP_map_plugin.so");
  gSystem->Load("libToyDP.so");
  gSystem->Load("libNet.so");
  gSystem->Load("libTree.so");
  gSystem->Load("libIOPoolInput.so");

  gSystem->Load("libAnalyzer_C.so");
  gSystem->Load("libAnalyzerDict_C.so");
  Analyzer* analyzer = new Analyzer(fname,nevents,maxPrintEvents,g4ModuleLabel,minEnergy,cformat);
  analyzer->begin();
  analyzer->analyze();
  analyzer->plot();
  analyzer->printOutCanvases();
  analyzer->write();
  return;
} 

