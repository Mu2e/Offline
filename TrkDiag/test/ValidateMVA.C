void ValidateMVA(char* filename="MVA.root"){
  TString curMacroPath(gROOT->GetMacroPath());
  gROOT->SetMacroPath(curMacroPath+":./:$ROOT_DIR/source/root/tmva/test/:");
  gROOT->LoadMacro("$ROOT_DIR/source/root/tmva/test/TMVAGui.C+");
  TMVAGui(filename);
}
