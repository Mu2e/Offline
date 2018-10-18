#include "TChain.h"
#include "TFileCollection.h"
TChain* LoadChain(const char* filelist,const char* treepath) {
  TFileCollection tfc("list","",filelist);
  TChain* tc = new TChain(treepath);
  tc->AddFileInfoList(tfc.GetList());
  return tc;
}

