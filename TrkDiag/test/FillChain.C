#include <iostream>
#include <fstream>
#include <string>
#include "TChain.h"
void FillChain(TChain*& chain,const char* files, const char* tname="RKFDeM/trkdiag") {
  cout << "filling chain for tree " << tname << " from files in" << files << endl;
  if(chain == 0)chain = new TChain(tname);
  ifstream ifs(files);
  if(ifs.is_open()){
    string file;
    while(getline(ifs,file)){
      cout << "adding file " << file << endl;
      chain->Add(file.c_str());
    }
    ifs.close();
  } else
    cout << "File " << files << " can't be opened, aborting" << endl;
}

