#include <iostream>
#include <fstream>
#include <string>
#include "TChain.h"
void FillChain(TChain* chain,const char* files, bool verbose=false){
  cout << "filling chain for tree " << chain->GetName() << " from files in " << files << endl;
  ifstream ifs(files);
  if(ifs.is_open()){
    string file;
    while(getline(ifs,file)){
      if(verbose)cout << "adding file " << file << endl;
      chain->Add(file.c_str());
    }
    ifs.close();
  } else
    cout << "File " << files << " can't be opened, aborting" << endl;
}

