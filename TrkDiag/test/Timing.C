#include "TH1F.h"
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>

void Timing(std::vector<std::string> names,double maxtime=0.1 ){
  std::vector<TH1F*> histos;
  TCanvas* tcan = new TCanvas("timing","timing",800,800);
  unsigned nlin = rint(ceil(sqrt(names.size())));
  std::cout << "nlin = "<< nlin << std::endl;
  tcan->Clear();
  tcan->Divide(nlin,nlin);
  unsigned ican(1);
  ifstream fs;
  for(size_t iname=0;iname<names.size();++iname){
    std::string title= names[iname]+" timing;seconds;events";
    histos.push_back(new TH1F(names[iname].c_str(),title.c_str(),100,0,maxtime));
    std::string file = names[iname]+".dat";
    std::cout << "Opening file " << file << std::endl;
    fs.open(file.c_str());
    double time;
    while((fs >> time) != 0){
    //  std::cout << "time = " << time << std::endl;
      histos[iname]->Fill(time);
    }
    tcan->cd(ican);++ican;
    histos[iname]->Draw();
    fs.close();
  }

}

