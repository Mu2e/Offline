#include "Offline/Mu2eG4/inc/scorerFTDTable.hh"
#include "CLHEP/Units/SystemOfUnits.h"
#include "cetlib_except/exception.h"

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

namespace mu2e{

  scorerFTDTable::scorerFTDTable(const std::string& filename, const std::string& method) :
     filename_(filename),
     method_(method)
  {
    initialize();
  }


  void scorerFTDTable::initialize()
  {
    std::ifstream file(filename_);
    if (!file.is_open()) throw cet::exception("BADINPUT")<<"scorerFTDTable: file "<<filename_
                                                         <<"does NOT exist \n"<< std::endl;
    std::string line,item;

    getline(file,line);
    getline(file,line);
    getline(file,line);

    bool found(false);
    int iColumn(0);
    std::stringstream linestream(line);
    while (linestream >> item){
      if (!item.compare("(MeV)")) continue;
      if (!item.compare(method_)) {found=true; break;}
      ++iColumn;
    }

    //Need to put proper warning message here
    if (!found) std::cout<<"Method "<<method_<<" not found, switching to ISO"<<std::endl;

    while (getline(file, line)){
      std::stringstream linestream(line);
      double number;
      linestream >> number;
      energies_.push_back(number);

      for (int i=1;i<=iColumn;++i) linestream >> number;
      coeffs_.push_back(number*1e-12*CLHEP::cm*CLHEP::cm/CLHEP::mm/CLHEP::mm);
    }

    if (coeffs_.empty())   throw cet::exception("BADINPUT")<<"scorerFTDTable: wrong formatting in file "
                                                           <<filename_<<".\n "<< std::endl;
  }


  void scorerFTDTable::print()
  {
    std::cout<<"Flux-to-effetive dose for file "<<filename_<<" with method "<<method_<<std::endl;
    for (size_t i=0;i<energies_.size();++i) std::cout<<energies_[i]<<" "<<coeffs_[i]<<std::endl;
  }


  double scorerFTDTable::evaluate(double energy)
  {
    if (energy < energies_.front()) return coeffs_.front();
    if (energy > energies_.back() ) return coeffs_.back();

    size_t idx(0);
    while (energies_[idx+1]<energy) ++idx;
    return (energy-energies_[idx])*(coeffs_[idx+1]-coeffs_[idx])/(energies_[idx+1]-energies_[idx]) + coeffs_[idx];
  }
}
