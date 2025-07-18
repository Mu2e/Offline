#include "Offline/Mu2eG4/inc/scorerFTDTable.hh"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "CLHEP/Units/SystemOfUnits.h"

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
                                                         <<"does NOT exist \n";
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

    if (!found){
      mf::LogWarning logm("G4");
      logm << "Method "<<method_<<" not found in scorerFTDTable for "<<filename_
           <<", switching to ISO"<<"\n";
    }

    while (getline(file, line)){
      std::stringstream linestream(line);
      double number(0);
      linestream >> number;
      energies_.push_back(number);

      // the fluence-to-dose conversion numbers in ICRP113 are given in pSv / cm2,
      // transform in Sv / mm2 to be compatible with Geant4
      double pico(1e-12);
      for (int i=1;i<=iColumn;++i) linestream >> number;
      coeffs_.push_back(number*pico*CLHEP::cm*CLHEP::cm/CLHEP::mm/CLHEP::mm);
    }

    if (coeffs_.empty()) throw cet::exception("BADINPUT")<<"scorerFTDTable: wrong formatting in file "
                                                         <<filename_<<".\n ";
  }


  void scorerFTDTable::print()
  {
    std::cout<<"Flux-to-effetive dose for file "<<filename_<<" with method "<<method_<<"\n";
    for (size_t i=0;i<energies_.size();++i) std::cout<<energies_[i]<<" "<<coeffs_[i]<<"\n";
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
