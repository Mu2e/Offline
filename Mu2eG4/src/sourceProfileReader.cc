#include "Offline/Mu2eG4/inc/sourceProfileReader.hh"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "CLHEP/Random/RandFlat.h"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"

#include "G4SystemOfUnits.hh"
#include <string>


namespace mu2e {

  sourceProfileReader::sourceProfileReader(const std::string& filename, CLHEP::HepRandomEngine& engine)
   : randFlat_(engine)
  {
     readTimeProfile(filename);
  }

  void sourceProfileReader::readTimeProfile(const std::string& filename)
  {
     bin_.clear();
     profile_.clear();

     ConfigFileLookupPolicy configFile;
     std::string file_name = configFile(filename);
     std::ifstream infile (file_name.c_str(), std::ios::in);
     if (!infile){
        throw cet::exception("INIT")<<"scorerDelayedRadiation::readTimeProfile "
                                    <<file_name<<" does not exist\n";
     }

     double bin(0.), flux(0.),rsum(0.);
     while (infile >> bin >> flux) {
       rsum += flux;
       bin_.push_back(bin * s);
       profile_.push_back(rsum);
     }

     if (bin_.size()<2) {
        throw cet::exception("INIT")<<"scorerDelayedRadiation::readTimeProfile "
                                    <<filename<<" has wrong format, less than 2 entries\n";
     }
  }

  double sourceProfileReader::GetSourceTime()
  {
    double rand = randFlat_.fire(0.0,1.0);

    size_t i(0);
    while (profile_[i] < rand) ++i;
    i = std::min(i,bin_.size()-1);  //protect against rounding errors

    return bin_[i] + randFlat_.fire(0.0,1.0)*(bin_[i+1]-bin_[i]);
  }


}
