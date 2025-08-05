#ifndef sourceProfileReader_h
#define sourceProfileReader_h 1
//
// Read source time profile for radiation studies
// Author BE
//
#include "CLHEP/Random/RandFlat.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"



namespace mu2e {

  class sourceProfileReader
  {
     public:
       sourceProfileReader(const std::string& filename, CLHEP::HepRandomEngine& engine);
       ~sourceProfileReader() = default;

       double  GetSourceTime();


     private:
       void readTimeProfile(const std::string& filename);

       CLHEP::RandFlat     randFlat_;
       std::vector<double> bin_{};
       std::vector<double> profile_{};

  };

}
#endif
