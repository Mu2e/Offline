#ifndef scorerFTDConverter_HH
#define scorerFTDConverter_HH

// Convert fluence to effective dose for most particle species using data
// taken from ICRP publication 116
//
#include "Offline/Mu2eG4/inc/scorerFTDTable.hh"

#include <vector>
#include <string>


namespace mu2e{

  class scorerFTDConverter
  {
     public:
        scorerFTDConverter(const std::string& method = "ISO");
        ~scorerFTDConverter() = default;

        void    print();
        double  evaluate(int pdgCode, double energy);

    private:
      scorerFTDTable photon;
      scorerFTDTable electron;
      scorerFTDTable positron;
      scorerFTDTable muminus;
      scorerFTDTable muplus;
      scorerFTDTable piminus;
      scorerFTDTable piplus;
      scorerFTDTable proton;
      scorerFTDTable neutron;
      scorerFTDTable helium;
  };

}
#endif
