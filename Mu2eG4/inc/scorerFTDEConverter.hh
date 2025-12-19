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
      scorerFTDTable photon_;
      scorerFTDTable electron_;
      scorerFTDTable positron_;
      scorerFTDTable muminus_;
      scorerFTDTable muplus_;
      scorerFTDTable piminus_;
      scorerFTDTable piplus_;
      scorerFTDTable proton_;
      scorerFTDTable neutron_;
  };

}
#endif




#ifndef scorerFTDConverter_HH
#define scorerFTDConverter_HH


namespace mu2e
{
  class scorerFTDConverter
  {
    public:
      virtual ~scorerFTDConverter() {};
      virtual void  evaluate(int pdgCode, double energy)= 0;
      virtual void  print() = 0;
  };
}

#endif
