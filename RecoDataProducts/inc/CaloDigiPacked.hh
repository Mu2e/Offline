// Original author G. Pezzullo

#ifndef RecoDataProducts_CaloDigiPacked_hh
#define RecoDataProducts_CaloDigiPacked_hh

#include <iostream>
#include <vector>


namespace mu2e {

  class CaloDigiPacked
  {

   
    public:

      CaloDigiPacked(): _caloDigiOutput() {}
      CaloDigiPacked(std::vector<int> &caloDigiOutput): _caloDigiOutput(caloDigiOutput) {}

      std::vector<int> const& output() const {return _caloDigiOutput;}
      

    private:

      std::vector<int> _caloDigiOutput;             
    

  };

 

}

#endif
