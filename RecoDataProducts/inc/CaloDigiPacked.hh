#ifndef RecoDataProducts_CaloDigiPacked_hh
#define RecoDataProducts_CaloDigiPacked_hh

// Original author G. Pezzullo

// C++ includes
#include <iostream>
#include <vector>

// Mu2e includes

namespace mu2e {

  struct CaloDigiPacked{

  public:

    CaloDigiPacked():
      _caloDigiOutput(0) {
    }

    CaloDigiPacked(const CaloDigiPacked &caloDigi):
      _caloDigiOutput(caloDigi.output()) {
    }

    CaloDigiPacked(std::vector<int> CaloDigiPackedOutput):
      _caloDigiOutput(CaloDigiPackedOutput) {
    }

    // Accessors
    std::vector<int>      output()  const { return _caloDigiOutput; }

    // Accept compiler generated versions of d'tor, copy c'tor, assignment operator.

    // Print contents of the object.
    void print( std::ostream& ost = std::cout, bool doEndl = true ) const;

  private:

    std::vector<int>     _caloDigiOutput;             
    

  };

 

} // namespace mu2e

#endif /* RecoDataProducts_CaloDigiPacked_hh */
