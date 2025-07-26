#ifndef scorerFTDTable_HH
#define scorerFTDTable_HH

// Store and extrapolate fluence-to-effective dose coefficients.
//
// File format taken from ICRP publication 116 supplemental material
//    https://www.icrp.org/publication.asp?id=icrp%20publication%20116
//

#include <vector>
#include <string>


namespace mu2e{

  class scorerFTDTable
  {
     public:
        scorerFTDTable(const std::string& filename, const std::string& method);
        ~scorerFTDTable() = default;

        void    initialize();
        void    print();
        double  evaluate(double energy);

    private:
       std::string          filename_;
       std::string          method_;
       std::vector<double>  energies_;
       std::vector<double>  coeffs_;
  };

}
#endif
