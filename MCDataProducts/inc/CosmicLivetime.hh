#ifndef MCDataProducts_CosmicLivetime_hh
#define MCDataProducts_CosmicLivetime_hh

// Original author Stefano Roberto Soleti

// C++ includes
#include <iostream>
#include <vector>
#include <math.h>

// Mu2e includes

namespace mu2e {

  struct CosmicLivetime{

  public:

    CosmicLivetime():
      _primaries(0),
      _area(0),
      _lowE(0),
      _highE(0),
      _fluxConstant(0)  {
    }

    CosmicLivetime( unsigned int primaries,
                    float area,
                    float lowE,
                    float highE,
                    float fluxConstant ):
      _primaries(primaries),
      _area(area),
      _lowE(lowE),
      _highE(highE),
      _fluxConstant(fluxConstant) {
        const float eslope = -2.7;

        const float EiToOneMinusGamma = pow(_lowE, 1 + eslope);
        const float EfToOneMinusGamma = pow(_highE, 1 + eslope);
        // http://pdg.lbl.gov/2018/reviews/rpp2018-rev-cosmic-rays.pdf eq. 29.2
        _livetime = _primaries / (M_PI * _area * _fluxConstant * (EfToOneMinusGamma - EiToOneMinusGamma) / (1. + eslope));
    }

    // Accessors
    unsigned int primaries() const { return _primaries; }
    float area()             const { return _area; }
    float lowE()             const { return _lowE;}
    float highE()            const { return _highE; }
    float fluxConstant()     const { return _fluxConstant; }
    float liveTime()         const { return _livetime; }


    // Accept compiler generated versions of d'tor, copy c'tor, assignment operator.

    // Print contents of the object.
    void print( std::ostream& ost = std::cout, bool doEndl = true ) const;

  private:

    unsigned int _primaries;
    float _area;             // m2
    float _lowE;             // (GeV)
    float _highE;            // (GeV)
    float _fluxConstant;
    float _livetime;         // s

  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   CosmicLivetime const& lt){
    lt.print(ost,false);
    return ost;
  }

} // namespace mu2e

#endif /* MCDataProducts_CosmicLivetime_hh */
