#ifndef Mu2eUtilities_BinnedSpectrum_hh
#define Mu2eUtilities_BinnedSpectrum_hh

//
// Base class to allow generic access to all the classes that define
// a momentum spectrum.
//
// $Id: BinnedSpectrum.hh,v 1.4 2014/04/25 17:26:42 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/04/25 17:26:42 $
//
// Original author Kyle Knoepfel 
//                 

// C++ includes
#include <assert.h>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "Mu2eUtilities/inc/Table.hh"

namespace mu2e {

  class BinnedSpectrum {

  public:

    BinnedSpectrum() : _binWidth(0.), _nBins(0) {}

    // To make accessible to CLHEP::RandGeneral and ROOT::TGraph,
    // return pointer to first vector entries
    double const * getPDF()         const { return &(*_spectrum.second.begin() ); }
    double         getPDF(size_t i) const { return _spectrum.second.at(i); }

    double const * getAbscissa()         const { return &(*_spectrum.first .begin() ); }
    double         getAbscissa(size_t i) const { return _spectrum.first .at(i); }

    size_t         getNbins()    const { return _nBins;    }
    double         getBinWidth() const { return _binWidth; }

    void           print()       const {
      for ( size_t i(0); i < _nBins ; i++ ) {
        std::cout << "  Abscissa: " << _spectrum.first.at(i)
                  << "  PDF: "      << _spectrum.second.at(i) << std::endl;
      }
    }

    template<class Shape, typename... Args>
    void initialize(double emin, double emax, double spectRes, Args... args) {

      _binWidth = spectRes;

      std::vector<double> pdf, abscissa;

      Shape s( std::forward<Args>(args)... );
      for ( double step = emin ; step <= emax ; step += spectRes ) {
        abscissa.push_back( step );
        pdf     .push_back( s.getWeight(step) );
      }

      assert( abscissa.size() == pdf.size() );
      _nBins    = abscissa.size();

      _spectrum = std::make_pair( std::move(abscissa), std::move(pdf) );

    }

    // To load data straight from a file, use the Table helper
    void initialize(const Table<2>& inputs) {
      _nBins    = inputs.getNrows();
      _spectrum.first.reserve(_nBins);
      _spectrum.second.reserve(_nBins);

      for(const auto& row : inputs.rawTable()) {
        _spectrum.first.emplace_back(row.first);
        _spectrum.second.emplace_back(row.second.at(0));
      }

      assert(_nBins > 1);
      _binWidth = (_spectrum.first.back() - _spectrum.first.front())/(_nBins - 1);
    }

  private:

    std::pair<std::vector<double>,std::vector<double>> _spectrum;
    double _binWidth;
    size_t _nBins;

  };

} // end of namespace mu2e

#endif /* Mu2eUtilities_BinnedSpectrum_hh */
