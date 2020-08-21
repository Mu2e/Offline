#ifndef Mu2eUtilities_BinnedSpectrum_hh
#define Mu2eUtilities_BinnedSpectrum_hh

//
// Base class to allow generic access to all the classes that define
// a momentum spectrum.
//
//
// Original author Kyle Knoepfel 
//                 

// C++ includes
#include <assert.h>
#include <stddef.h>
#include <algorithm>
#include <array>
#include <iostream>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"

#include "Mu2eUtilities/inc/Table.hh"

namespace fhicl { class ParameterSet; }

namespace mu2e {

  class BinnedSpectrum {

  public:

    BinnedSpectrum(const fhicl::ParameterSet& psphys);
    BinnedSpectrum() : _fixMax(false), _finalBin(false), _binWidth(0.), _nBins(0) {}

    // To make accessible to CLHEP::RandGeneral and ROOT::TGraph,
    // return pointer to first vector entries
    double const * getPDF()         const { return &(*_spectrum.second.begin() ); }
    double         getPDF(size_t i) const { return _spectrum.second.at(i); }

    double const * getAbscissa()         const { return &(*_spectrum.first .begin() ); }
    double         getAbscissa(size_t i) const { return _spectrum.first .at(i); }

    size_t         getNbins()    const { return _nBins;    }
    double         getBinWidth() const { return _binWidth; }
    double         getXMaxUnbinned() const { return _xmax_unbinned;}
    double         getXMax() const { return _xmax;}
    double         getXMin() const { return _xmin;}
    double         sample(double rand) {
      double temp = _xmin + (_xmax - _xmin) * rand;
      if (_finalBin && temp > _xmax-_binWidth){
        // for CE fix final bin
        // it is truncated so spread evenly across this bin
        double last_bin = _xmax - _binWidth;
        temp = last_bin + (temp-last_bin)/_binWidth * (_binWidth - (_xmax - _xmax_unbinned));
      }
      return temp;
    }

    void           print()       const {
      for ( size_t i(0); i < _nBins ; i++ ) {
        std::cout << "  Abscissa: " << _spectrum.first.at(i)
          << "  PDF: "      << _spectrum.second.at(i) << std::endl;
      }
    }

    template<class Shape, typename... Args>
      void initialize(double XMin, double XMax, double BinWidth, Args... args) {

      if(!(BinWidth > 0.) /*catches NaNs as well*/) {
        throw cet::exception("BADCONFIG")
          <<"BinnedSpectrum::initialize(): invalid binWidth = "<< BinWidth <<"\n";
      }

      _binWidth = BinWidth;

      std::vector<double> pdf, abscissa;

      Shape s( std::forward<Args>(args)... );
      if (_fixMax){
        for (double step = XMax-_binWidth/2.; step >= XMin + _binWidth/2.; step += _binWidth){
          abscissa.push_back( step );
          pdf     .push_back( s.getWeight(step) ); 
        }
        std::reverse(abscissa.begin(),abscissa.end());
        std::reverse(pdf.begin(),pdf.end());
      }else{
        for ( double step = XMin+_binWidth/2. ; step <= XMax-_binWidth/2.; step += _binWidth ) {	 
          abscissa.push_back( step );
          pdf     .push_back( s.getWeight(step) ); 
        }
      }
      if (_finalBin && abscissa[abscissa.size()-1] + _binWidth/2. < XMax){
        double step = abscissa[abscissa.size()-1] + _binWidth;
        abscissa.push_back(step);
        pdf.push_back(s.getWeight(step));
      }else{
        _finalBin = false;
      }

      assert( abscissa.size() == pdf.size() );
      _nBins    = abscissa.size();
      _xmin = abscissa[0] - _binWidth/2.;
      _xmax = abscissa[_nBins-1] + _binWidth/2.;
      _xmax_unbinned = XMax;

      _spectrum = std::make_pair( std::move(abscissa), std::move(pdf) );
    }

    // To load data straight from a file, use the Table helper
    // assumes bins are already centered not left edge
    void initialize(const Table<2>& inputs, bool binCentered, double elo=0, double ehi=0) {
      _nBins    = inputs.getNrows();
      _spectrum.first.reserve(_nBins);
      _spectrum.second.reserve(_nBins);

      _binWidth = inputs.rawTable()[1].first - inputs.rawTable()[0].first;

      for(const auto& row : inputs.rawTable()) {
        double x = row.first;
        double y = row.second.at(0);
        if (!binCentered){
          x += _binWidth/2.;
        }
        if (elo > 0 && x-_binWidth/2. < elo)
          continue;
        if (ehi > 0 && x+_binWidth/2. > ehi)
          break;
        _spectrum.first.emplace_back(x);
        _spectrum.second.emplace_back(y);
      }
      _nBins = _spectrum.first.size();

      assert(_nBins > 1);
      _binWidth = (_spectrum.first.back() - _spectrum.first.front())/(_nBins - 1);

      _xmin = _spectrum.first.front() - _binWidth/2.;
      _xmax = _spectrum.first.back() + _binWidth/2.;
    }

    void initialize(double constant){
      // hack to get CeEndpoint working
      _binWidth = 0;
      _xmin = constant;
      _xmax = constant;
      std::vector<double> abscissa, pdf;
      abscissa.push_back(constant);
      pdf.push_back(1);
      _nBins = 1;
      _spectrum = std::make_pair( std::move(abscissa), std::move(pdf) );
    }

  private:
    bool _fixMax;  // align bins so that top edge is exactly required xmax
    bool _finalBin; // hack for CeEndpoint, makes sure there is an extra truncated bin going all the way up to XMax

    // x values for Binned Spectrum is saved as BIN CENTERS not as left edge
    std::pair<std::vector<double>,std::vector<double>> _spectrum;

    double _xmin;			// low edge of bottom bin and high edge of top bin
    double _xmax;
    double _xmax_unbinned;              // xmax before truncation to bin edges
    double _binWidth;			// assumed to be constant
    size_t _nBins;

  };

} // end of namespace mu2e

#endif /* Mu2eUtilities_BinnedSpectrum_hh */
