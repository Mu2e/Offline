// This class handles a set of E/p distributions for different values
// of projected track path in the calorimeter.  Based on Stntuple's
// TEmuLogLH by Pasha Murat.
//
// Andrei Gaponenko, 2016

#ifndef ParticleID_inc_PIDLogLEp_hh
#define ParticleID_inc_PIDLogLEp_hh

#include <string>
#include <vector>

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include "fhiclcpp/ParameterSet.h"

#include "GeneralUtilities/inc/Binning.hh"
#include "GeneralUtilities/inc/NUBinning.hh"

namespace mu2e {

  class PIDLogLEp {
    static double binValueCutoff_;
    Binning epaxis_;
    NUBinning pathaxis_;
    // Outer container is vs path, inner is E/p slice
    std::vector<std::vector<double> > vals_;

  public:
    double value(double x, double y) const;

    // the smallest value() that is returned instead of log(0) if
    // the distribution vanishes at the given x.
    static double cutoff();

    struct Config {
      fhicl::Atom<std::string> inputFile {fhicl::Name("inputFile"),
          fhicl::Comment("File with a text representation of the distribution.")
          };

      fhicl::Sequence<double> pathBinBoundaries{fhicl::Name("pathBinBoundaries"),
          fhicl::Comment("Use E/p distributions for the given slices of projected track path length.")
          };
    };

    explicit PIDLogLEp(const Config& conf);
    explicit PIDLogLEp(const fhicl::ParameterSet& pset);
  };

}

#endif/*ParticleID_inc_PIDLogLEp_hh*/
