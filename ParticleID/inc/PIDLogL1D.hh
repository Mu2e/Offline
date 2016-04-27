// This class represents a 1D distribution used in particle ID
// determination.
//
// Andrei Gaponenko, 2016

#ifndef ParticleID_inc_PIDLogL1D_hh
#define ParticleID_inc_PIDLogL1D_hh

#include <string>
#include <vector>

#include "fhiclcpp/types/Atom.h"

#include "fhiclcpp/ParameterSet.h"

#include "GeneralUtilities/inc/Binning.hh"

namespace mu2e {

  class PIDLogL1D {
    static double binValueCutoff_;
    Binning axis_;
    std::vector<double> vals_;

  public:
    double value(double x) const;

    // the smallest value() that is returned instead of log(0) if
    // the distribution vanishes at the given x.
    static double cutoff();

    struct Config {
      fhicl::Atom<std::string> inputFile {fhicl::Name("inputFile"),
          fhicl::Comment("File with a text representation of the distribution")
          };
    };

    explicit PIDLogL1D(const Config& conf);
    explicit PIDLogL1D(const fhicl::ParameterSet& pset);
  };

}

#endif/*ParticleID_inc_PIDLogL1D_hh*/
