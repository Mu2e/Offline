#include "ParticleID/inc/PIDLogL1D.hh"

#include <set>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <numeric>

#include "fhiclcpp/types/Table.h"

#include <stdexcept>

#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"

namespace mu2e {

  double PIDLogL1D::binValueCutoff_ = 1.e-15;

  PIDLogL1D::PIDLogL1D(const fhicl::ParameterSet& pset)
    : PIDLogL1D{fhicl::Table<Config>{pset, std::set<std::string>{}}()}
  {}

  double PIDLogL1D::value(double dt) const {
    double res = binValueCutoff_;
    auto i = axis_.findBin(dt);
    if(i != Binning::nobin) {
      res = std::max(binValueCutoff_, vals_[i]);
    }
    return log(res);
  }

  double PIDLogL1D::cutoff() {
    return log(binValueCutoff_);
  }

  PIDLogL1D::PIDLogL1D(const Config& conf) {
    const std::string resolvedFileName = ConfigFileLookupPolicy()(conf.inputFile());
    std::ifstream infile(resolvedFileName);
    if(!infile.is_open()) {
      throw cet::exception("BADCONFIG")
        <<"PIDLogL1D(): Can not open file "<<resolvedFileName
        <<" for reading\n";
    }

    // Follow Stntuple's TEmuLogLH::ReadHistogram1D(const char* Fn, TH1F** Hist)
    // in parsing the file
    int linecounter=0;
    std::string line;
    while(std::getline(infile, line)) {

      if(!line.empty() && line[0]=='#') {
        // skip comment lines without counting
        continue;
      }

      ++linecounter;

      if(linecounter < 3) {
        // skip histogran title and histogram name lines
        continue;
      }
      else if(linecounter == 3) {
        // retrieve binning information
        std::istringstream is(line);
        std::string tmp;
        int nbins;        double xmin, xmax;
        if(!(is>>tmp) || ("nbx,xmin,xmax:" != tmp) ) {
          throw cet::exception("BADINPUT")
            <<"PIDLogL1D(): error retrieving binning info header from file "<<resolvedFileName<<"\n";
        }
        if(!(is>>nbins>>xmin>>xmax)) {
          throw cet::exception("BADINPUT")
            <<"PIDLogL1D(): error retrieving binning info data from file "<<resolvedFileName<<"\n";
        }
        axis_ = Binning(nbins,xmin,xmax);
      }
      else {
        std::istringstream is(line);
        double val;
        while(is>>val) vals_.emplace_back(val);
      }
    }

    if(axis_.nbins() != vals_.size()) {
      throw cet::exception("BADINPUT")
        <<"PIDLogL1D() error: nbins != the number of values provided: nbins = "
        <<axis_.nbins()<<", num values = "<<vals_.size()<<"\n";
    }

    // Normalize the histogram to unity integral
    double sum = std::accumulate(vals_.begin(), vals_.end(), 0.);
    std::for_each(vals_.begin(), vals_.end(), [sum](double &v){ v/=sum; });
  }

}
