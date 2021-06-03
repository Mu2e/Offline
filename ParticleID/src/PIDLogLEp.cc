#include "ParticleID/inc/PIDLogLEp.hh"

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

  double PIDLogLEp::binValueCutoff_ = 1.e-15;

  PIDLogLEp::PIDLogLEp(const fhicl::ParameterSet& pset)
    : PIDLogLEp{fhicl::Table<Config>{pset, std::set<std::string>{}}()}
  {}

  double PIDLogLEp::value(double ep, double path) const {
    double res = binValueCutoff_;

    auto ix = epaxis_.findBin(ep);
    auto iy = pathaxis_.findBin(path);

    if((ix != Binning::nobin)&&(iy != Binning::nobin)) {
      res = std::max(binValueCutoff_, vals_[iy][ix]);
    }

    return log(res);
  }

  double PIDLogLEp::cutoff() {
    return log(binValueCutoff_);
  }

  PIDLogLEp::PIDLogLEp(const Config& conf) {
    auto pathbb = conf.pathBinBoundaries();
    pathaxis_ = NUBinning(pathbb.begin(), pathbb.end());

    const std::string resolvedFileName = ConfigFileLookupPolicy()(conf.inputFile());
    std::ifstream infile(resolvedFileName);
    if(!infile.is_open()) {
      throw cet::exception("BADCONFIG")
        <<"PIDLogLEp(): Can not open file "<<resolvedFileName
        <<" for reading\n";
    }

    // Follow Stntuple's TEmuLogLH::ReadHistogramEp(const char* Fn, TH1F** Hist)
    // in parsing the file
    int linecounter=0;
    std::vector<double> buf;
    Binning tmppathaxis;
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
        int nxbins, nybins;
        double xmin, xmax, ymin, ymax;
        if(!(is>>tmp) || ("nbx,xmin,xmax,nby,ymin,ymax:" != tmp) ) {
          throw cet::exception("BADINPUT")
            <<"PIDLogLEp(): error retrieving binning info header from file "<<resolvedFileName<<"\n";
        }
        if(!(is>>nxbins>>xmin>>xmax>>nybins>>ymin>>ymax)) {
          throw cet::exception("BADINPUT")
            <<"PIDLogLEp(): error retrieving binning info data from file "<<resolvedFileName<<"\n";
        }
        epaxis_ = Binning(nybins,ymin,ymax);
        tmppathaxis = Binning(nxbins, xmin, xmax);
      }
      else {
        std::istringstream is(line);
        double val;
        while(is>>val) {
          buf.emplace_back(val);
        }
      }
    }

    // check consistency
    if(buf.size() != epaxis_.nbins() * tmppathaxis.nbins()) {
      throw cet::exception("BADINPUT")
        <<"PIDLogLEp() error: inconsistency in the number of provided bin values: expect "
        <<epaxis_.nbins()*tmppathaxis.nbins()<<", got "<<buf.size()<<"\n";
    }

    // Rebin data from buf into path length bins
    for(NUBinning::IndexType i=0; i<pathaxis_.nbins(); ++i) {
      vals_.emplace_back(std::vector<double>(epaxis_.nbins(), 0.));
    }
    for(Binning::IndexType tmpbin=0; tmpbin < tmppathaxis.nbins(); ++tmpbin) {
      NUBinning::IndexType pathbin = pathaxis_.findBin(tmppathaxis.binCenter(tmpbin));
      for(Binning::IndexType iebin=0; iebin < epaxis_.nbins(); ++iebin) {
        auto ibuf = tmpbin + iebin * tmppathaxis.nbins();
        vals_.at(pathbin).at(iebin) += buf.at(ibuf);
      }
    }

    // Normalize slice histograms
    for(auto& ephist: vals_) {
      double sum = std::accumulate(ephist.begin(), ephist.end(), 0.);
      std::for_each(ephist.begin(), ephist.end(), [sum](double &v){ v/=sum; });
    }

  }
}
