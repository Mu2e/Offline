// Andrei Gaponenko, 2012

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"

#include "GeometryService/inc/GeomHandle.hh"
#include "BFieldGeom/inc/BFieldManager.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include <iostream>

namespace mu2e {
  class BFieldTest01 : public art::EDAnalyzer {
    double x_;
    double y_;
    double zmin_;
    double zmax_;
    int    npoints_;

  public:
    explicit BFieldTest01(const fhicl::ParameterSet& pset)
      : x_(pset.get<double>("x"))
      , y_(pset.get<double>("y"))
      , zmin_(pset.get<double>("zmin"))
      , zmax_(pset.get<double>("zmax"))
      , npoints_(pset.get<int>("npoints"))
    {
      if(npoints_ < 1) {
        throw cet::exception("GEOM")
          << "Bad value in config file: npoints = "<<npoints_<<", should be >=1";
      }
    }

    void beginRun(const art::Run& run);
    void analyze(const art::Event&) {}
  };

  void BFieldTest01::beginRun(const art::Run& run) {
    GeomHandle<BFieldManager> bfmgr;

    double dz((zmax_-zmin_)/(npoints_ > 1 ? npoints_ - 1 : 1 /* dz not used for 1 point*/));
    for(int i=0; i<npoints_; ++i) {
      const double z(zmin_ + i*dz);
      CLHEP::Hep3Vector field = bfmgr->getBField(CLHEP::Hep3Vector(x_, y_, z));

      // Mu2e standard sw spews out printouts that are impossible to
      // suppress. They all seem to go to stdout, so'll use stderr
      // here to get non-polluted output for gnuplot.
      std::cerr<<x_<<" "<<y_<<" "<<z<<" "<<field[0]<<" "<<field[1]<<" "<<field[2]<<std::endl;

    }
  }

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::BFieldTest01);
