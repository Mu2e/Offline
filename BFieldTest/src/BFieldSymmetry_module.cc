// Rob Kutschke, 2015
//
// Check for the expected symmetry of the magnetic field:
//  B_x(x,-y,z) =  +B_x(x,-y,z)
//  B_y(x,-y,z) =  -B_y(x,-y,z)
//  B_z(x,-y,z) =  +B_z(x,-y,z)
//
// Can be used to check full maps that are supposed to be symmetric
// to ensure that they really are.
//
// Can also be used the measure the size of the asymmetry
// of a non-symmetric map.
//
// The work is done in the beginRun member function.
// The magnetic field map may depend on run number so it is
// not available at c'to time or beginJob time.
//

#include "BFieldGeom/inc/BFieldManager.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "SeedService/inc/SeedService.hh"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "TH1F.h"
#include "TNtuple.h"
#include "TString.h"

#include <cmath>
#include <iostream>
#include <string>

// Helper classes and functions in an anonymous namespace
namespace {

    // Helper to book and fill histograms and ntuple for one map.
    // Histograms and ntuple will be put into a subdirectory named for the map.
    struct BFSymmetryChecker {
        BFSymmetryChecker(mu2e::BFMap const& map, art::TFileService& tfs, int nbins);

        void fill(CLHEP::Hep3Vector const& point,
                  CLHEP::Hep3Vector const& b1,
                  CLHEP::Hep3Vector const& b2);

       private:
        // The map we are studying.
        mu2e::BFMap const& map_;

        // The histograms and ntuple we will fill for this map.
        // The hD* histograms are the log10 of the difference (or sum) of
        // field components; difference or sum is chosen to make the difference
        // close to zero.
        TH1F* hx_ = nullptr;
        TH1F* hy_ = nullptr;
        TH1F* hz_ = nullptr;
        TH1F* hDBx_ = nullptr;
        TH1F* hDBy_ = nullptr;
        TH1F* hDBz_ = nullptr;

        TNtuple* nt_ = nullptr;
    };

    // The content of one ntuple entry;
    struct NtEntry {
        float x;
        float y;
        float z;
        float bx1;
        float by1;
        float bz1;
        float bx2;
        float by2;
        float bz2;

        // Used to create the ntuple object.
        static TString names() { return "x:y:z:bx1:by1:bz1:bx2:by2:bz2"; }

        NtEntry(CLHEP::Hep3Vector const& p,
                CLHEP::Hep3Vector const& b1,
                CLHEP::Hep3Vector const& b2)
            : x(p.x()),
              y(p.y()),
              z(p.z()),
              bx1(b1.x()),
              by1(b1.y()),
              bz1(b1.z()),
              bx2(b2.x()),
              by2(b2.y()),
              bz2(b2.z()) {}

        // Need for calls to TNtuple::Fill.
        float const* address() const { return &x; }
    };

    BFSymmetryChecker::BFSymmetryChecker(mu2e::BFMap const& amap, art::TFileService& tfs, int nbins)
        : map_(amap) {
        // Make ROOT subdirectory to hold histograms
        art::TFileDirectory tfdir = tfs.mkdir(map_.getKey().c_str());

        // Create the histograms in the subdirectory
        hx_ = tfdir.make<TH1F>("hx", "Generated x value;[mm]", nbins, map_.xmin(), map_.xmax());
        hy_ = tfdir.make<TH1F>("hy", "Generated y value;[mm]", nbins, map_.ymin(), map_.ymax());
        hz_ = tfdir.make<TH1F>("hz", "Generated z value;[mm]", nbins, map_.zmin(), map_.zmax());
        hDBx_ = tfdir.make<TH1F>("hDBx", "Log10( Delta Bx), B in Tesla", 80, -40., 40.);
        hDBy_ = tfdir.make<TH1F>("hDBy", "Log10( Delta By), B in Tesla", 80, -40., 40.);
        hDBz_ = tfdir.make<TH1F>("hDBz", "Log10( Delta Bz), B in Tesla", 80, -40., 40.);

        nt_ = tfdir.make<TNtuple>("nt", "BField Symmetry Check", NtEntry::names());
        ;
    }

    // Helper function: return log10 of the argument.
    // Second argument is the value to return if the argument is zero.
    double safeLog10(double d, double valueIfZero = -39.) {
        if (d == 0.) {
            return valueIfZero;
        }
        return log10(std::abs(d));
    }

    void BFSymmetryChecker::fill(CLHEP::Hep3Vector const& point,
                                 CLHEP::Hep3Vector const& b1,
                                 CLHEP::Hep3Vector const& b2) {
        hx_->Fill(point.x());
        hy_->Fill(point.y());
        hz_->Fill(point.z());

        hDBx_->Fill(safeLog10(b1.x() - b2.x()));
        hDBy_->Fill(safeLog10(b1.y() + b2.y()));
        hDBz_->Fill(safeLog10(b1.z() - b2.z()));

        NtEntry nt(point, b1, b2);

        nt_->Fill(nt.address());
    }

    // Find the named map.
    mu2e::BFMap const& getMap(std::string const& mapName) {
        mu2e::GeomHandle<mu2e::BFieldManager> bfmgr;

        // Look for the map in the inner maps.
        auto const& imaps = bfmgr->getInnerMaps();
        for (auto const& map : imaps) {
            if (map->getKey() == mapName) {
                return *map;
            }
        }

        // Look for the map in the outer maps.
        auto const& omaps = bfmgr->getInnerMaps();
        for (auto const& map : omaps) {
            if (map->getKey() == mapName) {
                return *map;
            }
        }

        throw cet::exception("GEOM")
            << "BFieldSymmery: cannot find the map named: " << mapName << "\n";
    }

}  // end anonymous namespace

namespace mu2e {

    class BFieldSymmetry : public art::EDAnalyzer {
       public:
        explicit BFieldSymmetry(const fhicl::ParameterSet& pset);

        void beginRun(const art::Run& run) override;
        void analyze(const art::Event&) override {}

       private:
        // Names of maps to check
        std::vector<std::string> mapNames_;

        // Number of test points to draw.
        int nPoints_;

        // Number of bins in x, y, and z histograms
        int nBins_;

        // Uniform flat random distribution.
        CLHEP::RandFlat flat_;

        // Return a random point, distributed uniformly within the volume of the map;
        CLHEP::Hep3Vector fire(mu2e::BFMap const& map);
    };

}  // namespace mu2e

mu2e::BFieldSymmetry::BFieldSymmetry(const fhicl::ParameterSet& pset)
    : art::EDAnalyzer(pset),
      mapNames_(pset.get<std::vector<std::string>>("mapNames")),
      nPoints_(pset.get<int>("nPoints")),
      nBins_(pset.get<int>("nBins")),
      flat_(createEngine(art::ServiceHandle<mu2e::SeedService>()->getSeed())) {}

void mu2e::BFieldSymmetry::beginRun(const art::Run& run) {
    art::ServiceHandle<art::TFileService> tfs;

    for (auto const& name : mapNames_) {
        mu2e::BFMap const& map = getMap(name);

        // Book the histograms for this map.
        BFSymmetryChecker checker(map, *tfs, nBins_);

        // Draw requested number of points and do the check.
        for (int i = 0; i < nPoints_; ++i) {
            // A random point in the space of the map
            CLHEP::Hep3Vector p1 = fire(map);

            // The corresponding point, reflected in the plane y=0;
            CLHEP::Hep3Vector p2(p1.x(), -p1.y(), p1.z());

            // BField at both points
            CLHEP::Hep3Vector b1, b2;
            bool stat1 = map.getBFieldWithStatus(p1, b1);
            if (!stat1) {
                throw cet::exception("GEOM") << "BFieldSymmery: bad status from map " << name
                                             << " for the point: " << p1 << "\n";
            }

            bool stat2 = map.getBFieldWithStatus(p2, b2);
            if (!stat2) {
                throw cet::exception("GEOM") << "BFieldSymmery: bad status from map " << name
                                             << " for the point: " << p1 << "\n";
            }

            // Fill histograms.
            checker.fill(p1, b1, b2);
        }
    }
}

// Return a random point, uniformly distributed over the volume of the map.
CLHEP::Hep3Vector mu2e::BFieldSymmetry::fire(mu2e::BFMap const& map) {
    CLHEP::Hep3Vector val(flat_.fire(map.xmin(), map.xmax()), flat_.fire(map.ymin(), map.ymax()),
                          flat_.fire(map.zmin(), map.zmax()));

    return val;
}

DEFINE_ART_MODULE(mu2e::BFieldSymmetry);
