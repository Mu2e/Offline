// Andrei Gaponenko, 2012
// Major rewrite, Rob Kutschke, 2013
//
// Make scans along z and write out the value of the magnetic field at each point.
// The output is written to a file and each line in the file has the format:
//   x y z bx by bz
//
// These files can be easily read by gnuplot.
//
// One may make many such field scans in one job.
//
// The code looks for all parameter sets defined within the module paramter set.  Each of these sub
// parameter sets must have the format:
//
// name : { x : value y : value zmin : value zmax : value npoints : value }
//
// For each such parameter set, the code will create a file name "name.txt" in the current working
// directory.  The file will contain the above output for a scan along the line defined by the
// parameter set.
//
//
// Suggestions for future improvements:
//   1) For now, only scans along z are implemented but we can imagine adding scans along x or y
//      and scans along the on-axis coordinate through the bends.
//   2) Add TTree or TNtuple output
//

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"

#include "BFieldGeom/inc/BFieldManager.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Vector/ThreeVector.h"

#include <cstdio>
#include <iostream>
#include <vector>

#include "GeneralUtilities/inc/csv.hh"

namespace {

    // This struct has all of the information needed to make one scan.
    // It holds the information from one of the sub parameter sets.
    struct ScanConfig {
        ScanConfig(std::string const& aname, fhicl::ParameterSet const& modulePSet) : name(aname) {
            // Get the parameter set for this file from the module parameter set.
            fhicl::ParameterSet pset = modulePSet.get<fhicl::ParameterSet>(name);

            x = pset.get<double>("x");
            y = pset.get<double>("y");
            zmin = pset.get<double>("zmin");
            zmax = pset.get<double>("zmax");
            npoints = pset.get<int>("npoints");

            if (npoints < 1) {
                throw cet::exception("BFIELDTEST")
                    << "Bad value in sub parameter set named: " << name
                    << "npoints must be >=1 but it is: " << npoints << "\n";
            }
        }

        std::string name;
        double x;
        double y;
        double zmin;
        double zmax;
        int npoints;
    };

    struct CSVConfig {
        CSVConfig(std::string const& aname, fhicl::ParameterSet const& modulePSet) : name(aname) {
            // Get the parameter set for this file from the module parameter set.
            fhicl::ParameterSet pset = modulePSet.get<fhicl::ParameterSet>(name);

            csv_name = pset.get<std::string>("csv_name");
        }

        std::string name;
        std::string csv_name;
    };

    inline std::ostream& operator<<(std::ostream& ost, const ScanConfig& c) {
        ost << c.name << " x: " << c.x << " y: " << c.y << " zmin: " << c.zmin
            << " zmax: " << c.zmax << " npoints: " << c.npoints;
        return ost;
    }

    // Do a scan long the z axis.
    // void scanz(ScanConfig const& c, mu2e::BFieldManager const& bfmgr) {
    //     std::string outFile = c.name + ".txt";

    //     std::cout << "BFieldTest writing output file: " << outFile << std::endl;

    //     // Use c-style IO to get the format that I want.
    //     FILE* file = fopen(outFile.c_str(), "w");

    //     double dz((c.zmax - c.zmin) /
    //               (c.npoints > 1 ? c.npoints - 1 : 1 /* dz not used for 1 point*/));

    //     for (int i = 0; i < c.npoints; ++i) {
    //         const double z(c.zmin + i * dz);
    //         CLHEP::Hep3Vector field = bfmgr.getBField(CLHEP::Hep3Vector(c.x, c.y, z));

    //         std::fprintf(file, "%15.10f %15.10f %15.10f %15.10f %15.10f %15.10f\n", c.x, c.y, z,
    //                      field[0], field[1], field[2]);
    //     }
    // }

    void scanCSV(CSVConfig const& c, mu2e::BFieldManager const& bfmgr) {
        std::string outFile = c.name + ".txt";

        std::cout << "BFieldTest writing output file: " << outFile << std::endl;

        // Use c-style IO to get the format that I want.
        FILE* file = fopen(outFile.c_str(), "w");

        io::CSVReader<3> in(c.csv_name);
        in.set_header("x", "y", "z");
        double x, y, z;
        while (in.read_row(x, y, z)) {
            CLHEP::Hep3Vector field = bfmgr.getBField(CLHEP::Hep3Vector(x, y, z));

            std::fprintf(file, "%15.10f %15.10f %15.10f %15.10f %15.10f %15.10f\n", x, y, z,
                         field[0], field[1], field[2]);
        }
    };

};  // namespace

namespace mu2e {
    class BFieldTest01 : public art::EDAnalyzer {
       public:
        explicit BFieldTest01(const fhicl::ParameterSet& pset) : art::EDAnalyzer(pset) {
            std::vector<std::string> pskeys = pset.get_pset_names();

            if (pskeys.empty()) {
                throw cet::exception("BFIELDTEST") << "BFieldTest: there are no sub parameter "
                                                      "sets. Presume this is a serious error.\n";
            }

            // Read all of the sub parameter sets.
            scans_.reserve(pskeys.size());
            for (auto const& k : pskeys) {
                scans_.emplace_back(k, pset);
            }
        }

        void beginRun(const art::Run& run);
        void analyze(const art::Event&){};

       private:
        std::vector<CSVConfig> scans_;
    };

    void BFieldTest01::beginRun(const art::Run& run) {
        GeomHandle<BFieldManager> bfmgr;

        // Make the plots.
        for (auto const& c : scans_) {
            scanCSV(c, *bfmgr);
        }
    };

}  // namespace mu2e

DEFINE_ART_MODULE(mu2e::BFieldTest01);
