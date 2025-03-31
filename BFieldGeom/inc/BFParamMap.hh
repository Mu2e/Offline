#ifndef BFieldGeom_BFParamMap_hh
#define BFieldGeom_BFParamMap_hh
//
// Class to hold one magnetic field map. The map is defined by a parametric function.
// Units are: space point in mm, field values in tesla.
//
//
// Original Brian Pollack, based on work by Krzysztof Genser, Rob Kutschke, Julie Managan, Bob
// Bernstein.

//#include <iosfwd>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "Offline/BFieldGeom/inc/BFInterpolationStyle.hh"
#include "Offline/BFieldGeom/inc/BFMap.hh"
#include "Offline/BFieldGeom/inc/BFMapType.hh"
#include "Offline/BFieldGeom/inc/Container3D.hh"
#include "CLHEP/Vector/ThreeVector.h"

using namespace std;

namespace mu2e {
    class BFParamMap : public BFMap {
       public:
        friend class BFieldManagerMaker;

        BFParamMap(std::string const& filename,
                   double xmin,
                   double xmax,
                   double ymin,
                   double ymax,
                   double zmin,
                   double zmax,
                   BFMapType::enum_type atype,
                   double scale,
                   bool warnIfOutside = false)
            : BFMap(filename, xmin, xmax, ymin, ymax, zmin, zmax, atype, scale, warnIfOutside){};

        bool getBFieldWithStatus(const CLHEP::Hep3Vector&, CLHEP::Hep3Vector&) const override;

        bool isValid(const CLHEP::Hep3Vector& point) const override;

        void print(std::ostream& os) const override;

       private:
        // objects used to store the fit parameters
        int _ns = 0;
        int _ms = 0;
        double _Reff =0.;
        vector<vector<double> > _As;
        vector<vector<double> > _Bs;
        vector<double> _Ds;
        vector<vector<double> > _kms;
        mutable vector<vector<double> > _iv;
        mutable vector<vector<double> > _ivp;

        // pre calculate additional constants needed for eval
        void calcConstants();

        // evaluate the fit for a given point.
        bool evalFit(const CLHEP::Hep3Vector&, CLHEP::Hep3Vector&) const;
    };  // namespace mu2e

}  // end namespace mu2e

#endif /* BFieldGeom_BFParamMap_hh */
