#ifndef BFieldGeom_BFParamMap_hh
#define BFieldGeom_BFParamMap_hh
//
// Class to hold one magnetic field map. The map is defined by a parametric function.
// Units are: space point in mm, field values in tesla.
//
// $Id: BFParamMap.hh,v 1.00 2017/11/10 bpollack Exp $
// $Author: bpollack$
// $Date: 2017/11/10 $
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
#include "BFieldGeom/inc/BFInterpolationStyle.hh"
#include "BFieldGeom/inc/BFMap.hh"
#include "BFieldGeom/inc/BFMapType.hh"
#include "BFieldGeom/inc/Container3D.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "GeneralUtilities/inc/csv.hh"

using namespace std;

namespace mu2e {
    class BFParamMap : public BFMap {
       public:
        friend class BFieldManagerMaker;

        BFParamMap(std::string filename,
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

        ~BFParamMap(){};

        virtual bool getBFieldWithStatus(const CLHEP::Hep3Vector&, CLHEP::Hep3Vector&) const;

        virtual bool isValid(const CLHEP::Hep3Vector& point) const;

        virtual void print(std::ostream& os) const;

       private:
        // objects used to store the fit parameters
        int _ns;
        int _ms;
        double _Reff;
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
