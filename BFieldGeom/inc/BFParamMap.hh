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
#include <ostream>
#include <string>
#include "BFieldGeom/inc/BFInterpolationStyle.hh"
#include "BFieldGeom/inc/BFMap.hh"
#include "BFieldGeom/inc/BFMapType.hh"
#include "BFieldGeom/inc/Container3D.hh"
#include "BFieldGeom/inc/fiteval_c2.h"
#include "CLHEP/Vector/ThreeVector.h"

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
            : BFMap(filename, xmin, ymin, zmin, xmax, ymax, zmax, atype, scale, warnIfOutside){};
        //_fitFunc = new FitFunctionMaker2(
        //    "/mu2e/app/users/bpollack/BTrk/BTrk_working/Offline/BFieldGeom/test/"
        //    "Mau10_800mm_long.csv");
        //};

        ~BFParamMap(){};

        virtual bool getBFieldWithStatus(const CLHEP::Hep3Vector&, CLHEP::Hep3Vector&) const;

        virtual bool isValid(const CLHEP::Hep3Vector& point) const;

        virtual void print(std::ostream& os) const;

       private:
        FitFunctionMaker2* _fitFunc;

        // Functions used internally and by the code that populates the maps.

        // method to store the neighbors
        bool fitFunction(const CLHEP::Hep3Vector&, CLHEP::Hep3Vector&) const;
    };  // namespace mu2e

}  // end namespace mu2e

#endif /* BFieldGeom_BFParamMap_hh */
