#ifndef BFieldGeom_BFMap_hh
#define BFieldGeom_BFMap_hh
//
// Virtual class to hold one magnetic field map. The map is implementation is
// arbitrary, and can be grid-like, parametric, or other.
// Units are: space point in mm, field values in tesla.
//
//
// Original Rob Kutschke, based on work by Julie Managan and Bob Bernstein.
// Rewritten in part by Krzysztof Genser to save execution time
// Rewritten again by Brian Pollack to become pure-virtual base class for all types of BFMaps
//

#include <ostream>
#include <string>
#include "Offline/BFieldGeom/inc/BFInterpolationStyle.hh"
#include "Offline/BFieldGeom/inc/BFMapType.hh"
#include "Offline/BFieldGeom/inc/Container3D.hh"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {
    class BFMap {
       public:
        friend class BFieldManagerMaker;

        BFMap(std::string const& filename,
              double xmin,
              double xmax,
              double ymin,
              double ymax,
              double zmin,
              double zmax,
              BFMapType::enum_type atype,
              double scale,
              bool warnIfOutside = false)
            : _key(filename),
              _warnIfOutside(warnIfOutside),
              _xmin(xmin),
              _xmax(xmax),
              _ymin(ymin),
              _ymax(ymax),
              _zmin(zmin),
              _zmax(zmax),
              _type(atype),
              _scaleFactor(scale){};

        virtual ~BFMap() = default;
        BFMap( BFMap const&  ) = default;
        BFMap( BFMap&&       ) = default;
        BFMap& operator=(BFMap const&  ) = default;
        BFMap& operator=(BFMap&&       ) = default;

        // Accessors
        virtual bool getBFieldWithStatus(const CLHEP::Hep3Vector&, CLHEP::Hep3Vector&) const = 0;

        // Validity checker
        virtual bool isValid(const CLHEP::Hep3Vector& point) const = 0;

        double xmin() const { return _xmin; };
        double xmax() const { return _xmax; };
        double ymin() const { return _ymin; };
        double ymax() const { return _ymax; };
        double zmin() const { return _zmin; };
        double zmax() const { return _zmax; };

        BFMapType type() const { return _type; }

        const std::string& getKey() const { return _key; };

        virtual void print(std::ostream& os) const = 0;

       protected:
        // Filename, database key or other id information that describes
        // where this map came from.
        std::string _key;

        // If true, then print a warning message when a point is outside the region
        // in which the map is defined; else return a field with a value of (0.,0.,0.);
        // This does happen under normal operation of G4 so we should not warn by default.
        bool _warnIfOutside;

        // Min and Max values.
        double _xmin, _xmax, _ymin, _ymax, _zmin, _zmax;

        // G4BL, PM, or possible future types.
        BFMapType _type;

        // A scale factor applied overall.
        double _scaleFactor;

    };

}  // end namespace mu2e

#endif /* BFieldGeom_BFMap_hh */
