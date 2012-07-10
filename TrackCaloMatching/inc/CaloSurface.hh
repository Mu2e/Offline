//
//
//
// $Id: CaloSurface.hh,v 1.1 2012/07/10 00:02:19 gianipez Exp $
// $Author: gianipez $
// $Date: 2012/07/10 00:02:19 $
//
// Original author G. Pezzullo & G. Tassielli
//


#ifndef CALOSURFACE_HH
#define CALOSURFACE_HH


//#include "BaBar/include/DetectorModel/DetVolumeType.hh"
#include "BaBar/include/CLHEP/Geometry/HepPoint.h"
#include "BaBar/include/CLHEP/Geometry/Transformation.h"
#include "BaBar/include/DetectorModel/DetSurface.hh"

#include "CLHEP/Vector/ThreeVector.h"

#include <iostream>
class HepPoint;
class HepTransformation;

namespace mu2e{
class CaloSurface : public DetSurface {
public :
        //Construct from a transform.
        CaloSurface(const HepTransformation & trans,  double halfSide1 /*= 330.*/,  double halfSide2 /*= 1320.*/ , Hep3Vector norm);
        //DetSurface::DetSurface(trans);

        //CaloSurface(const HepTransformation & trans);

        //  copy constructor
        CaloSurface(const DetSurface& other);
        //Destructor
        ~CaloSurface ();

        //Set data members values for defining the surface sizes
        void SetHalfSide1(double side1);
        void SetHalfSide2(double side2);

        inline double GetHalfSide1() const ;
        inline double GetHalfSide2() const ;

        //  Copy function.  This includes the possibility of generating
        //  a family of surfaces from a given one, according to a parameter.
        //  The parameter represents the perpendicular 'size' of the new
        //  surface relative to the old.
        //
        DetSurface* copyOf(double fparam=0.0) const ;

        CaloSurface*copy();

        //  Return the direction normal to the surface closest to the
        //  input point (returned value is the closest distance to the
        //  surface).  The sign convention on the distance is that the HepPoint on the
        //  surface can be calculated as surfpoint = testpoint + dist*norm.
        //
        double normTo( const HepPoint& testpoint ,CLHEP::Hep3Vector& norm) const ;

        // same, including computation of the associated surface point
        double normalTo( const HepPoint&,CLHEP::Hep3Vector&,SurfacePoint&) const ;

        // first derivative of the distance to the surface from the given point
        // as a function of position along the given direction
        //
        double firstDeriv( const HepPoint&, const CLHEP::Hep3Vector& ) const ;



        //  Return the distance to the surface from a given point along
        //  the given direction.  The behavior is now explicit regarding whether
        //  the intersection forward, backward, or closest to the surface is desired
        //  (the old convention of forward is preserved as default).  The return
        //  value is now an enum, with success the same value as before, but with
        //  failure now specifying a variety of conditions.
        //  The sign convention of the returned distance is such that the HepPoint
        //  on the surface can be computed as
        //  surfpoint = testpoint + distance*dir.norm()
        //
        intertype distTo(const HepPoint& testpoint ,const CLHEP::Hep3Vector& dir, double& dist, intermode mode=closest) const ;

        // same, returning surface point.
        //intertype distanceTo( const HepPoint&,const CLHEP::Hep3Vector&, double&, SurfacePoint&,intermode mode=closest) const;

        // scalar curvature at a surface point.  This is the maximum of the 2-d curvature
        // for all points
        //
        double curvature( const SurfacePoint& ) const ;

        //  Special function to find the min/max perpendicular distance from the surface
        //  to a line segment defined by 2 points.  Used in DetSurfaceSet
        //
        void segmentMinMax(const HepPoint&,const HepPoint&,double&,double&) const ;

        // Functions that define the 'surface system of coordinates'
        //  space point for a point on the surface for the given surface coordinates
        HepPoint   spacePoint(const SurfacePoint&) const ;

        //  normal to the surface at a surface point for the given surface coordinates
        CLHEP::Hep3Vector normal(const SurfacePoint&) const ;

        //  Surface coordinate directions at a surface point
        CLHEP::Hep3Vector surfaceDirection(const SurfacePoint&,int) const;

        //  returns the surface coordinates for a given space point, with a condition
        //  flag if the point is not on the surface within tolerances.
        int surfacePoint(const HepPoint&, SurfacePoint&,double tol= 1.0/*[mm]*/ ) const ;

        //  Defines the coordinates as infinite or wrapped.  Wrapped coordinates must
        //  go from 0 to 2pi.
        bool wrappedCoordinate(int) const ;

        // superseed operator ==
        bool operator==( const CaloSurface& otherSurface ) const ;

        // printout
        void print(std::ostream& os) const;
        void printAll(std::ostream& os) const{}
private:

        double _HalfSide1;//along local X axis local frame
        double _HalfSide2;//along local z axis local frame
        Hep3Vector _norm;
};
}


#endif
