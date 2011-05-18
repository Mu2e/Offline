#ifndef CDFTrajectory_Trajectory_hh
#define CDFTrajectory_Trajectory_hh
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Normal3D.h"
#include "CLHEP/Geometry/Plane3D.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Geometry/Vector3D.h"
#include "CLHEP/Vector/ThreeVector.h"
#include <iostream>


namespace mu2e {

class Line;
class Trajectory  {

public:

  // Constructor
  inline Trajectory();
  // Destructor
  virtual ~Trajectory();

  class Location {



  public:


  Location() : __s(0.0),_position(0.0,0.0,0.0),_direction(0.0,0.0,0.0)
      , _dirIsValid(false)
      {}

//      Location(double s, Point3D &position, Vector3D &direction) :
//        __s(s), _position(position), _direction(direction)
//        , _dirIsValid(true)
//        {}

    // Involves an extra copy/destructor over the one above, I think.
    Location(double s, HepGeom::Point3D<double> position, HepGeom::Vector3D<double> direction) :
      __s(s), _position(position), _direction(direction)
      , _dirIsValid(true)
      {}

    Location(double s,
             double a, double b, double c,
             double x, double y, double z) :
      __s(s), _position(a,b,c), _direction(x,y,z)
      , _dirIsValid(true)
      {}

    //     //KR: a new constructor, see comment above
    //     Location(double s, Point3D position,
    // 	     const Trajectory* t) :
    //       __s(s), _position(position), _direction(0.0,0.0,0.0),
    //       _traj(t), _dirIsComputed(false)
    //       {}

    // Replaced by CG 1998/09/11. Also see above comments
    Location(double s, double a, double b, double c) :
      __s(s), _position(a,b,c)
      , _dirIsValid(false)
      {}

    ~Location() {}

//      void setLocation(double s, Point3D &position, Vector3D &direction) {
//        __s=s;
//        _position=position;
//        _direction=direction;
//        _dirIsValid=true;
//      }

    // Involves an extra copy/destructor over the one above, I think.
    void setLocation(double s, HepGeom::Point3D<double> position, HepGeom::Vector3D<double> direction) {
      __s=s;
      _position=position;
      _direction=direction;
      _dirIsValid=true;
    }

    void setLocation(double s,
                     double a, double b, double c,
                     double x, double y, double z) {
      __s=s;
      _position.setX(a); // should avoid any destructors / constructors this way
      _position.setY(b); // The set? method should be inlined in CLHEP
      _position.setZ(c);
      _direction.setX(x);
      _direction.setY(y);
      _direction.setZ(z);
      _dirIsValid=true;
    }

    void setLocation(double s, double a, double b, double c) {
      __s=s;
      _position.setX(a); // should avoid any destructors / constructors this way
      _position.setY(b); // The set? method should be inlined in CLHEP
      _position.setZ(c);
//       _direction.setX(0);
//       _direction.setY(0);
//       _direction.setZ(0);
      _dirIsValid=false;
    }

    const double & s(void) const { return __s; }

    const HepGeom::Point3D<double> & position(void) const { return _position; }

    inline const HepGeom::Vector3D<double> & direction(void) const;

    void print(std::ostream & os=std::cerr) const;

  private:
    double __s;
    HepGeom::Point3D<double> _position;
    HepGeom::Vector3D<double> _direction;

    mutable bool _dirIsValid;

  };

  // Get position as a function of pathlength
  virtual HepGeom::Point3D<double> getPosition(double s = 0.0) const = 0;

  // Get direction as a function of pathlength
  virtual HepGeom::Vector3D<double> getDirection(double s = 0.0) const = 0;

  // Get both position and location as a function of pathlength
  virtual void getLocation(Trajectory::Location & loc, double s = 0.0) const;

  // Get the tangent line
  virtual Line getTangent(double s=0.0) const;

  // Get pathlength at specifed distance from z-axis
  virtual double getPathLengthAtRhoEquals(double rho) const;

  // Get D0 to a point
  virtual double getPathLengthTo(const HepGeom::Point3D<double> &point) const;


  // Get pathlength to another trajectory
  virtual double getPathLengthTo(const Trajectory &traj) const;

  // Get D0 to another trajectory
  virtual double getDzeroTo(const HepGeom::Point3D<double> &point) const;

  // Get D0 to another trajectory
  virtual double getDzeroTo(const Trajectory &traj) const;

  // Return a new intersection with a plane.  Deletion is NOT the
  // responsibility of the caller.
  virtual Location * newIntersectionWith(const HepGeom::Plane3D<double> &plane) const;

  // Get the second derivative vector.
  virtual HepGeom::Vector3D<double> getSecondDerivative(double s = 0.0) const = 0;


protected:

  static const double ECOARSE; // various numbers defining search precision
  static const double COARSE ;
  static const double FINE   ;
  static const double EFINE  ;
  static const double EEFINE ;
  static const double COMP_PREC; // Minimum precision for comparisons
  static const double NUM_PREC; // Minimum precision for things which *should* be zero



};

#include "CDFTrajectory/inc/Trajectory.icc"
} // end namespace mu2e

#endif /* CDFTrajectory_Trajectory_hh */

