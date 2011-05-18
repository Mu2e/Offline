#ifndef CDFTrajectory_Line_hh
#define CDFTrajectory_Line_hh
#include "CDFTrajectory/inc/Trajectory.hh"
#include "CDFTrajectory/inc/Helix.hh"
namespace mu2e {
class Line : public Trajectory {

public:

  // Default constructor
  Line();

  // Copy constructor
  Line(const Line &right);

  // Construct from a point and a direction
  Line(const HepGeom::Point3D<double> &point, const HepGeom::Vector3D<double> &direction, bool checkNormalisation = true);

  // Destructor
  ~Line();

  // Assignment
  const Line & operator=(const Line &right);

  // Position
  virtual HepGeom::Point3D<double> getPosition(double s = 0.0) const;

  // Direction
  virtual HepGeom::Vector3D<double> getDirection(double s = 0.0) const;


  // Parameters.
  double getD0() const;
  double getZ0() const;
  double getPhi0() const;
  double getCotTheta() const;

  // Position and direction together
  virtual void getLocation(Trajectory::Location & loc, double s = 0.0) const;

  // Get pathlength at specifed distance from z-axis
  virtual double getPathLengthAtRhoEquals(double rho) const;

  // Second deiviative
  virtual HepGeom::Vector3D<double>  getSecondDerivative(double s = 0.0) const;

  //  virtual double getPathLengthTo(const Trajectory &traj)

  // This is to keep the base class methods from being hidden
  virtual double getPathLengthTo(const Trajectory &t) const {
                            return Trajectory::getPathLengthTo(t);}
  virtual double getDzeroTo     (const Trajectory &t) const {
                              return Trajectory::getDzeroTo(t);}
// Pathlength to point
  virtual double getPathLengthTo(const HepGeom::Point3D<double> &point) const;

  // DZero to point
  virtual double getDzeroTo(const HepGeom::Point3D<double> &point) const;

  // Pathlength to helix
  //  virtual double getPathLengthTo(const Helix &helix);

  // Pathlength to line
  virtual double getPathLengthTo(const Line & line);

  // DZero to line
  virtual double getDzeroTo(const Line &line);

  // Intersection with a plane
  virtual Trajectory::Location * newIntersectionWith(const HepGeom::Plane3D<double> &plane) const;

  private:

      double _sinTheta;
      double _cosTheta;
      double _z0;
      double _d0;
      double _phi0;
      double _sinPhi0;
      double _cosPhi0;

};



} // end namespace mu2e

#endif /* CDFTrajectory_Line_hh */


