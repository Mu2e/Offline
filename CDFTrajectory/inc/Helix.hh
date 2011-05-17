#ifndef CDFTrajectory_Helix_hh
#define CDFTrajectory_Helix_hh

#include "CDFTrajectory/inc/Trajectory.hh"
#include "CDFTrajectory/inc/Angle.hh"
//class HepVector;
//class HepPoint3D;
namespace mu2e {
class Helix : public Trajectory  {

public:

  //  Constructor
  inline Helix();

  // Construct from particle momentum, position, field, and charge.
  Helix(const HepGeom::Vector3D<double> &,
	const HepGeom::Point3D<double> &,
	double q, 
	double Field);

  // Copy Constructor
  inline Helix(const Helix &right);

  // Construct on five helix parameters
  inline Helix( 
        double cotTheta, 
        double curvature, 
        double z0, 
        double d0, 
        Angle  phi0
        );

  // Destructor
  virtual ~Helix();

  // Set Cot Theta
  inline void setCotTheta(double cotTheta);

  // Set Curvature
  inline void setCurvature(double curvature);

  // Set Z0 Parameter
  inline void setZ0(double z0);

  // Set D0 Parameter
  inline void setD0(double d0);

  // Set Phi0 Parameter, will be ranged from 0-2Pi
  inline void setPhi0(Angle phi0);

  // Assignment operator
  inline const Helix & operator=(const Helix &right);

  // Relational operators.
  bool operator==(const Helix & right) const;
  bool operator!=(const Helix & right) const;

  // Get Position as a function of (three-dimensional) path length
  virtual HepGeom::Point3D<double> getPosition(double s = 0.0) const;

  // Get Direction as a function of (three-dimensional) path length
  virtual HepGeom::Vector3D<double> getDirection(double s = 0.0) const;

  // Get the second derivative of the helix vs (three-dimensional) path length
  virtual HepGeom::Vector3D<double> getSecondDerivative(double s = 0.0) const;

  // Get both position and direction at once.
  virtual void getLocation(Trajectory::Location & loc, double s = 0.0) const;

  // Get pathlength at fixed rho=sqrt(x^2 + y^2)
  virtual double getPathLengthAtRhoEquals(double rho) const;
  
  //////////////////////////////////////////////////////////////////////////////////
  // KCDF: analytical computation of helix/plane intersection.
  //
  // What we really compute is the intersection of a line and 
  // a circle (projected helix) in the x-y plane.
  //
  //     >>>>>>>>>>  W A R N I N G    W A R N I N G    W A R N I N G  <<<<<<<<<<
  //     >                                                                     <
  //     > We assume the plane to be parallel or perpendicular                 <
  //     > to the z axis (i.e. the B-field),                                   <
  //     > since otherwise there is no analytical solution. (One would end up  <
  //     > with an equation of type cos(alpha) == alpha.)                      <
  //     > Although we know this assumption doesn´t hold exactly, we think it  <
  //     > is a reasonable first approximation.                                <
  //     > In cases when our assumption has to be dropped, one can use the     <
  //     > intersection point computed here as a *good* starting point of a    <
  //     > numerical method, or one uses another analytical method to compute  <
  //     > the intersection of the tangent line at the point and the plane.    <
  //     > We plan to use one of these approaches in the near future, but      <
  //     > this is NOT YET IMPLEMENTED!                                        <
  //     > For the time being, we invoke the old numerical                     < 
  //     > Trajectory::newIntersectionWith in such circumstances.              <
  //     >                                                                     <
  //     >>>>>>>>>>  W A R N I N G    W A R N I N G    W A R N I N G  <<<<<<<<<<
  //
  // Kurt Rinnert,  08/31/1998
  //////////////////////////////////////////////////////////////////////////////////
  Location* newIntersectionWith(const HepGeom::Plane3D<double> & plane) const;

  // Get certain parameters as a function of two-dimensional R.
  Angle  getPhiAtR(double r) const;
  double getZAtR(double r) const;
  double getL2DAtR(double r) const;
  double getCosAlphaAtR(double r) const;

  // Get signed 1/radius
  double getInverseRadius() const;

  // Get unsigned radius
  double getRadius() const;

  // Get the turning angle as a function of path length
  SignedAngle getTurningAngle(double s) const;

  // Get the Curvature
  double getCurvature() const;

  // Get helicity, positive for a counterclockwise helix
  double getHelicity() const;

  // Get cotangent of theta
  double getCotTheta() const;

  // Get phi0
  Angle getPhi0() const;

  // Get d0
  double getD0() const;

  // Get the Z0 parameter
  double getZ0() const;

  // Get sign of the z component of angular momentum about origin
  double getSignLz() const;

  // Get sines and cosines of Phi0 and Theta
  double getSinPhi0() const;
  double getCosPhi0() const;
  double getSinTheta() const;
  double getCosTheta() const;

  // Set Parameters.  
  inline void setParameters(const CLHEP::HepVector &p);
  /*
  // Get the parameters as a vector.
  inline const CLHEP::HepVector & getParameters() const;
  
  // Create a helix from a vector
  static Helix create(const CLHEP::HepVector & v);
  */
  // Return size of the parameter space (=5)
  static unsigned int getParameterSpaceSize();

private:


  // This is the Helix:
  double _cotTheta;
  double _curvature;
  double _z0;
  double _d0;
  Angle  _phi0;

  // This is the cache
  mutable bool   _isStale;
  mutable double _sinPhi0;
  mutable double _cosPhi0;
  mutable double _sinTheta;
  mutable double _cosTheta;
  mutable double _s;
  mutable double _aa;
  mutable double _ss;
  mutable double _cc;
  mutable CLHEP::HepVector  * _vParameters;
  mutable bool _centerIsValid; //needed by newIntersectionWith KR
  mutable double _m_x;
  mutable double _m_y;

  //needed whenever _sinPhi0, _cosPh0, _sinTheta, or _cosTheta is used.
  inline void _refreshCache() const;

  // neede whenever _ss or _cc are used.
  inline void _cacheSinesAndCosines(double s) const;
  
};

#include "CDFTrajectory/inc/Helix.icc"
} // end namespace mu2e
#endif /* CDFTrajectory_Helix_hh */


