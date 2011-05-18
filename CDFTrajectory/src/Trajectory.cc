
#include "CDFTrajectory/inc/Trajectory.hh"
#include "CDFTrajectory/inc/Line.hh"

#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
using std::ios_base;
using std::streamsize;
using std::endl;
using std::cerr;

namespace mu2e {
const double Trajectory::ECOARSE = 10.0;
const double Trajectory::COARSE  =  1.0;
const double Trajectory::FINE    =  0.01;
const double Trajectory::EFINE   =  0.00001;
const double Trajectory::EEFINE  =  0.000001;
const double Trajectory::COMP_PREC  =  EEFINE/10.;

const double Trajectory::NUM_PREC = 1.0e-15; // empirical

// Class Trajectory

Trajectory::~Trajectory()
{
}

// Return a Location object containing position and location data. This
// should be overriden in subclasses with obtimized versions -- the whole
// point is to cut down on calculations common to getPosition and
// getDirection.
void Trajectory::getLocation(Trajectory::Location & loc, double s) const {
  // If you find you're calling this default (and *slow*) routine, then
  // implement an optimized version for the trajectory concerned. See
  // Helix::getLocation(...) for an example.
  HepGeom::Point3D<double> pos=getPosition(s);
  CLHEP::Hep3Vector dir=getDirection(s);
  loc.setLocation(s,pos,dir);
}

Line Trajectory::getTangent(double s) const {
  Trajectory::Location location;
  getLocation(location,s);
  return Line(location.position(),location.direction());
}

double Trajectory::getPathLengthAtRhoEquals(double rho) const
{
  return 0;
}

double Trajectory::getPathLengthTo(const HepGeom::Point3D<double> &point) const {
  double s = 0;
  const int MAXTRIES=200;
  Trajectory::Location location;
  //  for (int i=0;i<10;i++) {
  for (int ntries=0 ; ntries<MAXTRIES ; ntries++) {
    getLocation(location,s);
    HepGeom::Point3D<double> xp = location.position();
    double DeltaS = (point -xp)*location.direction();
    //    if ((float) fabs(s+DeltaS)!=(float) fabs(s)) {
    if (std::abs(s+DeltaS)!=std::abs(s)) {
      s += DeltaS;
    }
    else {
      return s;
    }
  }
  return 0;
}




double Trajectory::getPathLengthTo(const Trajectory &traj) const {

  //-------------------------------------------------------------------------//
  //                                                                         //
  //   A bit more complicated than the basic method, {see                    //
  //   Line::getPathLengthTo(const Line &) const} because the procedure      //
  //   becomes iterative on arbitrary trajectories whose tangent vector      //
  //   varies with path length.  But basically the same calculation within a //
  //   loop.                                                                 //
  //                                                                         //
  //-------------------------------------------------------------------------//

  double s0=0,s1=0,s0Prev=0,s1Prev=0,dXPrev=std::numeric_limits<double>::max(),f=1.0;
  Trajectory::Location x0Loc, x1Loc;
  for (int i=0;i<20;i++) {
    getLocation(x0Loc,s0);
    traj.getLocation(x1Loc,s1);
    HepGeom::Point3D<double>   deltaX =x1Loc.position()-x0Loc.position();
    double       d0Dotd1=x0Loc.direction().dot(x1Loc.direction());
    double       delta0 = (d0Dotd1*deltaX.dot(x1Loc.direction())-
                           deltaX.dot(x0Loc.direction()))/(d0Dotd1*d0Dotd1-1);
    double delta1 = (-d0Dotd1*deltaX.dot(x0Loc.direction())+
                     deltaX.dot(x1Loc.direction()))/(d0Dotd1*d0Dotd1-1);

    if (std::abs(delta0)<EFINE && std::abs(delta1)<EFINE)
      {
        return s0;
      }
    double dxMag = deltaX.mag();

    // We have in fact increased.  Do not move forward, move back, reduce
    // Reduce the factor and try again..
    //
    if (dxMag>dXPrev) {
      f /= 2.0;
      s0 = s0Prev;
      s1 = s1Prev;
      continue;
    }

    s0+=(f*delta0);
    s1+=(f*delta1);
    s0Prev=s0;
    s1Prev=s1;
    dXPrev=dxMag;
  }
  return s0;
}


double Trajectory::getDzeroTo(const HepGeom::Point3D<double> &point) const {
  return 0;
}


// This is a signed quantity.  The sign is the *opposite* of the
// angular momentum of this trajectory about the POCA of the
// second trajectory, in the direction of the second
// trajectory. It's the same convention as the CDF D0, if you
// consider the D0 to be the signed POCA w.r.t a trajectory
// moving towards the positive z axis along the beam.
//
// Oh, subclasses!  You are exhorted to follow this convention!
// Exhort! Exhort! Exhort! Exhort! Exhort! Exhort! Exhort! Exhort!
double Trajectory::getDzeroTo(const Trajectory &traj) const {

  double s0=0,s1=0;
  Trajectory::Location x0Loc,x1Loc;
  for (int i=0;i<100;i++) {
    getLocation(x0Loc,s0);
    traj.getLocation(x1Loc,s1);
    HepGeom::Point3D<double>   deltaX=x1Loc.position()-x0Loc.position();
    double       d0Dotd1=x0Loc.direction().dot(x1Loc.direction());
    double       delta0 = (d0Dotd1*deltaX.dot(x1Loc.direction())
                           -deltaX.dot(x0Loc.direction()))/(d0Dotd1*d0Dotd1-1);
    double delta1 = (-d0Dotd1*deltaX.dot(x0Loc.direction())
                     +deltaX.dot(x1Loc.direction()))/(d0Dotd1*d0Dotd1-1);

    if (std::abs(delta0)<EFINE && std::abs(delta1)<EFINE)
      {
        CLHEP::Hep3Vector R=x0Loc.position()-x1Loc.position();
        CLHEP::Hep3Vector likeZ=x1Loc.direction();
        CLHEP::Hep3Vector rotR=CLHEP::Hep3Vector(x0Loc.direction()).cross(R);
        double sgn = rotR.dot(likeZ) > 0 ? 1.0:-1.0;
        return sgn*rotR.mag();
        //        return s0;
      }
    s0+=delta0;
    s1+=delta1;
  }

  return 0.0;

}

Trajectory::Location * Trajectory::newIntersectionWith(const HepGeom::Plane3D<double> &plane) const
{
  double deltaS=1.0,s=0.0;
  HepGeom::Vector3D<double> normal = plane.normal();
  Trajectory::Location *ploc = new Trajectory::Location ;
  for (int iteration=0;iteration<100;iteration++) {
    //    if (((float)fabs(s+deltaS)!=(float)fabs(s))){
    if (std::abs(s+deltaS)!=std::abs(s)){
      getLocation(*ploc,s);
      deltaS=((-(plane.distance(ploc->position())))/(normal.dot(ploc->direction())));
      s+=deltaS;
    }
    else {
      if (s>0) return ploc;
    }
  }
  delete ploc;
  return NULL;
}

void Trajectory::Location::print(std::ostream & os) const {

#ifdef DEFECT_OLD_IOSTREAM_HEADERS
  const iosbasefmtflags oldopts = os.flags();
#else
  const ios_base::fmtflags oldopts = os.flags();
#endif
  os.setf(ios_base::fixed,ios_base::floatfield); // set fixed format
  streamsize old_prec = os.precision(4); // save old and set new

  os << __s ;
  os << ": ";
  os << _position;
  os << ", ";
  os << _direction;
  os << std::endl;

  os.precision(old_prec);
  os.flags(oldopts);
}

// Backed out by CG 1998/09/10
//
// Storage of a pointer which may or may not disappear at some point in the
// future is too dangerous for the advantage it gives.
//
// Trajectory::Location::direction(void) has been re-inlined.
////////////////////////////////////////////////////////////////
//KCDF:
//
//We modified this to avoid computation of directions
//where they are never needed. (See Trajectory.hh)
//e.g. in Helix::newIntersectionWith(plane)
//
//For this purpose a private member _dirIsComputed (bool) and
//a pointer to a trajectory is introduced.
//The direction is computed only the first time it is accessed.
//(If at all!)
//
//When the old constructors are used Trajectory::Location
//behaves like before, i.e. nothing is computed.
//
//Kurt Rinnert 09/09/1998
////////////////////////////////////////////////////////////////
// const HepGeom::Vector3D<double>& Trajectory::Location::direction(void) const {
//   if(!_dirIsComputed){
//     _direction = _traj->getDirection(__s);
//     _dirIsComputed = true;
//   }
//   return _direction;
// }

// #ifdef TRACKING_DEBUG_WAS_ON
// #undef NO_TRACKING_DEBUG
// #undef TRACKING_DEBUG_WAS_ON
// #endif

} //namespace mu2e
