#include "CDFTrajectory/inc/Line.hh"

#include <math.h>

namespace mu2e {
double Line::getD0() const {
  return _d0;
}

double Line::getCotTheta() const {
  return _sinTheta==0 ? 0:_cosTheta/_sinTheta;
}

double Line::getPhi0() const {
  return _phi0;
}

double Line::getZ0() const {
  return _z0;
}

Line::Line()
  :_z0(0.0),
   _d0(0.0),
   _phi0(0.0)
{
  _sinTheta=1.0;
  _cosTheta=0.0;
  _sinPhi0=0.0;
  _cosPhi0=1.0;
}

Line::Line(const Line &right)
  :_sinTheta(right._sinTheta),
   _cosTheta(right._cosTheta),
   _z0(right._z0),
   _d0(right._d0),
   _phi0(right._phi0),
   _sinPhi0(right._sinPhi0),
   _cosPhi0(right._cosPhi0)
{
}

Line::Line(const HepGeom::Point3D<double> &point, const HepGeom::Vector3D<double> &direction_tmp, bool checkNormalisation) {

  // By default, renormalise the direction vector. Skip if checkNormalisation is set to false
  // for faster operation :
  HepGeom::Vector3D<double> direction(direction_tmp);
  if (checkNormalisation)
    {
      direction.setMag(1.0);
    }

  // This is the pathological case, parallel to z-axis.
  if (direction.perp2()<NUM_PREC)
    {
      _sinTheta=0.0;
      _cosTheta = direction.z()>0 ? 1.0:-1.0;
      _d0=point.perp();
      if (_d0 == 0) {
        _sinPhi0 = 0;
        _cosPhi0 = 0;
        _phi0    = 0;
      }
      else
        {
          if (fabs(_d0)>NUM_PREC)
            {
              _sinPhi0= -point.x()/_d0;
              _cosPhi0=  point.y()/_d0;
            }
          else
            {
              double perp = direction.perp();
              _sinPhi0 = direction.y()/perp;
              _cosPhi0 = direction.x()/perp;
            }
          _phi0=atan2(_sinPhi0,_cosPhi0);
        }
      _z0=0;
    }
  else if (fabs(direction.z())<NUM_PREC) // this is another one, perpendicular
    {
      _sinTheta = 1.0;
      _cosTheta = 0.0;
      _d0 = (direction.cross(point)).z();
      _phi0             = direction.phi();
      _sinPhi0          = sin(_phi0);
      _cosPhi0          = cos(_phi0);
      _z0 = point.z();
    }
  // This is the normal case
  else
    {
      CLHEP::Hep3Vector dPerp=direction;
      dPerp.setZ(0);
      dPerp.setMag(1.0);

      double dirz = direction.z();
      double theta = dirz>0 ? atan(direction.perp()/dirz) : M_PI + atan(direction.perp()/dirz);
      _sinTheta=sin(theta);
      _cosTheta=cos(theta);
      _d0               = (dPerp.cross(point)).z();
      _phi0             = direction.phi();
      _sinPhi0          = sin(_phi0);
      _cosPhi0          = cos(_phi0);

      double s                      = dPerp.dot(point);
      double sprime                 = s/_sinTheta;
      _z0                            = (point - sprime*direction).z();
#define  THESE_ARE_BOTH_BOGUS
#ifndef  THESE_ARE_BOTH_BOGUS
#ifdef ORIGINAL
      double x=-_d0*_sinPhi0;
      //    double y= _d0*_cosPhi0;
      double delta_x = x-point.x();
      //    double delta_y = y-point.y();
      double alpha_x = direction.x();
      //    double alpha_y = direction.y();
      double alpha_z = direction.z();
      _z0 = point.z() + (delta_x/alpha_x)*alpha_z;
#else
      double x=-_d0*_sinPhi0;
      double y= _d0*_cosPhi0;
      double delta_x = x-point.x();
      double delta_y = y-point.y();
      double trans = sqrt(delta_x*delta_x + delta_y*delta_y);
      _z0 = point.z() + trans*_cosTheta/_sinTheta;

#endif
#endif
    }


}


Line::~Line()
{
}


const Line & Line::operator=(const Line &right)
{
  if (this!=&right) {
    _sinTheta=right._sinTheta;
    _cosTheta=right._cosTheta;
    _z0=right._z0;
    _d0=right._d0;
    _phi0=right._phi0;
    _sinPhi0=right._sinPhi0;
    _cosPhi0=right._cosPhi0;
  }
  return *this;
}


HepGeom::Point3D<double> Line::getPosition(double s) const
{
  return HepGeom::Point3D<double>(-_d0*_sinPhi0+s*_cosPhi0*_sinTheta,_d0*_cosPhi0+s*_sinPhi0*_sinTheta,_z0+s*_cosTheta);
}

HepGeom::Vector3D<double> Line::getDirection(double ) const
{
  return HepGeom::Vector3D<double>(_cosPhi0*_sinTheta,_sinPhi0*_sinTheta,_cosTheta);
}

void Line::getLocation(Trajectory::Location & loc, double s) const {
  double cP0sT = _cosPhi0*_sinTheta, sP0sT = _sinPhi0*_sinTheta;
  loc.setLocation(s,
                  -_d0*_sinPhi0+s*cP0sT,
                  _d0*_cosPhi0+s*sP0sT,
                  _z0+s*_cosTheta,
                  cP0sT,sP0sT,_cosTheta);
}

double Line::getPathLengthAtRhoEquals(double rho) const {
  if ( _sinTheta != 0.0 ) {
    if ( rho > fabs( _d0 ) ) {
      return sqrt( rho*rho - _d0*_d0 ) / _sinTheta ;
    }
    else {
      return double(0);
    }
  }
  else {
    return double(0);
  }
}

Trajectory::Location * Line::newIntersectionWith(const HepGeom::Plane3D<double> &plane) const {
  if (plane.normal().dot(getDirection(0))!=0.0){
    double s=-(plane.distance(getPosition(0)))/(plane.normal().dot(getDirection(0)));
    Location * tmp = new Location();
    if (s>0)  {
      getLocation(*tmp,s);
      return tmp;
    }
    else {
      delete tmp;
    }
  }
  return nullptr;
}

HepGeom::Vector3D<double> Line::getSecondDerivative(double ) const{
  return HepGeom::Vector3D<double>(0,0,0);
}

//double Line::getDzeroTo(const HepGeom::Point3D<double> &point) const {
double Line::getPathLengthTo(const HepGeom::Point3D<double> &point) const {

  return (point - getPosition(0)).dot(getDirection(0));

}

//double Line::getPathLengthTo(const HepGeom::Point3D<double> &point) const {
double Line::getDzeroTo(const HepGeom::Point3D<double> &point) const {

  //  double opposite = getDzeroTo(point);
  double adjacent = getPathLengthTo(point);
  HepGeom::Point3D<double> hypvec = point-getPosition(0);
  double hypo = hypvec.mag();
  //  return sqrt(hypo*hypo-opposite*opposite);
  // Think this quantity is signed right. See Trajectory::getDzeroTo(const Trajectory &traj)
  //  return sqrt(hypo*hypo-adjacent*adjacent)*(getDirection(0).cross(hypvec)>0.)?-1.0:1.0;
  return sqrt(hypo*hypo-adjacent*adjacent);//*(getDirection(0).cross(hypvec)>0.)?-1.0:1.0;
}

// double Line::getPathLengthTo(const Helix &helix) {
//   // For speed, do a specific function for line/helix, and worry about the
//   // rest later.
//   // 24th Nov 1997. This is a *really* simple algorithm that assumes:
//   // a) the poca is along the positive helix-s direction;
//   // b) the poca is the first local minimum along the helix;
//   // Both these *should* be true for a sensible helix and a silicon strip
//   // line (r-phi or r-z).
//   // If Helix::getArcLengthAtRhoEquals were implemented, we'd have a far closer seed.

//   int ntries = 0;
//   const int MAXtries = 200;

//   // starting positions
//   HepGeom::Point3D<double> pointCurr;
//   HepGeom::Point3D<double> pointB = helix.getPosition(0);

//   double delta = 0.0;
//   double old_delta = 0.0;
//   HepGeom::Point3D<double> old_pointB; // may not need this

//   do {

//     old_delta=delta;

//     pointCurr=getPosition(getPathLengthTo(pointB));
//     old_pointB=pointB;
//     pointB=helix.getPosition(helix.getPathLengthTo(pointCurr));

//     delta = (pointB - pointCurr).mag();

//     if (ntries && (delta > old_delta)) {
//       // divide the distance between pointB and
//       // old_pointB into bits and work backwards
//       double start = getPathLengthTo(pointB);
//       double end = getPathLengthTo(old_pointB);
//       double step = EFINE*((end>start)?1.0:-1.0);
//       old_delta=delta;
//       for ( double posn = start+step; (end-posn)>step ; posn+=step ) {
//         delta = (pointB - pointCurr).mag();
//         if (delta < EEFINE || delta > old_delta) break;
//         old_delta=delta;
//       }

//       pointB = old_pointB;
//     }

//   } while (ntries++ < MAXtries && delta >= EEFINE);

//   if (ntries == MAXtries) {
//     std::cout << "Problem with Line.cc. MAXtries exceeded. Check or mail greenc@fnal.gov" << std::endl;
//   }

//   return getPathLengthTo(pointCurr);

// }

double Line::getPathLengthTo(const Line &line) {

  HepGeom::Point3D<double>   x1 = getPosition(0);
  HepGeom::Point3D<double>   x2 = line.getPosition(0);
  HepGeom::Point3D<double>  d1 = getDirection(0);
  HepGeom::Point3D<double>  d2 = line.getDirection(0);
  HepGeom::Point3D<double>   deltaX=x2-x1;
  double       d1Dotd2=d1.dot(d2);
  return (d1Dotd2*deltaX.dot(d2)-deltaX.dot(d1))/(d1Dotd2*d1Dotd2-1);

}


double Line::getDzeroTo(const Line &line) {

  // *think* this is signed right. See Trajectory::getDzeroTo(const Trajectory &traj)

  return getDirection(0).cross(line.getDirection(0)).unit().dot(line.getPosition(0)-getPosition(0));

}
}// end namespace mu2e
