#include "CDFTrajectory/inc/Helix.hh"
#include "CDFTrajectory/inc/Line.hh"
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Vector3D.h"
// #include "ErrorLogger_i/gERRLOG.hh"

#ifdef DEFECT_OLD_STDC_HEADERS
extern "C" {
#include <math.h>
}
#else
#include <cmath>
#endif
namespace mu2e {
bool Helix::operator==( const Helix & right) const {
  return _cotTheta==right._cotTheta &&
         _curvature==right._curvature &&
         _z0==right._z0 &&
         _d0==right._d0 &&
         _phi0==right._phi0;
}

bool Helix::operator!=( const Helix & right) const {
  return !((*this)==right);
}



Helix::~Helix()
{
  delete _vParameters;
}



HepGeom::Point3D<double> Helix::getPosition(double s) const
{
  _cacheSinesAndCosines(s);
  if (s==0.0 || _curvature==0.0) {
      return HepGeom::Point3D<double>(-_d0*_sinPhi0+s*_cosPhi0*_sinTheta,
                        _d0*_cosPhi0+s*_sinPhi0*_sinTheta,
                        _z0+s*_cosTheta);
  } else {
      return HepGeom::Point3D<double>((_cosPhi0*_ss-_sinPhi0*(2.0*_curvature*_d0+1.0-_cc))
                        /(2.0*_curvature),
                        (_sinPhi0*_ss+_cosPhi0*(2.0*_curvature*_d0+1.0-_cc))
                        /(2.0*_curvature),
                        _s*_cosTheta + _z0);
  }
}

HepGeom::Vector3D<double> Helix::getDirection(double s) const
{
  _cacheSinesAndCosines(s);
  if (s==0.0) {
      return HepGeom::Vector3D<double>(_cosPhi0*_sinTheta,_sinPhi0*_sinTheta,_cosTheta);
  }
  double   xtan     = _sinTheta*(_cosPhi0*_cc -_sinPhi0*_ss);
  double   ytan     = _sinTheta*(_cosPhi0*_ss +_sinPhi0*_cc);
  double ztan       = _cosTheta;
    return HepGeom::Vector3D<double>(xtan,ytan,ztan);
}

void Helix::getLocation(Trajectory::Location & loc, double s) const {

  _cacheSinesAndCosines(s);
  double cP0sT = _cosPhi0*_sinTheta, sP0sT = _sinPhi0*_sinTheta;
  if ( s && _curvature) {
    loc.setLocation(s,
                    (_cosPhi0*_ss-
                     _sinPhi0*(2.0*_curvature*_d0+1.0-_cc))
                    /(2.0*_curvature),
                    (_sinPhi0*_ss+
                     _cosPhi0*(2.0*_curvature*_d0+1.0-_cc))
                    /(2.0*_curvature),
                    s*_cosTheta + _z0,
                    cP0sT*_cc-sP0sT*_ss,
                    cP0sT*_ss+sP0sT*_cc,
                    _cosTheta
                    );
  } else {
    loc.setLocation(s,
                    -_d0*_sinPhi0+s*cP0sT,
                    _d0*_cosPhi0+s*sP0sT,
                    _z0+s*_cosTheta,
                    cP0sT,sP0sT,_cosTheta
                    );
  }
}

double Helix::getPathLengthAtRhoEquals(double rho) const
{
  return (getSinTheta()?(getL2DAtR(rho)/getSinTheta()):0.0);
}

double Helix::getInverseRadius() const
{
  return _curvature*2.0;
}

double Helix::getRadius() const
{
  return fabs(1.0/getInverseRadius());
}

SignedAngle Helix::getTurningAngle(double s) const
{
  return s/getRadius();
}

double Helix::getCurvature() const
{
  return _curvature;
}

double Helix::getHelicity() const
{
  return _curvature>0 ? 1.0 : -1.0 ;
}

double Helix::getCotTheta() const
{
  return _cotTheta;
}

Angle Helix::getPhi0() const
{
  return _phi0;
}

double Helix::getD0() const
{
  return _d0;
}

double Helix::getSignLz() const
{
  return _d0>0 ? -1.0 : 1.0 ;
}

double Helix::getZ0() const
{
  return _z0;
}


HepGeom::Vector3D<double> Helix::getSecondDerivative(double s) const{
  double phi1= _phi0+s*2.0*_curvature*_sinTheta;
  double sp1=sin(phi1);
  double   xsecond     = -_sinTheta*sp1*2.0*_curvature*_sinTheta;
  double   ysecond     = _sinTheta*sqrt(1.0-sp1*sp1)*2.0*_curvature*_sinTheta;
  return HepGeom::Vector3D<double>(xsecond,ysecond,0.0);
}

double Helix::getSinPhi0() const{
  _refreshCache();
  return _sinPhi0;
}
double Helix::getCosPhi0() const{
  _refreshCache();
  return _cosPhi0;
}
double Helix::getSinTheta() const{
  _refreshCache();
  return _sinTheta;
}
double Helix::getCosTheta() const{
  _refreshCache();
  return _cosTheta;
}

Angle Helix::getPhiAtR(double rho) const {
  double c = getCurvature();
  double d = getD0();
  double a = c/(1.0+2.0*c*d);
  double b = d - a*d*d;
  double arcsin = a*rho+b/rho;
  if (arcsin>1.0 || arcsin<-1.0) {
    //---- 3-Mar-2002 (RS):  low-level code should not call error logger
    // // Warning if significantly out of range.
    // double diff = abs( arcsin ) - 1.0;
    // if ( diff > 1.0E-8 )
    // {
    //   errlog.setSubroutine("Helix::getPhiAtR(rho)");
    //   errlog(ELwarning, "[rho out of bounds]")
    //         << " arcsin = " << arcsin << endmsg;
    // }
    //------------------------------------------------------------------
    return arcsin > 0. ? M_PI : -M_PI;
  }
  Angle phi = getPhi0() + asin(arcsin);
  return phi;
}

double Helix::getZAtR(double rho) const {
  return _z0 + getCotTheta()*getL2DAtR(rho);
}

double Helix::getL2DAtR(double rho) const {
  double L2D;

  double c = getCurvature();
  double d = getD0();
  // errlog.setSubroutine("Helix::getZAtR(rho)");

  if (c!=0.0) {
    double rad = (rho*rho-d*d)/(1.0+2.0*c*d);
    double rprime;
    if (rad<0.0) {
      //---- 3-Mar-2002 (RS):  low-level code should not call error logger
      // // Warning if significantly less than zero (=> testing at radius
      // // smaller than closest approach). Otherwise, assume it should be zero.
      // // Tolerance of 10^-8 => error of order 1 um at 300 GeV
      // if ( rad < -1.0E-8 ) {
      //   errlog(ELwarning, "[rho out of bounds]")
      //          << "@SUB=Helix::getL2DAtR(rho)"
      //          << " rho = " << rho
      //          << ", d0 = " << d << endmsg;
      // }
      //------------------------------------------------------------------
      rprime = 0.0;
    }
    else {
      rprime = sqrt( rad );
    }
    if (c*rprime > 1.0 || c*rprime < -1.0) {
      //---- 3-Mar-2002 (RS):  low-level code should not call error logger
      // // Warning if product is significantly out of range (=> testing at
      // // radius larger than apogee). Otherwise, assume it should be unity.
      // if ( abs( c*rprime ) - 1.0 > 1.0E-8 ) {
      //   errlog(ELwarning, "[workaround illegal asin argument]")
      //          << "@SUB=Helix::getL2DAtR(rho)"
      //          << " c = " << c
      //          << ", rprime = " << rprime
      //          << ", (c*rprime)-1 = " << (c * rprime ) - 1.0 << endmsg;
      // }
      //------------------------------------------------------------------
      L2D = c*rprime > 0. ? M_PI/c : -M_PI/c;
    }
    else
      L2D = asin(c*rprime)/c;
  } else {
    //     L2D = rho;
    double rad = rho*rho-d*d;
    double rprime;
    if (rad<0.0) rprime = 0.0;
    else rprime = sqrt( rad );

    L2D = rprime;
  }
  return L2D;
}

double Helix::getCosAlphaAtR(double rho) const {
  double c = getCurvature();
  double d = getD0();
  double a = c/(1.0+2.0*c*d);
  //AML: not necessary? not used: double b = d - a*d*d;

  double phi = getPhiAtR(rho);

  double x = rho*cos(phi);
  double y = rho*sin(phi);

  double dir_x0 = (1.0+2.0*c*d) * (1.0-2.0*a*y);
  double dir_y0 = (1.0+2.0*c*d) * 2.0*a*x;

  double phi0 = getPhi0();
  double dir_x = dir_x0*cos(phi0) - dir_y0*sin(phi0);
  double dir_y = dir_x0*sin(phi0) + dir_y0*cos(phi0);

  double cosalpha = (x*dir_x + y*dir_y)/rho;
  return cosalpha;
}

///////////////////////////////////////////////////////////////////////////////////////
// KCDF: analytical computation of helix/plane intersection.
//
// What we really compute is the intersection of a line and a circle (projected helix)
// in the x-y plane.
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
////////////////////////////////////////////////////////////////////////////////////////
Trajectory::Location* Helix::newIntersectionWith (const HepGeom::Plane3D<double> & plane) const
{
//-------------------------------------------------------------
// (Weiming 8-Aug-2001) Commented out analytical calculation
// See comments below.
//-------------------------------------------------------------
//
//   if (!getSinTheta()) return NULL; // fastest way out of a screwy situation.
//
//   double alpha0, alpha1, alpha2, deltaAlpha1, deltaAlpha2; //angles measured around center of circle
//   double PHI0; //global phi coordinate
//   double n_x, n_y; //defines plane orientation
//   double p_x, p_y; //point in plane
//   double r; //radius of circle
//   double s_x1, s_y1; //first intersection point
//   double s_x2, s_y2; //second intersection point
//   double A, B, C, B2_4AC, SQRTB2_4AC; //some quantities to simplify notation
//   double t1, t2; //line parameters
//
//   double sign = getHelicity();
//
//   HepGeom::Vector3D<double> normal = plane.normal();
//
//   n_x = normal.x();
//   n_y = normal.y();
//
//   if (fabs(normal.z()) > 1.0E-6) {
//
//     //see WARNING above...
//     if((fabs(n_x) > 1.0E-6) || (fabs(n_y) > 1.0E-6)){
//       return Trajectory::newIntersectionWith(plane);
//
//     //the perpendicular case
//     }else{
//       r = 0.5 / fabs(_curvature); // 1 / | 2 kappa |
//
//       PHI0 = _phi0 + sign * 0.5 * M_PI;
//       alpha0 = _phi0 - sign * M_PI * 0.5;
//
//       // compute center of circle (only if not already done)
//       if(!_centerIsValid){
// 	_m_x = (sign*_d0 + r) * cos(PHI0);
// 	_m_y = (sign*_d0 + r) * sin(PHI0);
// 	_centerIsValid = true;
//       }
//
//       double s_z = plane.point().z();
//       double s = (s_z - _z0)/(_cotTheta);
//       double alpha = sign*s/r + alpha0;
//
//       double s_x = _m_x + r*cos(alpha);
//       double s_y = _m_y + r*sin(alpha);
//
//       return new Trajectory::Location(s/getSinTheta(),
//                                       s_x, s_y, s_z);
//
//     }
//   }
//
//   r = 0.5 / fabs(_curvature); // 1 / | 2 kappa |
//
//   p_x = plane.point().x();
//   p_y = plane.point().y();
//
//   alpha0 = _phi0 - sign * M_PI * 0.5;
//   if(alpha0 < 0.0){
//     alpha0 += 2.0*M_PI;
//   }else if(alpha0 > 2.0*M_PI){
//     alpha0 -= 2.0*M_PI;
//   };
//
//   PHI0 = _phi0 + sign * 0.5 * M_PI;
//
//   // compute center of circle (only if not already done)
//   if(!_centerIsValid){
//     _m_x = (sign*_d0 + r) * cos(PHI0);
//     _m_y = (sign*_d0 + r) * sin(PHI0);
//     _centerIsValid = true;
//   }
//
//   // coefficients of quadratic equation
//   A = n_x*n_x + n_y*n_y;
//   B = 2.0*(n_y*p_x - n_y*_m_x - n_x*p_y + n_x*_m_y);
//   C = (p_x - _m_x)*(p_x - _m_x) + (p_y - _m_y)*(p_y - _m_y) - r*r;
//
//   B2_4AC = B*B - 4.0*A*C;
//
//   //ther may be no solution...
//   if (B2_4AC >= 0.0){
//     SQRTB2_4AC = sqrt(B2_4AC);
//
//     // the two solutions
//     t1 =  (-B + SQRTB2_4AC)/(2.0*A);
//     t2 =  (-B - SQRTB2_4AC)/(2.0*A);
//
//     // compute first intersection
//     s_x1 = p_x + t1*n_y;
//     s_y1 = p_y - t1*n_x;
//     double deltaX1 = s_x1 - _m_x;
//     double deltaY1 = s_y1 - _m_y;
//     alpha1 = ( deltaX1 ?
//                atan2(deltaY1,deltaX1):
//                ( deltaY1 > 0.0 ? M_PI*0.5 : -M_PI*0.5));
//     alpha1 = ( alpha1 < 0 ? alpha1 + 2*M_PI : alpha1 );
//     deltaAlpha1 = sign*(alpha1 - alpha0); //note the sign!
//
//     // compute second intersection
//     s_x2 = p_x + t2*n_y;
//     s_y2 = p_y - t2*n_x;
//     double deltaX2 = s_x2 - _m_x;
//     double deltaY2 = s_y2 - _m_y;
//     alpha2 = ( deltaX2 ?
//                atan2(deltaY2,deltaX2):
//                ( deltaY2 > 0.0 ? M_PI*0.5 : -M_PI*0.5));
//     alpha2 = ( alpha2 < 0 ? alpha2 + 2*M_PI : alpha2 );
//     deltaAlpha2 = sign*(alpha2 - alpha0); //note the sign!
//
//     // We choose the solution with the smallest arclength measured from the perihel in
//     // the correct direction. This boils down to the difference in alpha.
//
//     if((deltaAlpha2 - deltaAlpha1) > 0.0 ) {
//       if(deltaAlpha1 > 0.0) {
//         double s = r*deltaAlpha1;
//         return new Trajectory::Location(s/getSinTheta(), s_x1, s_y1, _z0 + _cotTheta*s);
//       }else if(deltaAlpha2 > 0.0) {
//         double s = r*deltaAlpha2;
//         return new Trajectory::Location(s/getSinTheta(), s_x2, s_y2, _z0 + _cotTheta*s);
//       }else{
//         double s = r*(deltaAlpha1 + 2.0*M_PI);
//         return new Trajectory::Location(s/getSinTheta(), s_x1, s_y1, _z0 + _cotTheta*s);
//       }
//     }else{
//       if(deltaAlpha2 > 0.0) {
//         double s = r*deltaAlpha2;
//         return new Trajectory::Location(s/getSinTheta(), s_x2, s_y2, _z0 + _cotTheta*s);
//       }else if(deltaAlpha1 > 0.0) {
//         double s = r*deltaAlpha1;
//         return new Trajectory::Location(s/getSinTheta(), s_x1, s_y1, _z0 + _cotTheta*s);
//       }else{
//         double s = r*(deltaAlpha2 + 2.0*M_PI);
//         return new Trajectory::Location(s/getSinTheta(), s_x2, s_y2, _z0 + _cotTheta*s);
//       }
//     }
//   }else{
//     return NULL;
//   }

//***************************************************************************
// Weiming (8-Aug-2001):
// copied from Trajectory.cc. It seems has much better quality of projections
//***************************************************************************

  if (!getSinTheta()) return NULL; // fastest way out of a screwy situation.

  double deltaS=1.0,s=0.0;
  s = fabs(plane.d())/getSinTheta();
  HepGeom::Vector3D<double> normal = plane.normal();
  Trajectory::Location *ploc = new Trajectory::Location ;
  for (int iteration=0;iteration<100;iteration++) {
    if (fabs(deltaS)>0.0001){
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


Helix::Helix(const HepGeom::Vector3D<double> & MomentumGev,
	     const HepGeom::Point3D<double>  & PositionCm,
	     double q,
	     double BFieldTesla) {


  double CotTheta,W,Z0,D0;
  Angle Phi0;
  if (BFieldTesla != 0.0 && q != 0.0) {
    double CurvatureConstant=0.0029979;
    double Helicity       = -1.0 * fabs(BFieldTesla) * fabs(q) / (BFieldTesla*q);
    double Radius         = fabs(MomentumGev.perp()/(CurvatureConstant*BFieldTesla*q));

    if(Radius==0.0) {
      W = HUGE_VAL;
      CotTheta = 0.0;
    } else {
      W = Helicity/Radius;
      CotTheta = MomentumGev.z()/MomentumGev.perp();
    }
    Angle phi1      = MomentumGev.phi();
    double x        = PositionCm.x(),        y        = PositionCm.y(),    z = PositionCm.z();
    double sinPhi1  = sin(phi1),             cosPhi1  = cos(phi1);
    double gamma    = atan((x*cosPhi1 + y*sinPhi1)/(x*sinPhi1-y*cosPhi1 -1/W));
    Phi0            = phi1+gamma;
    D0              = ((1/W + y*cosPhi1 - x*sinPhi1) /cos(gamma) - 1/W);
    Z0              = z + gamma*CotTheta/W;
  }
  else {
    Line  auxLine(PositionCm,MomentumGev.unit());
    Z0       = auxLine.getZ0();
    Phi0     = auxLine.getPhi0();
    CotTheta = auxLine.getCotTheta();
    D0       = auxLine.getD0();
    W        = 0.0;
    //    Hep3Vector direction          = MomentumGev.unit();
    //    Hep3Vector projectedDirection = Hep3Vector(direction.x(),direction.y(),0.0).unit();
    //    double s                      = projectedDirection.dot(PositionCm);
    //    double sprime                 = s/sin(direction.theta());
    //    Z0                            = (PositionCm - sprime*direction).z();
    //    Phi0                          = MomentumGev.phi();
    //    CotTheta                      = MomentumGev.z()/MomentumGev.perp();
    //    W                             = 0.0;
    //    D0                           = (PositionCm.y()*cos(Phi0) - PositionCm.x()*sin(Phi0));
  }

  // The line commented out below doesn't appear to work (Msm and Dsw Jul 2000)
  //  Helix(CotTheta,W/2,Z0,D0,Phi0);
  // So replace it with the contents of the constructor
   _cotTheta=CotTheta;
   _curvature=W/2;
   _z0=Z0;
   _d0=D0;
   _phi0=Phi0;
   _isStale=1;
   _s=-999.999;
   _aa =-999.999;
   _ss =-999.999;
   _cc =-999.999;
   _sinPhi0 = 1.0;
   _cosPhi0 = 1.0;
   _sinTheta = 1.0;
   _cosTheta = 1.0;
   _vParameters=NULL;
   _centerIsValid=false;
   _m_x=0.0;
   _m_y=0.0;
}

} //namespace mu2e
