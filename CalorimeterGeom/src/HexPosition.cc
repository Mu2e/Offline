//
// $Id: HexPosition.cc,v 1.2 2013/03/05 20:33:25 aluca Exp $
// $Author: aluca $
// $Date: 2013/03/05 20:33:25 $
//
// Hold information about position of a hexagonal cell
//
// Original author B Echenard - P. Ongmongkolkul
//


// C++ includes
#include <iostream>
#include <cmath>
#include "CLHEP/Vector/TwoVector.h"

// Mu2e includes
#include "CalorimeterGeom/inc/HexPosition.hh"




namespace mu2e {


     HexPosition::HexPosition():
     _l(0),_k(0),_x(0),_y(0),_dmin(1000000),_dmax(0),_scale(1){}


     HexPosition::HexPosition(int l, int k, double scale): 
     _l(l),_k(k),_x(0),_y(0),_dmin(1000000),_dmax(0),_scale(scale)
     {      
	_x = (_l+_k)*std::sqrt(3)/2.0;
	_y = (_l-_k)/2.0;
	calculateDistancesToOrigin();        
     }


     void HexPosition::calculateDistancesToOrigin(void)
     {    

	 std::pair<double,double> p1(_x-0.288675135 ,_y-0.5);
	 std::pair<double,double> p2(_x+0.288675135 ,_y-0.5);
	 std::pair<double,double> p3(_x+0.577350269 ,_y);
	 std::pair<double,double> p4(_x+0.288675135 ,_y+0.5);
	 std::pair<double,double> p5(_x-0.288675135 ,_y+0.5);
	 std::pair<double,double> p6(_x-0.577350269 ,_y);

	 double dn[6]={0};    
	 dn[0] = calcDistToLine(p1,p2);
	 dn[1] = calcDistToLine(p2,p3);
	 dn[2] = calcDistToLine(p3,p4);
	 dn[3] = calcDistToLine(p4,p5);
	 dn[4] = calcDistToLine(p5,p6);
	 dn[5] = calcDistToLine(p6,p1);

	 double dx[6]={0};
	 dx[0]  = calcNorm(_x+0.577350269 ,_y);
	 dx[1]  = calcNorm(_x-0.577350269 ,_y);    
	 dx[2]  = calcNorm(_x+0.288675135 ,_y+0.5);
	 dx[3]  = calcNorm(_x+0.288675135 ,_y-0.5);
	 dx[4]  = calcNorm(_x-0.288675135 ,_y+0.5);
	 dx[5]  = calcNorm(_x-0.288675135 ,_y-0.5);

	 for (int i=0;i<6;++i) {
	    if (dn[i]<_dmin) _dmin=dn[i];
	    if (dx[i]>_dmax) _dmax=dx[i];
	 } 
     }


     double HexPosition::calcDistToLine(std::pair<double,double>& a, std::pair<double,double>& b) const
     {
	 double t = -3.0*( (b.first-a.first)*a.first + (b.second-a.second)*a.second);
	 if (t<0) return calcNorm(a.first,a.second);
	 if (t>1) return calcNorm(b.first,b.second);
	 return calcNorm( a.first + t*(b.first-a.first), a.second + t*(b.second-a.second));
     }


     int HexPosition::ringNo(void) const
     {
	if (_l*_k>0) return std::abs(_l+_k); 
	return std::max(std::abs(_l),std::abs(_k));
     } 

     double HexPosition::calcNorm(double x, double y) const
     {
	return std::sqrt(x*x+y*y);
     } 

     double HexPosition::distanceTo(double x, double y) const
     {
	return _scale*calcNorm(x-_x,y-_y);
     } 

}
