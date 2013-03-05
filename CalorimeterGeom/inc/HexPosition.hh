#ifndef CalorimeterGeom_HexPosition_hh
#define CalorimeterGeom_HexPosition_hh
//
// $Id: HexPosition.hh,v 1.2 2013/03/05 20:33:25 aluca Exp $
// $Author: aluca $
// $Date: 2013/03/05 20:33:25 $
//
// Hold information about position of a hexagonal cell
//
// Original author B Echenard - P. Ongmongkolkul
//

#include <vector>
#include "CLHEP/Vector/TwoVector.h"


namespace mu2e {


     class HexPosition {


	 public:

	   HexPosition();
	   HexPosition(int l, int k, double scale=1.0);
 
	   int l(void) const                        {return _l;}
	   int k(void) const                        {return _k;}
	   int ringNo(void) const;

	   double distMin(void) const               {return _scale*_dmin;}
	   double distMax(void) const               {return _scale*_dmax;}
	   CLHEP::Hep2Vector XY(void) const            {return CLHEP::Hep2Vector(_scale*_x,_scale*_y);}
	   double distanceTo(double x, double y) const;


	 private:

	   int _l;
	   int _k;

	   double _x;
	   double _y;
	   double _dmin;
	   double _dmax;
	   double _scale;

	   void   calculateDistancesToOrigin(void);
	   double calcDistToLine(std::pair<double,double>& a, std::pair<double,double>& b) const;
	   double calcNorm(double x,double y) const;


     };


} //namespace mu2e


#endif /* CalorimeterGeom_HexPosition_hh */
