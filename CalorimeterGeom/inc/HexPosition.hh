#ifndef CalorimeterGeom_HexPosition_hh
#define CalorimeterGeom_HexPosition_hh

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
 
	   int getl(void) const                        {return _l;}
	   int getk(void) const                        {return _k;}
	   int getRingNo(void) const;                   

	   double getDistMin(void) const               {return _scale*_dmin;}
	   double getDistMax(void) const               {return _scale*_dmax;}   
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
