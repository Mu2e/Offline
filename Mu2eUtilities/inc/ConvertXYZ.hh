#ifndef _MU2E_UTILITIES_CONVERTXYZ_HH
#define _MU2E_UTILITIES_CONVERTXYZ_HH

#include "DataProducts/inc/XYZVec.hh"
#include "BTrk/BbrGeom/HepPoint.h"

using namespace mu2e;

HepPoint ConvertToHepPoint(XYZVec vec1){
	HepPoint point(vec1.x(),vec1.y(),vec1.z());
	return point;

}
Hep3Vector ConvertToHep3Vector(XYZVec vec1){
	Hep3Vector vec(vec1.x(),vec1.y(),vec1.z());
	return vec;

}
XYZVec ConvertToXYZ(Hep3Vector vec1){
	XYZVec XYZ(vec1.x(), vec1.y(), vec1.z());
	return XYZ;

}

XYZVec ConvertToXYZ(HepPoint vec1){
	XYZVec XYZ(vec1.x(), vec1.y(), vec1.z());
	return XYZ;

}

#endif
