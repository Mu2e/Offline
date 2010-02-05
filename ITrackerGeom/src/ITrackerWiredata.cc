#include "ITrackerGeom/inc/ITrackerWiredata.hh"

ITrackerWiredata::ITrackerWiredata(){
	PosMatrix = new TObjArray(0);
	NcelLayer = 0;
	epsilon = 0x0;
	alfa = 0x0;
	radius_z0 = 0x0;
}

ITrackerWiredata::~ITrackerWiredata(){
	PosMatrix->Delete();

	NcelLayer = 0;
	if (epsilon) delete [] epsilon;
	if (alfa) delete [] alfa;
	if (radius_z0) delete [] radius_z0;
}
