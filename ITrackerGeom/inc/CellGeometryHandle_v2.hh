#ifndef CELLGEOMETRYHANDLE_V2_HH
#define CELLGEOMETRYHANDLE_V2_HH

#include "ITrackerGeom/inc/CellGeometryHandle.hh"
#include "ITrackerGeom/inc/ITracker.hh"

namespace mu2e {

class CellGeometryHandle_v2 : public CellGeometryHandle{

	friend class ITrackerMaker;

protected:
	CellGeometryHandle_v2(ITracker *itr=0x0);

public:

	~CellGeometryHandle_v2();

	virtual void  SelectCell(int SupLayer, int CelLayer, int Cell);

protected:
    const ITracker *_itr;

};

}

#endif /* CELLGEOMETRYHANDLE_V2_HH */
