#ifndef CELLGEOMETRYHANDLE_V3_HH
#define CELLGEOMETRYHANDLE_V3_HH

#include "ITrackerGeom/inc/CellGeometryHandle.hh"
#include "ITrackerGeom/inc/Cell.hh"
#include "ITrackerGeom/inc/ITracker.hh"

#include "CLHEP/Geometry/Transform3D.h"
#include "CLHEP/Geometry/Point3D.h"

namespace mu2e {

class CellGeometryHandle_v3 : public CellGeometryHandle{

        friend class ITrackerMaker;

protected:
        CellGeometryHandle_v3(ITracker *itr=0x0);

public:

        ~CellGeometryHandle_v3();

        virtual void  SelectCell(int SupLayer, int CelLayer, int Cell);
        virtual void  SelectWireDet(unsigned long det);
        virtual unsigned long computeDet(int SupLayer, int CelLayer, int Cell);

protected:
    const ITracker *_itr;
    int _nLayer;

};

}

#endif /* CELLGEOMETRYHANDLE_V3_HH */
