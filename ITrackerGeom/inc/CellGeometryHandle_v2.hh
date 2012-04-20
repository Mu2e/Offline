#ifndef ITrackerGeom_CellGeometryHandle_v2_hh
#define ITrackerGeom_CellGeometryHandle_v2_hh

#include "ITrackerGeom/inc/CellGeometryHandle.hh"
#include "ITrackerGeom/inc/ITracker.hh"

namespace mu2e {

class CellGeometryHandle_v2 : public CellGeometryHandle{

        friend class ITrackerMaker;

protected:
        CellGeometryHandle_v2(ITracker *itr=0x0);

public:

        ~CellGeometryHandle_v2();

        virtual void  SelectCell(int SupLayer, int CelLayer, int Cell, bool isUpstrm=false);
        virtual void  SelectCell(int absRadID, int Cell, bool isUpstrm=false);
        virtual int computeAbsRadID(int SupLayer, int CelLayer, bool isUpstrm=false);

protected:
        const ITracker *_itr;
        void  SelectCell();
};

}

#endif /* ITrackerGeom_CellGeometryHandle_v2_hh */
