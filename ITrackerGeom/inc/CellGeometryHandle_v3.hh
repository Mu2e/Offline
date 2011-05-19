#ifndef ITrackerGeom_CellGeometryHandle_v3_hh
#define ITrackerGeom_CellGeometryHandle_v3_hh

#include "ITrackerGeom/inc/CellGeometryHandle.hh"
#include "ITrackerGeom/inc/ITracker.hh"

namespace mu2e {

class CellGeometryHandle_v3 : public CellGeometryHandle{

    friend class ITrackerMaker;

protected:
    CellGeometryHandle_v3(ITracker *itr=0x0);

public:

    ~CellGeometryHandle_v3();

    virtual void  SelectCell(int SupLayer, int CelLayer, int Cell);
    virtual void  SelectCellDet(unsigned long det);
    virtual unsigned long computeDet(int SupLayer, int CelLayer, int Cell);

protected:
    const ITracker *_itr;
    int _nLayer;

private:
    // no copying:
    CellGeometryHandle_v3( CellGeometryHandle_v3 const & );
    void operator = ( CellGeometryHandle_v3 const & );

};

}

#endif /* ITrackerGeom_CellGeometryHandle_v3_hh */
