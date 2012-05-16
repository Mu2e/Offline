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

    virtual void  SelectCell(int SupLayer, int CelLayer, int Cell, bool isUpstrm=false);
    virtual void  SelectCellDet(unsigned long det);
    virtual void  SelectCell(int absRadID, int Cell, bool isUpstrm=false);
    virtual unsigned long computeDet(int SupLayer, int CelLayer, int Cell, bool isUpstrm=false);
    virtual bool  canIntersectInZ(float &zCorss, float &distWires, unsigned long compDet) const;
    virtual bool  canIntersectInZ(float &zCorss, float &distWires, int compAbsRadID, int compICell, bool compIsUpstrm=false) const;
    virtual bool  canIntersectInZ(float &zCorss, float &distWires, int compSupLayer, int compCelLayer, int compCell, bool compIsUpstrm=false) const;

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
