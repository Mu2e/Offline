// interface to manage the geometries of the ITracker cells
//
// $Id: CellGeometryHandle_v3_DBL.hh,v 1.3 2012/09/25 10:08:30 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/09/25 10:08:30 $
//
// Original author G. Tassielli
//

#ifndef ITrackerGeom_CellGeometryHandle_v3_DBL_hh
#define ITrackerGeom_CellGeometryHandle_v3_DBL_hh

#include "ITrackerGeom/inc/CellGeometryHandle_v3.hh"
#include "ITrackerGeom/inc/ITracker.hh"

namespace mu2e {

class CellGeometryHandle_v3_DBL : public CellGeometryHandle_v3{

    friend class ITrackerMaker;

protected:
    CellGeometryHandle_v3_DBL(ITracker *itr=0x0);

public:

    ~CellGeometryHandle_v3_DBL();

    virtual void  SelectCell(int SupLayer, int CelLayer, int Cell, bool isUpstrm=false);
    virtual void  SelectCellDet(unsigned long det);
    virtual void  SelectCell(int absRadID, int Cell, bool isUpstrm=false);
    virtual unsigned long computeDet(int SupLayer, int CelLayer, int Cell, bool isUpstrm=false);
    virtual int computeAbsRadID(int SupLayer, int CelLayer, bool isUpstrm=false);
    virtual const CLHEP::Hep3Vector& GetCellCenter() const;
    virtual const CLHEP::Hep3Vector& GetCellDirection() const;
    virtual double GetCellHalfLength() const;
    virtual void  SelectComp_Cell(int compSupLayer, int compCelLayer, int compICell, bool compIsUpStream=false);
    virtual void  SelectComp_CellDet(unsigned long compDet);
    virtual void  SelectComp_Cell(int compIAbsRadID, int compICell, bool compIsUpStream=false);
    virtual bool  canIntersectInZ(float &zCorss, float &distWires) const;

private:
    // no copying:
    CellGeometryHandle_v3_DBL( CellGeometryHandle_v3_DBL const & );
    void operator = ( CellGeometryHandle_v3_DBL const & );
    void  AdjustMatrix();
    int _nSLayer;
    double _DnStrmDeadWireLngt;
    double _UpStrmDeadWireLngt;
    double _tmpCellHalfLength;
    CLHEP::Hep3Vector _tmpDirection;
    CLHEP::Hep3Vector _tmpMidPoint;
    float wCntPos[3];

    inline void  Comp_AdjustMatrix();
    void  Comp_WirePosAtLength(float length, float *pos);
    bool _Comp_isUpStream;
    double _Comp_DnStrmDeadWireLngt;
    double _Comp_UpStrmDeadWireLngt;
    double _Comp_tmpCellHalfLength;
    CLHEP::Hep3Vector _Comp_tmpDirection;
    CLHEP::Hep3Vector _Comp_tmpMidPoint;
    float Comp_wCntPos[3];

};

}

#endif /* ITrackerGeom_CellGeometryHandle_v3_DBL_hh */
