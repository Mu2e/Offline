// interface to manage the geometries of the ITracker cells
//
// $Id: CellGeometryHandle_v3.hh,v 1.9 2012/09/25 10:06:54 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/09/25 10:06:54 $
//
// Original author G. Tassielli
//

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
    virtual void  SelectComp_Cell(int compSupLayer, int compCelLayer, int compICell, bool compIsUpStream=false);
    virtual void  SelectComp_CellDet(unsigned long compDet);
    virtual void  SelectComp_Cell(int compIAbsRadID, int compICell, bool compIsUpStream=false);
    virtual bool  canIntersectInZ(float &zCorss, float &distWires) const;
    virtual bool  canIntersectInZ(float &zCorss, float &distWires, unsigned long compDet) const;
    virtual bool  canIntersectInZ(float &zCorss, float &distWires, int compAbsRadID, int compICell, bool compIsUpstrm=false) const;
    virtual bool  canIntersectInZ(float &zCorss, float &distWires, int compSupLayer, int compCelLayer, int compCell, bool compIsUpstrm=false) const;

    virtual double CrossingPathOnFieldWires(CLHEP::Hep3Vector const &point, CLHEP::Hep3Vector const &dir,
                                            FWireSide &sideFlag, CLHEP::Hep3Vector &fwPca, double tolerance=0.002) const;
    virtual double CrossingPathOnSenseWires(CLHEP::Hep3Vector const &point, CLHEP::Hep3Vector const &dir,
                                            CLHEP::Hep3Vector &swPca, double tolerance=0.002) const;

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
