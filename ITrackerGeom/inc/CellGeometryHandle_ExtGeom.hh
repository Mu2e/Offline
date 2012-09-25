// interface to manage the geometries of the ITracker cells
//
// $Id: CellGeometryHandle_ExtGeom.hh,v 1.6 2012/09/25 10:08:30 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/09/25 10:08:30 $
//
// Original author G. Tassielli
//

#ifndef ITrackerGeom_CellGeometryHandle_ExtGeom_hh
#define ITrackerGeom_CellGeometryHandle_ExtGeom_hh

#include "ITrackerGeom/inc/CellGeometryHandle.hh"
#include "ITrackerGeom/inc/ITrackerWireposition.hh"

namespace mu2e {

class CellGeometryHandle_ExtGeom : public CellGeometryHandle {

        friend class ITrackerMaker;

protected:
        CellGeometryHandle_ExtGeom(const char *WireDataFile = "ITrackerWireData.root");

public:

        ~CellGeometryHandle_ExtGeom();

    virtual void SelectCell(int SupLayer, int CelLayer, int Cell);
    virtual void SelectCellDet(unsigned long det);// Det Method
    virtual void Global2Local(double *global, double *local);
    virtual void Local2Global(double *local, double *global);
    virtual void WirePosAtEndcap(float *right, float *left);
    virtual void WirePosAtZ(float z, float *pos);
    virtual float GetWireAlfa();
    virtual float GetWireEpsilon();
    virtual float GetCellRad();
    virtual float GetCellInsideRad();

    virtual double DistFromWire(double *global);
    virtual double DistFromWireCenter(double *global);

    const TGeoHMatrix *GetGeoMatrix() { return _iTwire->GetGeoMatrix(); }

  private:

    ITrackerWireposition *_iTwire;

    Double_t tmpGlobal[3];
    Double_t tmpLocal[3];
    Float_t  tmpRight[3];
    Float_t  tmpLeft[3];
    Float_t  tmpPosAtZ[3];

};

}

#endif /* ITrackerGeom_CellGeometryHandle_ExtGeom_hh */
