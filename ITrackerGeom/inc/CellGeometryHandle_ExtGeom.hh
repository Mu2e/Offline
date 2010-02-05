#ifndef CELLGEOMETRYHANDLEEXTGEOM_HH
#define CELLGEOMETRYHANDLEEXTGEOM_HH

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
    virtual void SelectWireDet(unsigned long det);// Det Method
    virtual void Global2Local(double *global, double *local);
    virtual void Local2Global(double *local, double *global);
    virtual void WirePosAtEndcap(float *right, float *left);
    virtual void WirePosAtZ(float z, float *pos);
    virtual float GetWireAlfa();
    virtual float GetWireEpsilon();
    virtual float GetCellRad();

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

#endif /* CELLGEOMETRYHANDLEEXTGEOM_HH */
