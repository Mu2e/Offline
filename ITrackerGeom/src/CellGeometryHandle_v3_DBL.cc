#include "ITrackerGeom/inc/CellGeometryHandle_v3_DBL.hh"

namespace mu2e {

CellGeometryHandle_v3_DBL::CellGeometryHandle_v3_DBL(ITracker *itr):CellGeometryHandle_v3(itr){
        _nSLayer=itr->nSuperLayers();
}


CellGeometryHandle_v3_DBL::~CellGeometryHandle_v3_DBL() {
}

inline void  CellGeometryHandle_v3_DBL::AdjustMatrix() {
        if (_isDownStream) {
                _DnStrmDeadWireLngt=_itr->zZonesLimits()[1]/cos(_cell->getWire()->getEpsilon());
                _tmpCellHalfLength=_cell->getHalfLength()-_DnStrmDeadWireLngt;
                _tmpCellHalfLength*=0.5;
                WirePosAtLength(_tmpCellHalfLength+_DnStrmDeadWireLngt,wCntPos);
                _tmpMidPoint.set(wCntPos[0],wCntPos[1],wCntPos[2]);
                _tmpDirection.set(_cell->getDirection().x(),_cell->getDirection().y(),_cell->getDirection().z());
        }
        else if (_isUpStream) {
                _UpStrmDeadWireLngt=_itr->zZonesLimits()[0]/cos(_cell->getWire()->getEpsilon());
                _tmpCellHalfLength=_cell->getHalfLength()+_UpStrmDeadWireLngt;
                _tmpCellHalfLength*=0.5;
                WirePosAtLength(-_tmpCellHalfLength+_UpStrmDeadWireLngt,wCntPos);
                _tmpMidPoint.set(wCntPos[0],wCntPos[1],wCntPos[2]);
                _tmpDirection.set(-_cell->getDirection().x(),-_cell->getDirection().y(),-_cell->getDirection().z());
        }
}

void CellGeometryHandle_v3_DBL::SelectCell(int SupLayer, int CelLayer, int Cell, bool isUpstrm) {
        if (SupLayer>=_nSLayer || isUpstrm){
                _isUpStream=true;
                _isDownStream=false;
                SupLayer=SupLayer%_nSLayer;
        }
        else {
                _isUpStream=false;
                _isDownStream=true;
        }
        CellGeometryHandle_v3::SelectCell(SupLayer,CelLayer,Cell);
        AdjustMatrix();
}

void CellGeometryHandle_v3_DBL::SelectCellDet(unsigned long  det) {
        // Return the SuperLayer
        int fSuperLayer=(int)(det*0.0001);

        //Return the Wire
        int fWire=(int)((det)-((fSuperLayer)*10000));

        fSuperLayer--;
        if (fSuperLayer>=_nSLayer){
                _isUpStream=true;
                _isDownStream=false;
                fSuperLayer=fSuperLayer%_nSLayer;
        }
        else {
                _isUpStream=false;
                _isDownStream=true;
        }

        //Call the upper method
        CellGeometryHandle_v3::SelectCell(fSuperLayer,0,fWire);
        AdjustMatrix();
}

void  CellGeometryHandle_v3_DBL::SelectCell(int absRadID, int Cell, bool isUpstrm) {
        SelectCell(absRadID,0,Cell,isUpstrm);
}

unsigned long CellGeometryHandle_v3_DBL::computeDet(int SupLayer, int CelLayer, int Cell, bool isUpstrm) {
        if (SupLayer<_nSLayer && isUpstrm) {
                SupLayer+=_nSLayer;
        }
        return CellGeometryHandle_v3::computeDet(SupLayer,CelLayer,Cell);
}

int CellGeometryHandle_v3_DBL::computeAbsRadID(int SupLayer, int CelLayer, bool isUpstrm) {
        if (SupLayer>=_nSLayer || isUpstrm){
                SupLayer=SupLayer%_nSLayer;
        }
        int absRadID = SupLayer;
        return absRadID;
}

const CLHEP::Hep3Vector& CellGeometryHandle_v3_DBL::GetCellCenter() const {
        return _tmpMidPoint;
}

const CLHEP::Hep3Vector& CellGeometryHandle_v3_DBL::GetCellDirection() const {
        return _tmpDirection;
}

double CellGeometryHandle_v3_DBL::GetCellHalfLength() const {
        return _tmpCellHalfLength;
}

} // namespace mu2e
