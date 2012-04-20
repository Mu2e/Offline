#include "ITrackerGeom/inc/CellGeometryHandle_v2.hh"

namespace mu2e {

CellGeometryHandle_v2::CellGeometryHandle_v2(ITracker *itr) {
        _itr=itr;
}

CellGeometryHandle_v2::~CellGeometryHandle_v2() {
}

void CellGeometryHandle_v2::SelectCell() {
        SuperLayer *sl=_itr->getSuperLayer(_fSuperLayer);
        int SelectedCellLayer = _fLayer*2;
        if (_fLayer==0) SelectedCellLayer+=2;
        if (_fSuperLayer==0) SelectedCellLayer++;
        _itl=sl->getLayer(SelectedCellLayer);
        _cell=_itl->getCell(_fWire);
        _cell->_tmpMidPoint  = _cell->_senseWire.get()->getMidPoint();
        _cell->_tmpDirection = _cell->_senseWire.get()->getDirection();
 //        _cell = _itr->getSuperLayer(_fSuperLayer)->getLayer(CelLayer)->getCell(Cell);
        _matrx = _cell->getWire()->get3DTransfrom();
        _invmatrx = _cell->getWire()->get3DInvTransfrom();
}

void CellGeometryHandle_v2::SelectCell(int SupLayer, int CelLayer, int Cell, bool isUpstrm) {
        _fSuperLayer=SupLayer;
        _fLayer=CelLayer;
        _fWire=Cell;
        _absRadID = SupLayer*_itr->nRing()+CelLayer;
        SelectCell();
}

void  CellGeometryHandle_v2::SelectCell(int absRadID, int Cell, bool isUpstrm) {
        _absRadID=absRadID;
        _fSuperLayer=absRadID/_itr->nRing();
        _fLayer=absRadID-_itr->nRing()*_fSuperLayer;
        _fWire=Cell;
        SelectCell();
}

int CellGeometryHandle_v2::computeAbsRadID(int SupLayer, int CelLayer, bool isUpstrm) {
        int absRadID = SupLayer*_itr->nRing()+CelLayer;
        return absRadID;
}

} // namespace mu2e
