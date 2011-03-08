#include "ITrackerGeom/inc/CellGeometryHandle_v2.hh"

namespace mu2e {

CellGeometryHandle_v2::CellGeometryHandle_v2(ITracker *itr) {
        _itr=itr;
}

CellGeometryHandle_v2::~CellGeometryHandle_v2() {
}

void CellGeometryHandle_v2::SelectCell(int SupLayer, int CelLayer, int Cell) {
        SuperLayer *sl=_itr->getSuperLayer(SupLayer);
        int SelectedCellLayer = CelLayer*2;
        if (CelLayer==0) SelectedCellLayer+=2;
        if (SupLayer==0) SelectedCellLayer++;
        /*boost::shared_ptr<ITLayer> */_itl=sl->getLayer(SelectedCellLayer);
        _cell=_itl->getCell(Cell);
        _cell->_tmpMidPoint  = _cell->_senseWire.get()->getMidPoint();
        _cell->_tmpDirection = _cell->_senseWire.get()->getDirection();
 //        _cell = _itr->getSuperLayer(SupLayer)->getLayer(CelLayer)->getCell(Cell);
        _matrx = _cell->getWire()->get3DTransfrom();
        _invmatrx = _cell->getWire()->get3DInvTransfrom();
        _fSuperLayer=SupLayer;
        _fLayer=CelLayer;
        _fWire=Cell;
}

} // namespace mu2e
