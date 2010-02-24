#include "ITrackerGeom/inc/CellGeometryHandle_v3.hh"

namespace mu2e {

CellGeometryHandle_v3::CellGeometryHandle_v3(ITracker *itr) {
	_itr = itr;
	_nLayer = _itr->nSuperLayers();
	_nLayer--;
}

CellGeometryHandle_v3::~CellGeometryHandle_v3() {
}

void CellGeometryHandle_v3::SelectCell(int SupLayer, int CelLayer, int Cell) {
	SuperLayer *sl=_itr->getSuperLayer(SupLayer);
	if (SupLayer==_nLayer && CelLayer==1) CelLayer=0;
	if (CelLayer>0)	throw cms::Exception("GEOM")<< "The requested cell layer is not allowed in the Square cell geometry"<<std::endl;
	int SelectedCellLayer = 1;
	if (SupLayer==0) SelectedCellLayer++;
	_itl=sl->getLayer(SelectedCellLayer);
	_cell=_itl->getCell(Cell);
	_matrx = _cell->getWire()->get3DTransfrom();
	_invmatrx = _cell->getWire()->get3DInvTransfrom();
	_fSuperLayer=SupLayer;
	_fLayer=CelLayer;
	_fWire=Cell;
}

void CellGeometryHandle_v3::SelectWireDet(unsigned long  det) {
	// Return the SuperLayer
	int fSuperLayer=(int)(det*0.0001);

	//Return the Wire
	int fWire=(int)((det)-((fSuperLayer)*10000));

	fSuperLayer--;

	//Call the upper method
	SelectCell(fSuperLayer,0,fWire);
}

unsigned long CellGeometryHandle_v3::computeDet(int SupLayer, int CelLayer, int Cell) {
	unsigned long det = ((SupLayer+1)*10000)+(Cell);
	return det;
}

} // namespace mu2e
