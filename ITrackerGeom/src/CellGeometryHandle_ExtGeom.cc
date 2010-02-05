#include "ITrackerGeom/inc/CellGeometryHandle_ExtGeom.hh"


namespace mu2e {

CellGeometryHandle_ExtGeom::CellGeometryHandle_ExtGeom(const char *WireDataFile) {
	_iTwire = new ITrackerWireposition(false,WireDataFile);

}

CellGeometryHandle_ExtGeom::~CellGeometryHandle_ExtGeom() {
	delete _iTwire;
}

void CellGeometryHandle_ExtGeom::SelectCell(int SupLayer, int CelLayer, int Cell){
	_iTwire->SelectWire(SupLayer,CelLayer,Cell);
	_fSuperLayer=_iTwire->GetSuperLayer();
	_fLayer=_iTwire->GetCelRing();
	_fWire=_iTwire->GetWire();
}

void CellGeometryHandle_ExtGeom::SelectWireDet(unsigned long  det)
{
	_iTwire->SelectWireDet(det);
	_fSuperLayer=_iTwire->GetSuperLayer();
	_fLayer=_iTwire->GetCelRing();
	_fWire=_iTwire->GetWire();
}

void CellGeometryHandle_ExtGeom::Global2Local(double *global, double *local)
{
	tmpGlobal[0]=global[0];
	tmpGlobal[1]=global[1];
	tmpGlobal[2]=global[2];
	_iTwire->Global2Local(tmpGlobal,tmpLocal);
	local[0]=tmpLocal[0];
	local[1]=tmpLocal[1];
	local[2]=tmpLocal[2];
}

void CellGeometryHandle_ExtGeom::Local2Global(double *local, double *global)
{
	local[0]=tmpLocal[0];
	local[1]=tmpLocal[1];
	local[2]=tmpLocal[2];
	_iTwire->Local2Global(tmpLocal,tmpGlobal);
	tmpGlobal[0]=global[0];
	tmpGlobal[1]=global[1];
	tmpGlobal[2]=global[2];
}

void CellGeometryHandle_ExtGeom::WirePosAtEndcap(float *right, float *left)
{
	_iTwire->WirePosAtEndcap(tmpRight,tmpLeft);
	right[0]=tmpRight[0];
	right[1]=tmpRight[1];
	right[2]=tmpRight[2];
	left[0]=tmpLeft[0];
	left[1]=tmpLeft[1];
	left[2]=tmpLeft[2];
}

void CellGeometryHandle_ExtGeom::WirePosAtZ(float z, float *pos)
{
	_iTwire->WirePosAtZ(z,tmpPosAtZ);
	pos[0]=tmpPosAtZ[0];
	pos[1]=tmpPosAtZ[1];
	pos[2]=tmpPosAtZ[2];
}

float CellGeometryHandle_ExtGeom::GetWireAlfa()
{
	return _iTwire->GetWireAlfa();
}

float CellGeometryHandle_ExtGeom::GetWireEpsilon()
{
	return _iTwire->GetWireEpsilon();
}

float CellGeometryHandle_ExtGeom::GetCellRad()
{
	return _iTwire->GetLayerRad();
}


double CellGeometryHandle_ExtGeom::DistFromWireCenter(double *global)
{
	tmpGlobal[0]=global[0];
	tmpGlobal[1]=global[1];
	tmpGlobal[2]=global[2];
	return _iTwire->DistFromWireCenter(tmpGlobal);
}

double CellGeometryHandle_ExtGeom::DistFromWire(double *global)
{
	tmpGlobal[0]=global[0];
	tmpGlobal[1]=global[1];
	tmpGlobal[2]=global[2];
	return _iTwire->DistFromWire(tmpGlobal);
}

} // namespace mu2e
