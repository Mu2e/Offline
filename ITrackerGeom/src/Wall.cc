// Framework includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Mu2e includes
#include "ITrackerGeom/inc/Wall.hh"

#ifndef __CINT__ 

namespace mu2e {

Wall::Wall(){
	_type			=	Wall::undefined;
	_pRmin			=	0.0;
	_pRmax			=	0.0;
	_pSPhi			=	0.0;
	_pDPhi			=	0.0;
	_pSTheta		=	0.0;
	_pDTheta		=	0.0;
	_pDz			=	0.0;
	_name			=	"";
//	_c				=	CLHEP::Hep3Vector(0.0,0.0,0.0);
	_nShells		=	0;
	_totalThickness	=	0.0;
}

Wall::Wall(Walltype wt){
	_type			=	wt;
	_pRmin			=	0.0;
	_pRmax			=	0.0;
	_pSPhi			=	0.0;
	_pDPhi			=	0.0;
	_pSTheta		=	0.0;
	_pDTheta		=	0.0;
	_pDz			=	0.0;
	_name			=	"";
//	_c				=	CLHEP::Hep3Vector(0.0,0.0,0.0);
	_nShells		=	0;
	_totalThickness	=	0.0;
}

Wall::Wall(const Wall &wl) {
	if (this!=&wl) {
		_type=wl.getType();
		_nShells=wl.getNShells();
		_totalThickness=wl.getTotalThickness();
		_pRmin=wl.getRmin();
		_pRmax=wl.getRmax();
		_pSPhi=wl.getSPhi();
		_pDPhi=wl.getDPhi();
		_pSTheta=wl.getSTheta();
		_pDTheta=wl.getDTheta();
		_pDz=wl.getDz();
		_name=wl.getName();
//		_c=wl.getC();
		_pos=wl.getPos();
		_materialsName=wl.getMaterialsName();
		_thicknesses=wl.getThicknesses();
	}
}

void Wall::addMaterials(int &wShellNumber, std::vector<std::string> *wShellsMatName, std::vector<double> *wShellsThicknesses) throw(cms::Exception) {
	_nShells=wShellNumber;
	if ( _nShells!=wShellsMatName->size() && _nShells!=wShellsThicknesses->size() )
		throw cms::Exception("GEOM")<< "Error in Configuration file! There is a disagreement between the vectors dimensions of a ITracker wall.\n";

	_materialsName.reset(wShellsMatName);
	_thicknesses.reset(wShellsThicknesses);
	_totalThickness = 0.0;
	for (int is = 0; is < _nShells; ++is) {
		_totalThickness += wShellsThicknesses->at(is);
	}
}

Wall& Wall::operator=(const Wall &wl) {
	if (this!=&wl) {
		_type=wl.getType();
		_nShells=wl.getNShells();
		_totalThickness=wl.getTotalThickness();
		_pRmin=wl.getRmin();
		_pRmax=wl.getRmax();
		_pSPhi=wl.getSPhi();
		_pDPhi=wl.getDPhi();
		_pSTheta=wl.getSTheta();
		_pDTheta=wl.getDTheta();
		_pDz=wl.getDz();
		_name=wl.getName();
//		_c=wl.getC();
		_pos=wl.getPos();
		_materialsName=wl.getMaterialsName();
		_thicknesses=wl.getThicknesses();
	}
}

} // namespace mu2e

#endif
