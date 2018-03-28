//--------------------------------------------------------------------------
// Name:
//   Mu2eDetectorModel: top-level detector model for Mu2e 
//
//        Copyright (C) 2015    Lawrence Berkeley Laboratory
// Author List:
//      Dave Brown 23 Dec 2015
//------------------------------------------------------------------------
#include "Mu2eBTrk/inc/Mu2eDetectorModel.hh"
#include "BTrk/DetectorModel/DetMaterial.hh"
#include "BTrk/MatEnv/MatDBInfo.hh"
#include "cetlib/coded_exception.h"
#include "TTrackerGeom/inc/TTracker.hh"

using namespace std;
namespace mu2e {
 
  Mu2eDetectorModel::Mu2eDetectorModel(fhicl::ParameterSet const& pset, TTracker const& ttracker) :
    _strawtype(0), 
    _gasmatname(pset.get<string>("StrawGasMaterialName","straw-gas")),
    _wallmatname(pset.get<string>("StrawWallMaterialName","straw-wall")),
    _wirematname(pset.get<string>("StrawWireMaterialName","straw-wire"))
    {
// find the materials needed for the tracker elements
/// Material information, BaBar style
    static MatDBInfo mat;
    const DetMaterial* gasmat = mat.findDetMaterial(_gasmatname.c_str());
    if(gasmat == 0)
      throw cet::exception("RECO")<<"mu2e::Mu2eDetectorModel: no material with name " 
      << _gasmatname << std::endl;
    const DetMaterial* wallmat = mat.findDetMaterial(_wallmatname.c_str());
    if(wallmat == 0)
      throw cet::exception("RECO")<<"mu2e::Mu2eDetectorModel: no material with name " 
      << _wallmatname << std::endl;
    // for now the wire type isn't set or used: FIXME!!!
    const DetMaterial* wiremat(0);
// overwrite the scattering fraction.  Sadly I must cast-off const for this
    double scatfrac = pset.get<double>("DahlLynchScatteringFraction",0.995);
    const_cast<DetMaterial*>(gasmat)->setScatterFraction(scatfrac);
    const_cast<DetMaterial*>(wallmat)->setScatterFraction(scatfrac);
// parameters for elements
    double tol = pset.get<double>("IntersectionTolerance",0.001);
// the offset displaces the element from the wire, which avoids problems when
// computing POCA on the fit trajectory.
    double offset = pset.get<double>("StrawElementOffset",0.25);
// maximum (fractional) radius to allow an intersection.  This avoids creating
// intersections with crazy long paths through the straw wall
    double rfrac = pset.get<double>("MaximumIntersectionRadiusFraction",0.96);
// construct the type.  This is reused by all straws
    _strawtype = new DetStrawType(gasmat,wallmat,wiremat,offset,tol,rfrac);
// loop over Planes
    // for(auto plane : ttracker.getPlanes()){ // tmp till we process cd3
    for ( size_t i=0; i!= ttracker.nPlanes(); ++i){
      const auto& plane = ttracker.getPlane(i);
  // loop over panels
      for(auto& panel : plane.getPanels()){
  // loop over straws
        for (const auto& straw : panel.getStrawPointers()) {
  // build the straw elements from this
          DetStrawElem* elem = new DetStrawElem(_strawtype,straw);
// push this into the map
          _strawmap[straw->id()] = elem;
	}
      }
    }
  }

  const DetStrawElem* Mu2eDetectorModel::strawElem(Straw const& straw) const{
    return strawElem(straw.id());
  }
  
  const DetStrawElem* Mu2eDetectorModel::strawElem(StrawId const& istraw) const{
    const DetStrawElem* retval(0);
    auto ifnd = _strawmap.find(istraw);
    if(ifnd != _strawmap.end())
      retval = ifnd->second;
    else
      throw cet::exception("RECO")<<"mu2e::Mu2eDetectorModel: no element associated to straw " 
      << istraw << std::endl;
    return retval;
  }

  Mu2eDetectorModel::~Mu2eDetectorModel() {
    delete _strawtype;
    for(auto istraw : _strawmap) {
      delete istraw.second;
    }
  }

}
