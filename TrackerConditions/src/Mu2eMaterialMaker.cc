
#include <iostream>

#include "Offline/TrackerConditions/inc/Mu2eMaterial.hh"
#include "Offline/TrackerConditions/inc/Mu2eMaterialMaker.hh"
#include "Offline/BTrkLegacy/inc/ExternalInfo.hh"


using namespace std;

namespace mu2e {


  Mu2eMaterial::ptr_t Mu2eMaterialMaker::fromFcl() {


    // make shared_ptr to the Mu2eDetector object on the heap
    Mu2eMaterial::ptr_t ptr = make_shared<Mu2eMaterial>();

    // file finder tells BTrk where to get data files
    ptr->_fileFinder = make_unique<FileFinder>(_config.elements(),
        _config.isotopes(),_config.materials());

    // particle info tells BTrk particle definitions
    // from our ParticleDataList
    ptr->_particleInfo = make_unique<ParticleInfo>();

    // this points BTrk static variables to our
    // implementations of file list and particle definitions
    ExternalInfo::set( ptr->_fileFinder.get() );
    ExternalInfo::set( ptr->_particleInfo.get() );

    // now construct the BTrk version of the description
    // of straw material

    string gasmatname = _config.strawGasMaterialName();
    string wallmatname = _config.strawWallMaterialName();
    string wirematname = _config.strawWireMaterialName();

    // MatDBInfo holds pointers to material in BTrk,
    // so needs to be long-lived
    MatDBInfo& mat = ptr->_mat;

    // this forces a pointer into the MatDBInfo cache,
    // the pointer also exists inside Btrk singletons
    const DetMaterial* gasmat = mat.findDetMaterial(gasmatname.c_str());
    if(gasmat == 0) {
      throw cet::exception("RECO_NO_GAS_MATERIAL")
        << "mu2e::Mu2eDetector: no material with name "
        << gasmatname << std::endl;
    }

    const DetMaterial* wallmat = mat.findDetMaterial(wallmatname.c_str());
    if(wallmat == 0) {
      throw cet::exception("RECO_NO_WALL_MATERIAL")
        <<"mu2e::Mu2eDetector: no material with name "
        << wallmatname << std::endl;
    }

    // for now the wire type isn't set or used: FIXME!!!
    const DetMaterial* wiremat(nullptr);

    // overwrite the scattering fraction.  Sadly I must cast-off const for this
    double scatfrac = _config.dahlLynchScatteringFraction();
    const_cast<DetMaterial*>(gasmat)->setScatterFraction(scatfrac);
    const_cast<DetMaterial*>(wallmat)->setScatterFraction(scatfrac);

    // the offset displaces the element from the wire, which avoids
    // problems when computing POCA on the fit trajectory.
    double offset = _config.strawElementOffset();

    // parameters for elements
    double tol = _config.intersectionTolerance();

    // maximum (fractional) radius to allow an intersection.  This avoids
    // creating intersections with crazy long paths through the straw wall
    double rfrac = _config.maximumIntersectionRadiusFraction();

    // construct the type.  This is reused by all straws
    // memory owned by Mu2eMaterial
    ptr->_strawtype = std::make_unique<DetStrawType>(
        gasmat,wallmat,wiremat,offset,tol,rfrac);

    return ptr;
  }


} // namespace mu2e
