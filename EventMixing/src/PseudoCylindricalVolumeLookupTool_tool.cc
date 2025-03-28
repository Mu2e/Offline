// Ed Callaghan
// Stop-gap tool to fake a volume lookup, only valid for protonabs1 (IPA)
// November 2024

#include "Offline/EventMixing/inc/PseudoCylindricalVolumeLookupTool.hh"

namespace mu2e{
  PseudoCylindricalVolumeLookupTool::PseudoCylindricalVolumeLookupTool(const Parameters& config):
      _ipa_name(config().ipa_name()),
      _st_name(config().st_name()),
      _other_name(config().other_name()),
      _initialized(false){
    /**/
  }

  PseudoCylindricalVolumeLookupTool::~PseudoCylindricalVolumeLookupTool() = default;

  std::string PseudoCylindricalVolumeLookupTool::Volume(const CLHEP::Hep3Vector& position_mu2e){
    // deferred initialization, to make use of GeometryService
    if (!_initialized){
      this->initialize();
    }

    // by default, position does not lie anywhere interesting
    std::string rv = _other_name;

    // translate position into det. coordinates, where IPA is centered radially
    const auto position = (*_frame)->toDetector(position_mu2e);

    // ipa
    // manual inspection of cylindrical geo. object, in det. coordinates
    const auto& part = (*_ipa)->part(0);
    auto center_mu2e = part.center();
    auto center = (*_frame)->toDetector(center_mu2e);
    auto znom = center.z();
    auto zlo = znom - part.halfLength();
    auto zhi = znom + part.halfLength();
    auto rlo = std::min(part.innerRadiusAtStart(),
                        part.innerRadiusAtEnd());
    auto rhi = std::max(part.outerRadiusAtStart(),
                        part.outerRadiusAtEnd());

    if (this->in_cylindrical_shell(position, zlo, zhi, rlo, rhi)){
      rv = _ipa_name;
      return rv;
    }

    // stopping target
    for (int i = 0 ; i < (*_st)->nFoils() ; i++){
      const auto& foil = (*_st)->foil(i);
      center = foil.centerInDetectorSystem();
      znom = center.z();
      zlo = znom - foil.halfThickness();
      zhi = znom + foil.halfThickness();
      rlo = foil.rIn();
      rhi = foil.rOut();
      if (this->in_cylindrical_shell(position, zlo, zhi, rlo, rhi)){
        rv = _st_name;
        return rv;
      }
    }

    return rv;
  }

  std::string PseudoCylindricalVolumeLookupTool::StartVolume(const SimParticle& particle){
    // deferred initialization, to make use of GeometryService
    if (!_initialized){
      this->initialize();
    }

    // query volume of starting position
    const auto& position = particle.startPosition();
    auto rv = this->Volume(position);
    return rv;
  }

  // deferred initialization is necessary to make use of GeometryService
  void PseudoCylindricalVolumeLookupTool::initialize(){
    _frame = std::make_unique< GeomHandle<DetectorSystem> >();
    _ipa = std::make_unique< GeomHandle<ProtonAbsorber> >();
    _st = std::make_unique< GeomHandle<StoppingTarget> >();
    _initialized = true;
  }

  bool PseudoCylindricalVolumeLookupTool::in_cylindrical_shell(
                                              const CLHEP::Hep3Vector& position,
                                              double zlo, double zhi,
                                              double rlo, double rhi){
    // test whether positions lies within cylidrical shell
    // nested to defer the more-expensive radius calculation
    bool rv = false;
    const auto z = position.z();
    if ((zlo < z) && (z < zhi)){
      const auto r = position.perp();
      if ((rlo < r) && (r < rhi)){
        rv = true;
      }
    }
    return rv;
  }
} // namespace mu2e

DEFINE_ART_CLASS_TOOL(mu2e::PseudoCylindricalVolumeLookupTool)
