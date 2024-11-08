// Ed Callaghan
// Stop-gap tool to fake a volume lookup, only valid for protonabs1 (IPA)
// November 2024

#include "Offline/EventMixing/inc/InnerProtonAbsorberPseudoVolumeLookupTool.hh"

namespace mu2e{
  InnerProtonAbsorberPseudoVolumeLookupTool::InnerProtonAbsorberPseudoVolumeLookupTool(const Parameters& config):
      _ipa_name(config().ipa_name()),
      _other_name(config().other_name()),
      _initialized(false){
    /**/
  }

  InnerProtonAbsorberPseudoVolumeLookupTool::~InnerProtonAbsorberPseudoVolumeLookupTool() = default;

  std::string InnerProtonAbsorberPseudoVolumeLookupTool::Volume(const CLHEP::Hep3Vector& position_mu2e){
    // deferred initialization, to make use of GeometryService
    if (!_initialized){
      this->initialize();
    }

    // be default, position does not lie in IPA
    std::string rv = _other_name;

    // translate position into det. coordinates, where IPA is centered radially
    const auto position = (*_frame)->toDetector(position_mu2e);

    // manual inspection of cylindrical geo. object, in det. coordinates
    const auto& part = (*_ipa)->part(0);
    const auto center_mu2e = part.center();
    const auto center = (*_frame)->toDetector(center_mu2e);
    const auto znom = center.z();
    const auto zlo = znom - part.halfLength();
    const auto zhi = znom + part.halfLength();

    // test whether positions lies within cylidrical shell
    // nested to defer the more-expensive radius calculation
    const auto z = position.z();
    if ((zlo < z) && (z < zhi)){
      const auto r = position.perp();
      const auto rlo = std::min(part.innerRadiusAtStart(),
                                part.innerRadiusAtEnd());
      const auto rhi = std::max(part.outerRadiusAtStart(),
                                part.outerRadiusAtEnd());
      if ((rlo < r) && (r < rhi)){
        rv = _ipa_name;
      }
    }

    return rv;
  }

  std::string InnerProtonAbsorberPseudoVolumeLookupTool::StartVolume(const SimParticle& particle){
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
  void InnerProtonAbsorberPseudoVolumeLookupTool::initialize(){
    _frame = std::make_unique< GeomHandle<DetectorSystem> >();
    _ipa = std::make_unique< GeomHandle<MECOStyleProtonAbsorber> >();
    _initialized = true;
  }
} // namespace mu2e

DEFINE_ART_CLASS_TOOL(mu2e::InnerProtonAbsorberPseudoVolumeLookupTool)
