// Andrei Gaponenko, 2016

#include <string>

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Utilities/InputTag.h"

#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "CutAndCountAnalysis/inc/CutAndCountAnalysis.hh"

namespace mu2e {
  class CutAndCountAnalyzer : public art::EDAnalyzer {
  public:

    explicit CutAndCountAnalyzer(const fhicl::ParameterSet& pset);
    virtual void analyze(const art::Event& event) override;

  private:
    CutAndCountAnalysis an_;
  };

  //================================================================
  CutAndCountAnalyzer::CutAndCountAnalyzer(const fhicl::ParameterSet& pset)
    : art::EDAnalyzer(pset)
    , an_(pset, *art::ServiceHandle<art::TFileService>())
  {}

  //================================================================
  void CutAndCountAnalyzer::analyze(const art::Event& event) {
    an_.accepted(event);
  }

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::CutAndCountAnalyzer);
