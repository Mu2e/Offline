//
// Create zero-suppressed STMDigis from unsuppressed STMDigis
//
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"

#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include <utility>
#include <algorithm>
// root
#include "TH1F.h"
#include "TF1.h"

#include "Offline/DataProducts/inc/STMTypes.hh"
#include "Offline/RecoDataProducts/inc/STMDigi.hh"
#include "Offline/RecoDataProducts/inc/STMHit.hh"
#include "Offline/STMReco/inc/MWDAlg.hh"

using namespace std;
using CLHEP::Hep3Vector;
namespace mu2e {

  class STMMovingWindowDeconvolution : public art::EDProducer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<art::InputTag> stmDigisTag{ Name("stmDigisTag"), Comment("InputTag for STMDigiCollection")};
        fhicl::Atom<double> M{ Name("M"), Comment("Input tag for number of channels")};
        fhicl::Atom<double> L{ Name("L"), Comment("Input tage for the L parameter")};
        fhicl::Atom<double> tau{ Name("tau"), Comment("Input tag for the RC time constant")};
        fhicl::Atom<double> nsigma_cut{ Name("nsigma_cut"), Comment("Input tag for the the nsigma cut parameter")};
        fhicl::Atom<double> thresholdgrad{ Name("thresholdgrad"), Comment("Input tag for the threshold gradient")};
        fhicl::Atom<double> fADC{ Name("fADC"), Comment("Input tag for the ADC frequency")};
        fhicl::Atom<int> cut_mode{ Name("cut_mode"), Comment("Input tag for the cut mode")};
        fhicl::Atom<double> fixed_cut_parameter{ Name("fixed_cut_parameter"), Comment("Input tag for the fixed cut parameter")};
      };
      using Parameters = art::EDProducer::Table<Config>;
      explicit STMMovingWindowDeconvolution(const Parameters& conf);

    private:
    void beginJob() override;
      void produce(art::Event& e) override;

    art::InputTag _stmDigisTag;
    mu2e::MWDAlg _mwd;
  };

  STMMovingWindowDeconvolution::STMMovingWindowDeconvolution(const Parameters& config )  :
    art::EDProducer{config}
    ,_stmDigisTag(config().stmDigisTag())
    ,_mwd(config().M(),config().L(),config().tau(),config().nsigma_cut(),config().thresholdgrad(),config().fADC(),config().cut_mode(),config().fixed_cut_parameter())
  {
    consumes<STMDigiCollection>(_stmDigisTag);
    produces<STMDigiCollection>();
  }

  void STMMovingWindowDeconvolution::beginJob() {
  }
    void STMMovingWindowDeconvolution::produce(art::Event& event) {
    // create output
    unique_ptr<STMDigiCollection> outputSTMDigis(new STMDigiCollection);
    auto digisHandle = event.getValidHandle<STMDigiCollection>(_stmDigisTag);

    //    std::cout << _mwd.print() << std::endl;
    for (const auto& digi : *digisHandle) {
      mu2e::data dat;
      dat.adc = &digi.adcs()[0];
      dat.t0 = digi.trigTime();
      dat.nadc = digi.adcs().size();
      _mwd.mwd_algorithm(&dat);
      //      _mwd.mwd_algorithm(&digi.adcs()[0], digi.adcs().size());
      std::vector<double> baseline =  _mwd.calculate_baseline();
            std::cout << "AE: baseline = " << baseline.at(0) << ", sigma = " << baseline.at(1) << std::endl;

      auto peaks = _mwd.find_peaks(baseline.at(0), baseline.at(1), 0);

      for (int i_peak = 0; i_peak < peaks->npeaks; ++i_peak) {
        std::vector<int16_t> mwd_adcs;
        mwd_adcs.push_back(peaks->peak_heights.at(i_peak));

        uint32_t extra = ((uint32_t)baseline.at(1) << 16) | ((uint32_t)baseline.at(0)); // RMS << Mean
        STMDigi stm_digi(digi.trigNum(), STMTrigType(digi.trigType().mode(), digi.trigType().channel(), STMDataType::kMWD), digi.trigTime()+(peaks->peak_times.at(i_peak))*1e3, 0, extra, STMDigiFlag::kOK, mwd_adcs);
        outputSTMDigis->push_back(stm_digi);
      }
    }
    std::cout << "AE: No. MWD Digis = " << outputSTMDigis->size() << std::endl;
    event.put(std::move(outputSTMDigis));
  }
}

DEFINE_ART_MODULE(mu2e::STMMovingWindowDeconvolution)
