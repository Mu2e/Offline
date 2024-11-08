// Ed Callaghan
// Tool-based filter to selectively remove dts data products
// October 2024

// stl
#include <memory>
#include <vector>

// art
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Utilities/make_tool.h"

// cetlib_except
#include "cetlib_except/exception.h"

// fhiclcpp
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/OptionalDelegatedParameter.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Name.h"

// mu2e
#include "Offline/EventMixing/inc/DetectorStepSelectionTool.hh"
#include "Offline/EventMixing/inc/UniversalDetectorStepSelectionTool.hh"

namespace mu2e{
  class MixingFilter: public art::EDProducer{
    public:
      struct Config{
        fhicl::Atom<art::InputTag> tag{
          fhicl::Name("MixingModule"),
          fhicl::Comment("art::InputTag to locate data products to filter")
        };
        fhicl::OptionalDelegatedParameter selection{
          fhicl::Name("selection"),
          fhicl::Comment("Configuration for tool to select detector steps")
        };
      };

      using Parameters = art::EDProducer::Table<Config>;
      MixingFilter(const Parameters&);

      void produce(art::Event&);

    protected:
      // use active selection, as opposed to rejection, because we must
      // instantiate a new collection either way, so rejection e.g.
      // via std::remove_if would entail a second pass
      using Selection = DetectorStepSelectionTool;
      std::unique_ptr<Selection> _selection;
      art::InputTag _tag;

      // templated to reduce code surface
      template<typename T>
      void register_dependencies();
      template<typename T>
      void filter_into(const art::Event&, T&);
      template<typename T>
      void filter_into_steplike(const T&, T&);

    private:
      /**/
  };

  MixingFilter::MixingFilter(const Parameters& config):
      art::EDProducer(config),
      _tag(config().tag()){

    // by default, apply a trivial filter which removes nothing
    fhicl::ParameterSet selection_config;
    if (config().selection.hasValue()){
      auto table = config().selection;
      selection_config = table.get_if_present<fhicl::ParameterSet>().value();
    }
    else{
      selection_config.put("tool_type", "UniversalDetectorStepSelectionTool");
    }
    _selection = art::make_tool<Selection>(selection_config);

    // framework hooks
    this->register_dependencies<CaloShowerStepCollection>();
    this->register_dependencies<CrvStepCollection>();
    this->register_dependencies<StrawGasStepCollection>();
  }

  void MixingFilter::produce(art::Event& event){
    auto csss = std::make_unique<CaloShowerStepCollection>();
    this->filter_into(event, *csss);
    event.put(std::move(csss));

    auto crvs = std::make_unique<CrvStepCollection>();
    this->filter_into(event, *crvs);
    event.put(std::move(crvs));

    auto sgss = std::make_unique<StrawGasStepCollection>();
    this->filter_into(event, *sgss);
    event.put(std::move(sgss));
  }

  template<typename T>
  void MixingFilter::register_dependencies(){
    this->consumes<T>(_tag);
    this->produces<T>();
  }

  template<typename T>
  void MixingFilter::filter_into(const art::Event& event, T& out){
    auto handle = event.getValidHandle<T>(_tag);
    this->filter_into_steplike(*handle, out);
  }

  template<typename T>
  void MixingFilter::filter_into_steplike(const T& in, T& out){
    auto selection = [this] (const T::value_type& step) {
      bool rv = _selection->Select(step);
      return rv;
    };
    out.reserve(in.size());
    auto it = std::copy_if(in.begin(), in.end(), out.begin(), selection);
    size_t size = std::distance(out.begin(), it);
    out.resize(size);
  }
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::MixingFilter)
