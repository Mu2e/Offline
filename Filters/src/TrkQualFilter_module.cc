
// ======================================================================
//
// TrkQualFilter_module: allows filtering on the TrkQual 
// when given a cut efficiency
//
// ======================================================================

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Principal/Handle.h"

// Mu2e includes.
#include "RecoDataProducts/inc/TrkQual.hh"
#include "Mu2eUtilities/inc/MVATools.hh"
#include "DbService/inc/DbHandle.hh"
#include "DbTables/inc/TrkQualDb.hh"


#include <iostream>
#include <string>

using namespace std;
namespace mu2e {

  class TrkQualFilter : public art::EDFilter {
  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Atom<art::InputTag> inputTag{Name("inputTag"), Comment("Input tag for TrkQualCollection")};
      fhicl::Atom<float> effRequest{Name("effRequest"), Comment("Cut efficiency requested")};
    };

    using Parameters = art::EDFilter::Table<Config>;
    explicit TrkQualFilter(const Parameters& conf);

  private:
    bool beginRun(art::Run& run) override;
    bool endRun(art::Run& run) override;
    bool filter(art::Event& event) override;

    mu2e::DbHandle<mu2e::TrkQualDb> _trkQualDb;

    art::InputTag _inputTag;
    float _effRequest;
  };

  TrkQualFilter::TrkQualFilter(const Parameters& conf):
    art::EDFilter{conf},
    _inputTag(conf().inputTag()),
    _effRequest(conf().effRequest())
    { 

    }

  bool TrkQualFilter::beginRun(art::Run& run) {
    return true;
  }
  bool TrkQualFilter::endRun(art::Run& run) {
    return true;
  }

  bool TrkQualFilter::filter(art::Event& event) {

    float trkQualCut = -1;
    auto const& trkQualTable = _trkQualDb.get(event.id());
    for (const auto& i_row : trkQualTable.rows()) {
      if (i_row.mvaname() == "TrkQual") {
	const auto& calib = trkQualTable.getCalib(i_row.idx());
	for (const auto& i_pair : calib) {
	  if (i_pair.first == _effRequest) {
	    trkQualCut = i_pair.second;
	    break;
	  }
	}
      }
    }

    if (trkQualCut < 0) {
      throw cet::exception("TrkQualFilter") << "trkQualCut is less than 0 (value = " << trkQualCut << ")" << std::endl;
    }
    
    auto trkQualCollsH = event.getValidHandle<TrkQualCollection>(_inputTag);

    bool pass = false;
    for(const auto& i_trkQual : *trkQualCollsH) {
      std::cout << "TrkQual = " << i_trkQual.MVAOutput() << std::endl;
      
      if (i_trkQual.MVAOutput() >= trkQualCut) {
      	pass = true;
      }
    }

    std::cout << "Passed? " << pass << std::endl;

    if (pass) {
      return true;
    }
    return false;
  }
}

using mu2e::TrkQualFilter;
DEFINE_ART_MODULE(TrkQualFilter);
