
//
// This module looks for early tracker digis (before the flash blanking)
// and uses them to estimate the proton bunch time of this event WRT
// the DAQ clock.  This is used mainly as a diagnostic.
//
// Original author David Brown, 10/1/2020 LBNL
//
// framework
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/EDProducer.h"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/RecoDataProducts/inc/ProtonBunchTime.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "TH1F.h"
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/error_of.hpp>
#include <boost/accumulators/statistics/error_of_mean.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>

#include <iostream>
#include <limits>

namespace mu2e {
  using namespace TrkTypes;

  class ProtonBunchTimeFromStrawDigis : public art::EDProducer
  {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      struct Config {
        fhicl::Atom<int> debug{ Name("DebugLevel"), Comment("Debug level"), 0};
        fhicl::Atom<int> diag{ Name("DiagLevel"), Comment("Diag level"), 0};
        fhicl::Atom<int> print{ Name("PrintLevel"), Comment("Print level"), 0};
        fhicl::Atom<unsigned> nbins{ Name("NBins"), Comment("Number of bins in time histogram"), 30};
        fhicl::Atom<unsigned> minnhits{ Name("MinNHits"), Comment("Minimum number of hits needed to use result"), 100};
        fhicl::Atom<float> maxtime { Name("MaxDigiTime"), Comment("Maximum time to accept the digi as 'early' (nsec)"),200.0 };
        fhicl::Atom<float> meantime { Name("MeanTimeOffset"), Comment("Mean time offset of 'early' hits (nsec)"),99.5}; // this should come from a database FIXME!
        fhicl::Atom<float> phasewindow { Name("PhaseWindow"), Comment("Width of Event Marker Phase Window (nsec)"),25.0};
        fhicl::Atom<art::InputTag> strawDigiTag{ Name("StrawDigiTag"), Comment("StrawDigi producer"),"makeSD" };
      };
      using Parameters = art::EDProducer::Table<Config>;
      explicit ProtonBunchTimeFromStrawDigis(Parameters const& config);
      virtual ~ProtonBunchTimeFromStrawDigis();
      void produce( art::Event& e) override;
      void beginJob ( ) override;
    private:
      int debug_, diag_, print_;
      unsigned nbins_, minnhits_;
      float maxtime_, meantime_;
      float phasemean_, phaseerr_;
      art::InputTag ewmtag_, sdtag_;
      ProditionsHandle<StrawResponse> _strawResponse_h;
  };

  ProtonBunchTimeFromStrawDigis::ProtonBunchTimeFromStrawDigis(Parameters const& config) :
    EDProducer(config),
    debug_(config().debug()),
    diag_(config().diag()),
    print_(config().print()),
    nbins_(config().nbins()),
    minnhits_(config().minnhits()),
    maxtime_(config().maxtime()),
    meantime_(config().meantime()),
    sdtag_(config().strawDigiTag())
  {
    float phasewindow = config().phasewindow();
    phasemean_ = phasewindow/2.0;
    phaseerr_ = phasewindow/sqrt(12.0);
    consumes<StrawDigiCollection>(sdtag_);
    produces<ProtonBunchTime>();
  }

  ProtonBunchTimeFromStrawDigis::~ProtonBunchTimeFromStrawDigis() {}

  void ProtonBunchTimeFromStrawDigis::beginJob(){
  }

  void ProtonBunchTimeFromStrawDigis::produce(art::Event& event) {
    using namespace boost::accumulators;
    auto const& srep = _strawResponse_h.get(event.id());
    auto sdH = event.getValidHandle<StrawDigiCollection>(sdtag_);
    auto const& sdcol(*sdH);
    std::unique_ptr<ProtonBunchTime> pbto(new ProtonBunchTime());
    // statistics accumulators
    //
    accumulator_set<float, stats<tag::mean, tag::error_of<tag::mean> >> meanacc;

    TH1F* timeplot(0);
    if(debug_ > 1 && event.id().event() < 1000) {
      art::ServiceHandle<art::TFileService> tfs;
      char hname[100];
      snprintf(hname,100,"PBT%i",event.id().event());
      timeplot = tfs->make<TH1F>(hname,"time spectrum;nsec",nbins_,0.0,maxtime_);
    }

    for (size_t isd=0;isd<sdcol.size();++isd) {
      const StrawDigi& digi = sdcol[isd];
      // correct to physical time
      TDCTimes times;
      srep.calibrateTimes(digi.TDC(),times,digi.strawId());
      // take the early time
      float mintime = std::min(times[StrawEnd::hv],times[StrawEnd::cal]);
      if(mintime < maxtime_){
        // This can be improved by using TOT to approximately correct for the propagation time, and maybe peak/ped to correct
        // for time slewing TODO!
        meanacc(mintime);
        // early hit: accumulate the histogram
        if(print_ > 1) std::cout << "Early hit time " << mintime << " SID " << digi.strawId() << std::endl;
        if(timeplot)timeplot->Fill(mintime);
      }
    }
    // if there aren't enought hits, set to 0
    unsigned nhits = extract_result<tag::count>(meanacc);
    if(nhits < minnhits_){
      pbto->pbtime_ = phasemean_;
      pbto->pbterr_ = phaseerr_;
    } else {
      pbto->pbtime_ = extract_result<tag::mean>(meanacc) - meantime_;
      pbto->pbterr_ = extract_result<tag::error_of<tag::mean>>(meanacc);
    }
    event.put(std::move(pbto));
  }
}

using mu2e::ProtonBunchTimeFromStrawDigis;
DEFINE_ART_MODULE(ProtonBunchTimeFromStrawDigis)

