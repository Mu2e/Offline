//
// A module to convert simulated CRV (Wideband) hits into a ROOT file
// with a structure similar to the one used for the Wideband standalone analysis.
//
//
// Original Author: Ralf Ehrlich

#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "Offline/DataProducts/inc/CRVId.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/CRVConditions/inc/CRVOrdinal.hh"

#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/MCDataProducts/inc/CrvStep.hh"
#include "Offline/MCDataProducts/inc/CrvCoincidenceClusterMC.hh"
#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"
#include "Offline/RecoDataProducts/inc/CrvCoincidenceCluster.hh"
#include "Offline/RecoDataProducts/inc/ProtonBunchTime.hh"

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"

#include <string>

#include <TMath.h>
#include <TNtuple.h>
#include <TH1F.h>
#include <TF1.h>
#include <TFitResult.h>

namespace mu2e
{
  class CrvTimingStudies : public art::EDAnalyzer
  {
    public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct Config
    {
      fhicl::Atom<std::string> crvRecoPulsesModuleLabel{ Name("crvRecoPulseModuleLabel"), Comment("CrvRecoPulse Label")};
      fhicl::Atom<double>      PEthreshold{ Name("PEthreshold"), Comment("PE threshold")};
    };

    using Parameters = art::EDAnalyzer::Table<Config>;
    explicit CrvTimingStudies(const Parameters& conf);
    ~CrvTimingStudies() override;
    void analyze(const art::Event& e) override;
    void beginJob() override;
    void endJob() override;

    private:
    std::string   _crvRecoPulsesModuleLabel;
    double        _PEthreshold;

    std::map<std::pair<int,int>,TH1F*>  _histTimeDiffs;

    ProditionsHandle<CRVOrdinal> _crvChannelMap_h;
  };

  CrvTimingStudies::CrvTimingStudies(const Parameters& conf) :
    art::EDAnalyzer{conf},
    _crvRecoPulsesModuleLabel(conf().crvRecoPulsesModuleLabel()),
    _PEthreshold(conf().PEthreshold())
  {
  }

  CrvTimingStudies::~CrvTimingStudies()
  {
  }

  void CrvTimingStudies::beginJob()
  {
  }

  void CrvTimingStudies::endJob()
  {
  }

  void CrvTimingStudies::analyze(const art::Event& event)
  {
    art::Handle<CrvRecoPulseCollection> crvRecoPulseCollection;
    if(!event.getByLabel(_crvRecoPulsesModuleLabel,"",crvRecoPulseCollection)) return;

    auto const& crvChannelMap = _crvChannelMap_h.get(event.id());

    std::map<uint16_t,std::vector<double> > fpgaTimes;
    for(auto crvRecoPulse=crvRecoPulseCollection->begin(); crvRecoPulse!=crvRecoPulseCollection->end(); ++crvRecoPulse)
    {
      if(crvRecoPulse->GetPEs()<_PEthreshold) continue;  //ignore pulses below a threshold
      if(crvRecoPulse->GetRecoPulseFlags().test(CrvRecoPulseFlagEnums::failedFit)) continue;  //ignore pulses with a failed fit

      auto barIndex         = crvRecoPulse->GetScintillatorBarIndex();
      int  SiPM             = crvRecoPulse->GetSiPMNumber();
      double recoPulseTime  = crvRecoPulse->GetPulseTime();

      uint16_t offlineChannel = barIndex.asUint()*CRVId::nChanPerBar + SiPM;
      CRVROC   onlineChannel  = crvChannelMap.online(offlineChannel);
      uint16_t ROC            = onlineChannel.ROC();
      uint16_t feb            = onlineChannel.FEB();
      uint16_t febChannel     = onlineChannel.FEBchannel();
ROC--;
feb--;
      uint16_t fpgaIndex      = (ROC*CRVId::nFEBPerROC*CRVId::nChanPerFEB+feb*CRVId::nChanPerFEB+febChannel)/(CRVId::nChanPerFEB/CRVId::nFPGAPerFEB);
//uint16_t fpgaIndex      = (ROC*CRVId::nFEBPerROC*CRVId::nChanPerFEB+feb*CRVId::nChanPerFEB+febChannel);

      fpgaTimes[fpgaIndex].push_back(recoPulseTime);
    }

    std::map<uint16_t,double> fpgaAverageTimes;
    for(auto fpga=fpgaTimes.begin(); fpga!=fpgaTimes.end(); ++fpga)
    {
      const auto &times=fpga->second;
      double averageTime=0;
      for(size_t i=0; i<times.size(); ++i) averageTime+=times.at(i);
      averageTime/=times.size();
      fpgaAverageTimes[fpga->first]=averageTime;
    }

    //compare time differences between different FPGAs
    for(auto fpga1=fpgaAverageTimes.begin(); fpga1!=fpgaAverageTimes.end(); fpga1++)
    for(auto fpga2=fpgaAverageTimes.begin(); fpga2!=fpgaAverageTimes.end(); fpga2++)
    {
      if(fpga1->first>=fpga2->first) continue; //don't compare with itself (=) and avoid comparing the same FPGAs twice (>)

      std::pair<int,int> histIndex(fpga1->first,fpga2->first);
      if(_histTimeDiffs.find(histIndex)==_histTimeDiffs.end())
      {
        art::ServiceHandle<art::TFileService> tfs;
        _histTimeDiffs[histIndex] = tfs->make<TH1F>(Form("fpgaTimeDiff_%i_%i",fpga1->first,fpga2->first),
                                                    Form("Time Diffs between FGPAs %i and %i;time difference [ns];Counts",fpga1->first,fpga2->first),
                                                    100,-50,50);
      }

      double timeDiff=fpga1->second-fpga2->second;
      _histTimeDiffs[histIndex]->Fill(timeDiff);
    }

  } // end analyze

} // end namespace mu2e

using mu2e::CrvTimingStudies;
DEFINE_ART_MODULE(CrvTimingStudies)
