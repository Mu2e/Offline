//
// Analyzer module to plot STM waveforms
//
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"

#include <utility>
#include <numeric>
#include <algorithm>
#include <cctype>
// root
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"

#include "Offline/RecoDataProducts/inc/STMWaveformDigi.hh"
#include "Offline/DataProducts/inc/STMChannel.hh"
#include "Offline/Mu2eUtilities/inc/STMUtils.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/STMConditions/inc/STMEnergyCalib.hh"

namespace mu2e {

  class PlotSTMWaveformDigis : public art::EDAnalyzer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<art::InputTag> stmWaveformDigisTag{ Name("stmWaveformDigisTag"), Comment("InputTag for STMWaveformDigiCollection")};
        fhicl::Atom<bool> subtractPedestal{ Name("subtractPedestal"), Comment("True/False whether to subtract the pedestal before plotting")};
        fhicl::Atom<std::string> xAxis{ Name("xAxis"), Comment("Choice of x-axis unit: \"sample_number\", \"adcs_time\", or \"event_time\"")} ;
        fhicl::Atom<int> verbosityLevel{ Name("verbosityLevel"), Comment("Verbosity level")};
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit PlotSTMWaveformDigis(const Parameters& conf);
     
    private:
    void beginJob() override;//For _hist
    void endJob() override; //For printing counter
    void analyze(const art::Event& e) override;
    
    TH1F* _hist; //Hist for WaveLength
    int ZeroLengthCount = 0;
    art::InputTag _stmWaveformDigisTag;

    art::ProductToken<STMWaveformDigiCollection> _stmWaveformDigisToken;
    bool _subtractPedestal;
    std::string _xAxis;
    int _verbosityLevel;
    ProditionsHandle<STMEnergyCalib> _stmEnergyCalib_h;
    STMChannel _channel;
  };

  PlotSTMWaveformDigis::PlotSTMWaveformDigis(const Parameters& config )  :
    art::EDAnalyzer{config},
    _stmWaveformDigisTag{config().stmWaveformDigisTag()},
    _stmWaveformDigisToken(consumes<STMWaveformDigiCollection>(config().stmWaveformDigisTag())),
    _subtractPedestal(config().subtractPedestal()),
    _xAxis(config().xAxis()),
    _verbosityLevel(config().verbosityLevel()),
    _channel(STMChannel::findByName("HPGe")) // FIXME: don't hardcode this probably don't want to do what we had before and try to infer it from the art::InputTag like this "STMUtils::getChannel(config().stmWaveformDigisTag()))"

  { }

  void PlotSTMWaveformDigis::beginJob(){

    art::ServiceHandle<art::TFileService> tfs;
    std::string X = std::string(_stmWaveformDigisTag.instance());
    std::transform(X.begin(),X.end(),X.begin(), toupper);
    std::string hWaveLength_title = "Waveform Lengths for " + X + " Pulses";
    _hist = tfs->make<TH1F>("hWaveLength", hWaveLength_title.c_str() ,1000,0,1000);
    
    }

  void PlotSTMWaveformDigis::endJob(){
    std::cout<< " Zero length Waveforms Count =  "<< ZeroLengthCount <<std::endl;
    }

  void PlotSTMWaveformDigis::analyze(const art::Event& event) {

    art::ServiceHandle<art::TFileService> tfs;
    auto waveformsHandle = event.getValidHandle(_stmWaveformDigisToken);
    std::stringstream histname, histtitle;
    int count = 0;
    STMEnergyCalib const& stmEnergyCalib = _stmEnergyCalib_h.get(event.id()); // get prodition
    const auto pedestal = stmEnergyCalib.pedestal(_channel);
    if (_verbosityLevel > 0) {
      std::cout << _channel.name() << " Pedestal = " << pedestal << std::endl;
    }


    const auto nsPerCt = stmEnergyCalib.nsPerCt(_channel);
    std::cout<<"size = "<<waveformsHandle->size()<<std::endl;
    for (const auto& waveform : *waveformsHandle) {

      histname.str("");
      histname << "evt" << event.event() << "_waveform" << count;
      histtitle.str("");
      histtitle << "Event " << event.event() << " Waveform " << count << " (" << _channel.name() << ")";

      if (waveform.adcs().size() == 0){
	  ZeroLengthCount++;
      } else {
	_hist->Fill(waveform.adcs().size()); //_hist was created outside so there should be no problem here

	std::stringstream histname2;//New histogram to include the offset waveform from trigTimeOffset
	histname2 << histname.str() << "_offset";
	
        Binning binning = STMUtils::getBinning(waveform, _xAxis, nsPerCt);
        TH1F* hWaveform = tfs->make<TH1F>(histname.str().c_str(), histtitle.str().c_str(), binning.nbins(),binning.low(),binning.high());
	TH1F* hWaveformOffset = tfs->make<TH1F>(histname2.str().c_str(),histtitle.str().c_str(),binning.nbins(),binning.low(),binning.high());//Exact copy except for name

	hWaveform->GetYaxis()->SetTitle("ADCs");
	hWaveform->GetXaxis()->SetTitle("Sample Number");

	hWaveformOffset->GetYaxis()->SetTitle("ADCs");
	hWaveformOffset->GetXaxis()->SetTitle("Sample Number");
	
	for(size_t i_adc = 0; i_adc < waveform.adcs().size();++i_adc){

          const auto adc = waveform.adcs().at(i_adc);
	  auto content = adc;
	  if (_subtractPedestal) {
	    content -= pedestal;
	  }  
	  hWaveform->SetBinContent(i_adc+1,content);
	  hWaveformOffset->SetBinContent(i_adc+1,content); //Exact copy
	}

	//For offset waveforms
	auto zs_offset = waveform.trigTimeOffset(); //Grabs stored offset
	int n_bins = hWaveformOffset->GetNbinsX(); //Grabs already contained nbins from waveform
	double xmin = hWaveformOffset->GetXaxis()->GetXmin();// gets xmin from waveform
	double xmax = hWaveformOffset->GetXaxis()->GetXmax();//gets xman from waveform

	hWaveformOffset->SetBins(n_bins, xmin + zs_offset, xmax + zs_offset);// shifts the xmin and xmax by offset, keeps numbers of bins
	  
	  
	}
      ++count;  
    }

  }
}
    

DEFINE_ART_MODULE(mu2e::PlotSTMWaveformDigis)
