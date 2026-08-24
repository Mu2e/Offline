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
#include <stdbool.h>
#include <map>

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
        fhicl::Atom<art::InputTag> stmWaveformDigisTag{ Name("stmWaveformDigisTag"),
            Comment("InputTag for STMWaveformDigiCollection")};
        fhicl::Atom<bool> subtractPedestal{ Name("subtractPedestal"),
            Comment("True/False whether to subtract the pedestal before plotting")};
        fhicl::Atom<std::string> xAxis{ Name("xAxis"),
            Comment("Choice of x-axis unit: \"sample_number\", \"waveform_time\", or \"event_time\"")} ;
        fhicl::Atom<int> verbosityLevel{ Name("verbosityLevel"),
            Comment("Verbosity level")};
        fhicl::Atom<bool> plotZSWithoutOffset{Name("plotZSWithoutOffset"),
            Comment("Whether to plot ZS without the trig timne offset"),
            false};
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit PlotSTMWaveformDigis(const Parameters& conf);

    private:
    void beginJob() override;//For _hist
    void endJob() override; //For printing counter
    void analyze(const art::Event& e) override;

    TH1F* _hist; //Hist for WaveLength
    int _zeroLengthCount = 0;
    art::InputTag _stmWaveformDigisTag;
    bool _plotZSWithoutOffset{false};

    //art::ProductToken<STMWaveformDigiCollection> _stmWaveformDigisToken;
    art::ProductToken<STMWaveformDigiCollectionMap> _stmWaveformDigisMapToken;
    bool _subtractPedestal;
    std::string _xAxis;
    int _verbosityLevel;
    ProditionsHandle<STMEnergyCalib> _stmEnergyCalib_h; //might have to change
    STMChannel _channel;
  };

  PlotSTMWaveformDigis::PlotSTMWaveformDigis(const Parameters& config )  :
    art::EDAnalyzer{config},
    _stmWaveformDigisTag{config().stmWaveformDigisTag()},
    //_stmWaveformDigisToken(consumes<STMWaveformDigiCollection>(config().stmWaveformDigisTag())),
    _plotZSWithoutOffset(config().plotZSWithoutOffset()),
    _stmWaveformDigisMapToken(consumes<STMWaveformDigiCollectionMap>(config().stmWaveformDigisTag())),
    _subtractPedestal(config().subtractPedestal()),
    _xAxis(config().xAxis()),
    _verbosityLevel(config().verbosityLevel()),
    _channel(STMUtils::getChannel(config().stmWaveformDigisTag()))

  { }

  void PlotSTMWaveformDigis::beginJob(){

    art::ServiceHandle<art::TFileService> tfs;
    std::string X = std::string(_stmWaveformDigisTag.instance()); //Gets instance name from fcl
    std::transform(X.begin(),X.end(),X.begin(), toupper); //Raises uppercase of DigiTag
    std::string hWaveLength_title = "Waveform Lengths for " + X + " Pulses"; //Builds title
    _hist = tfs->make<TH1F>("hWaveLength", hWaveLength_title.c_str() ,1000,0,1000); //makes the histogram hWaveLength

    }

  void PlotSTMWaveformDigis::endJob(){
    if (_verbosityLevel > 1){
      std::cout << " Zero length Waveforms Count =  " << _zeroLengthCount << std::endl;
    }
  }

  void PlotSTMWaveformDigis::analyze(const art::Event& event) {

    // Boolean for whether we get a match for zsHPGe or LaBr
    const std::string instance = std::string(_stmWaveformDigisTag.instance());
    const bool zsInstance = (instance == "zsHPGe" || instance == "zsLaBr");// Don't need offset waveforms for raw waveformdigis
    const bool rawInstance = (instance == "rawHPGe" || instance == "rawLaBr");

    if (!zsInstance){
        if (_verbosityLevel > 1){
            std::cout << "Instance : " << instance << " , not ZS ==> No Offset Waveform will be created" << std::endl;
        }
    }

    // STM Energy Calib
    art::ServiceHandle<art::TFileService> tfs;
    STMEnergyCalib const& stmEnergyCalib = _stmEnergyCalib_h.get(event.id()); // get prodition
    const auto pedestal = stmEnergyCalib.pedestal(_channel);
    if (_verbosityLevel > 0){
        std::cout << _channel.name() << " Pedestal = " << pedestal << std::endl;
    }

    const auto nsPerCt = stmEnergyCalib.nsPerCt(_channel);

    // Get handle for waveform handle with map
    // outer loop is needed now
    auto waveformsMapHandle = event.getValidHandle(_stmWaveformDigisMapToken);
    const auto& waveformMap = *waveformsMapHandle;
    // loop through the maps
    for (const auto& mu2e_evt : waveformMap){
        // Could get event information here
        const auto& eventHeader = mu2e_evt.first;
        // prodition may go here, not sure yet
        // Get waveform Digi collection information here
        const auto& waveforms = mu2e_evt.second; // treat this like a handle
        if (_verbosityLevel > 1) {
            std::cout << "Waveform size = " << waveforms.size() << std::endl;
        }

        // reset to zero when looping through new eventwindow tag
        std::stringstream histname, histtitle;
        std::stringstream histname2, histtitle2;
        std::stringstream histname3, histtitle3;
        int count = 0;

        // Second loop where rest of information goes
        for (const auto& waveform : waveforms){
            if (waveform.adcs().size() == 0){
                ++_zeroLengthCount;
            } else {
                // None empty waveforms go in here
                _hist->Fill(waveform.adcs().size());

                // Binning
                Binning binning = STMUtils::getBinning(waveform, _xAxis, nsPerCt); //nanosecondPetCount
                TH1F* hWaveform = nullptr;
                TH1F* hWaveformOffset = nullptr; // Standby
                TH1F* hWaveform3 = nullptr;

                if (rawInstance){
                    // histname in art
                    histname.str("");
                    //histname << "evt" << event.event() << "_waveform" << count;
                    histname << "ewt" << eventHeader.eventWindowTag() << "_waveform" << count << "_event" << event.event();
                    // histitle goes in the plot
                    histtitle.str("");
                    histtitle <<" EWT "<< eventHeader.eventWindowTag() << " Waveform " << count << " (" << _channel.name() << ")";

                    hWaveform = tfs->make<TH1F>(histname.str().c_str(), histtitle.str().c_str(),
                    binning.nbins(), binning.low(), binning.high());

                    // Get _xAxis for raw waveforms
                    hWaveform->GetYaxis()->SetTitle("ADCs");
                    if (_xAxis == "sample_number"){
                        hWaveform->GetXaxis()->SetTitle("Sample Number");
                    } else if (_xAxis == "waveform_time"){
                        hWaveform->GetXaxis()->SetTitle("Waveform Time [nsec]");
                    } else if (_xAxis == "event_time"){
                        hWaveform->GetXaxis()->SetTitle("Event Time [nsec]");
                    }

                } // end of raw instance

                // If the plotZS waveform is turned on
                if (zsInstance){
                    const auto zs_offset = waveform.trigTimeOffset(); // Grab stored offset

                    // descriptor
                    histname2.str("");
                    histname2 << "ewt"<< eventHeader.eventWindowTag() << "_waveform" << count << "_offset" << zs_offset << "_event" << event.event();
                    // title for plot
                    histtitle2.str("");
                    histtitle2 << " EWT " << eventHeader.eventWindowTag() << " " << instance << " Waveform " << count << " Offset " << zs_offset << " (" << _channel.name() << ")";
                    // make hist here
                    hWaveformOffset = tfs->make<TH1F>(histname2.str().c_str(), histtitle2.str().c_str(),
                    binning.nbins(), binning.low() + zs_offset, binning.high() + zs_offset); // Shifting bins using offset

                    // Get _xAxis for zsWaveforms
                    hWaveformOffset->GetYaxis()->SetTitle("ADCs");
                    if (_xAxis == "sample_number") {
                        hWaveformOffset->GetXaxis()->SetTitle("Sample Number (+Offset)");
                    } else if (_xAxis == "waveform_time") {
                        hWaveformOffset->GetXaxis()->SetTitle("Waveform Time [nsec] (+Offset)");
                    } else if (_xAxis == "event_time") {
                        hWaveformOffset->GetXaxis()->SetTitle("Event Time [nsec] (+Offset)");
                    }
                    // make unshfited set of plots
                    if (_plotZSWithoutOffset){
                        // here we plot for zs without the offset
                        // descriptor
                        histname3.str("");
                        histname3 << "ewt"<< eventHeader.eventWindowTag() << "_waveform" << count << "_offset" << zs_offset << "_event" << event.event() << "_unshifted";
                        // title for plot
                        histtitle3.str("");
                        histtitle3 << " EWT " << eventHeader.eventWindowTag() << " Unshifted " << instance << " Waveform " << " (" << _channel.name() << ")";
                        // make hist here
                        hWaveform3 = tfs->make<TH1F>(histname3.str().c_str(), histtitle3.str().c_str(),
                        binning.nbins(), binning.low(), binning.high()); // do not shift bins here

                        // Get _xAxis for zsWaveforms
                        hWaveform3->GetYaxis()->SetTitle("ADCs");
                        if (_xAxis == "sample_number") {
                            hWaveform3->GetXaxis()->SetTitle("Sample Number");
                        } else if (_xAxis == "waveform_time") {
                            hWaveform3->GetXaxis()->SetTitle("Waveform Time [nsec]");
                        } else if (_xAxis == "event_time") {
                            hWaveform3->GetXaxis()->SetTitle("Event Time [nsec]");
                        }
                    } // end of unshifted plots

                } // End of ZS instance

                // Fill waveforms here on plots that exist
                for (size_t i_adc = 0; i_adc < waveform.adcs().size(); ++i_adc){
                    const auto adc = waveform.adcs().at(i_adc);
                    auto content = adc; // y-axis
                    if (_subtractPedestal){
                        content -= pedestal;
                    }
                    if (rawInstance){
                        hWaveform->SetBinContent(i_adc + 1, content);
                    }
                    if (zsInstance){
                        hWaveformOffset->SetBinContent(i_adc + 1, content);
                        if (_plotZSWithoutOffset) {
                            hWaveform3 ->SetBinContent(i_adc + 1, content);
                        }
                    }
                    } // end of fill here
                } // else
                ++count;
            } // waveform loop
        }// map loop
    }// analyzer
} // nameSpace

DEFINE_ART_MODULE(mu2e::PlotSTMWaveformDigis)
