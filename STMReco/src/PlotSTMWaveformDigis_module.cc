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
        fhicl::Atom<art::InputTag> stmWaveformDigisMapTag{ Name("stmWaveformDigisMapTag"),
            Comment("InputTag for STMWaveformDigiCollectionMap")};
        fhicl::Atom<std::string> waveformType{ Name("waveformType"),
            Comment("Type of waveform to plot: \"raw\" or \"zs\"")};
        fhicl::Atom<bool> subtractPedestal{ Name("subtractPedestal"),
            Comment("True/False whether to subtract the pedestal before plotting")};
        fhicl::Atom<std::string> xAxis{ Name("xAxis"),
            Comment("Choice of x-axis unit: \"sample_number\", \"waveform_time\", or \"event_time\"")} ;
        fhicl::Atom<int> verbosityLevel{ Name("verbosityLevel"),
            Comment("Verbosity level")};
        fhicl::Atom<bool> plotZSWithoutOffset{Name("plotZSWithoutOffset"),
            Comment("Whether to plot ZS without the trig time offset"),
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
    art::InputTag _stmWaveformDigisMapTag;
    bool _plotZSWithoutOffset{false};

    //art::ProductToken<STMWaveformDigiCollection> _stmWaveformDigisToken;
    art::ProductToken<STMWaveformDigiCollectionMap> _stmWaveformDigisMapToken;
    std::string _waveformType;
    bool _subtractPedestal;
    std::string _xAxis;
    int _verbosityLevel;
    ProditionsHandle<STMEnergyCalib> _stmEnergyCalib_h; //might have to change
    STMChannel _channel;
  };

  PlotSTMWaveformDigis::PlotSTMWaveformDigis(const Parameters& config )  :
    art::EDAnalyzer{config},
    _stmWaveformDigisMapTag{config().stmWaveformDigisMapTag()},
    _plotZSWithoutOffset(config().plotZSWithoutOffset()),
    _stmWaveformDigisMapToken(consumes<STMWaveformDigiCollectionMap>(config().stmWaveformDigisMapTag())),
    _waveformType(config().waveformType()),
    _subtractPedestal(config().subtractPedestal()),
    _xAxis(config().xAxis()),
    _verbosityLevel(config().verbosityLevel()),
    _channel(STMUtils::getChannel(config().stmWaveformDigisMapTag()))

  {
    if (_waveformType != "raw" && _waveformType != "zs") {
        throw cet::exception("Configuration", "Invalid waveformType: " + _waveformType + ". Must be \"raw\" or \"zs\".");
    }
    if (_waveformType == "zs" && _plotZSWithoutOffset && _xAxis == "event_time" ) {
        throw cet::exception("Configuration", "plotZSWithoutOffset cannot be true when xAxis is \"event_time\". Use \"sample_number\" or \"waveform_time\" instead.");
    }
  }

  void PlotSTMWaveformDigis::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    std::string X = std::string(_stmWaveformDigisMapTag.instance()); //Gets instance name from fcl
    std::transform(X.begin(),X.end(),X.begin(), toupper); //Raises uppercase of DigiTag
    std::string hWaveLength_title = "Waveform Lengths for " + X + " Pulses"; //Builds title
    _hist = tfs->make<TH1F>("hWaveLength", hWaveLength_title.c_str() ,1000,0,1000); //makes the histogram hWaveLength
    }

  void PlotSTMWaveformDigis::endJob() {
    if (_verbosityLevel > 1){
      std::cout << " Zero length Waveforms Count =  " << _zeroLengthCount << std::endl;
    }
  }

  void PlotSTMWaveformDigis::analyze(const art::Event& event) {

    // Boolean for whether we get a match for zsHPGe or LaBr
    const std::string instance = std::string(_stmWaveformDigisMapTag.instance());
    const bool isZS = (_waveformType == "zs");
    const bool isRaw = (_waveformType == "raw");

    if (!isZS){
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

        std::stringstream histnameRaw, histtitleRaw;
        std::stringstream histnameZSOffset, histtitleZSOffset;
        std::stringstream histnameZSUnshifted, histtitleZSUnshifted;
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
                TH1F* hWaveformRaw = nullptr;
                TH1F* hWaveformZSOffset = nullptr; // Standby
                TH1F* hWaveformZSUnshifted = nullptr;

                if (isRaw){
                    // histname in art
                    histnameRaw.str("");
                    //histname << "evt" << event.event() << "_waveform" << count;
                    histnameRaw << "ewt" << eventHeader.eventWindowTag() << "_waveform" << count << "_event" << event.event();
                    // histitle goes in the plot
                    histtitleRaw.str("");
                    histtitleRaw <<" EWT "<< eventHeader.eventWindowTag() << " Waveform " << count << " (" << _channel.name() << ")";

                    hWaveformRaw = tfs->make<TH1F>(histnameRaw.str().c_str(), histtitleRaw.str().c_str(),
                    binning.nbins(), binning.low(), binning.high());

                    // Get _xAxis for raw waveforms
                    hWaveformRaw->GetYaxis()->SetTitle("ADCs");
                    if (_xAxis == "sample_number"){
                        hWaveformRaw->GetXaxis()->SetTitle("Sample Number");
                    } else if (_xAxis == "waveform_time"){
                        hWaveformRaw->GetXaxis()->SetTitle("Waveform Time [nsec]");
                    } else if (_xAxis == "event_time"){
                        hWaveformRaw->GetXaxis()->SetTitle("Event Time [nsec]");
                    }

                } // end of raw instance

                // If the plotZS waveform is turned on
                if (isZS){
                    const auto zs_offset = waveform.trigTimeOffset(); // Grab stored offset

                    // descriptor
                    histnameZSOffset.str("");
                    histnameZSOffset << "ewt"<< eventHeader.eventWindowTag() << "_waveform" << count << "_offset" << zs_offset << "_event" << event.event();
                    // title for plot
                    histtitleZSOffset.str("");
                    histtitleZSOffset << " EWT " << eventHeader.eventWindowTag() << " " << instance << " Waveform " << count << " Offset " << zs_offset << " (" << _channel.name() << ")";

                    //
                    double plotOffset = 0.0;
                    if (_xAxis == "sample_number") {
                        plotOffset = zs_offset;
                    } else if (_xAxis == "waveform_time") {
                        plotOffset = zs_offset * nsPerCt;
                    } else if (_xAxis == "event_time") {
                        plotOffset = 0;
                    }

                    // make hist here
                    hWaveformZSOffset = tfs->make<TH1F>(histnameZSOffset.str().c_str(), histtitleZSOffset.str().c_str(),
                    binning.nbins(), binning.low() + plotOffset, binning.high() + plotOffset); // Shifting bins using plotOffset

                    // Get _xAxis for zsWaveforms
                    hWaveformZSOffset->GetYaxis()->SetTitle("ADCs");
                    if (_xAxis == "sample_number") {
                        hWaveformZSOffset->GetXaxis()->SetTitle("Sample Number (+Offset)");
                    } else if (_xAxis == "waveform_time") {
                        hWaveformZSOffset->GetXaxis()->SetTitle("Waveform Time [nsec] (+Offset)");
                    } else if (_xAxis == "event_time") {
                        hWaveformZSOffset->GetXaxis()->SetTitle("Event Time [nsec]");
                    }
                    // make unshfited set of plots
                    if (_plotZSWithoutOffset && _xAxis != "event_time") {
                        // here we plot for zs without the offset
                        // descriptor
                        histnameZSUnshifted.str("");
                        histnameZSUnshifted << "ewt"<< eventHeader.eventWindowTag() << "_waveform" << count << "_offset" << zs_offset << "_event" << event.event() << "_unshifted";
                        // title for plot
                        histtitleZSUnshifted.str("");
                        histtitleZSUnshifted << " EWT " << eventHeader.eventWindowTag() << " Unshifted " << instance << " Waveform " << " (" << _channel.name() << ")";

                        // make hist here
                        hWaveformZSUnshifted = tfs->make<TH1F>(histnameZSUnshifted.str().c_str(), histtitleZSUnshifted.str().c_str(),
                        binning.nbins(), binning.low(), binning.high()); // do not shift bins here

                        // Get _xAxis for zsWaveforms
                        hWaveformZSUnshifted->GetYaxis()->SetTitle("ADCs");
                        if (_xAxis == "sample_number") {
                            hWaveformZSUnshifted->GetXaxis()->SetTitle("Sample Number");
                        } else if (_xAxis == "waveform_time") {
                            hWaveformZSUnshifted->GetXaxis()->SetTitle("Waveform Time [nsec]");
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
                  if (isRaw){
                    hWaveformRaw->SetBinContent(i_adc + 1, content);
                  }
                  if (isZS){
                    hWaveformZSOffset->SetBinContent(i_adc + 1, content);
                    if (_plotZSWithoutOffset && _xAxis != "event_time") {
                      hWaveformZSUnshifted->SetBinContent(i_adc + 1, content);
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
