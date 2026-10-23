//
// A module to plot time differences between different FPGAs of FEBs
// to find the time offsets of all channels.
//
//
// Original Author: Ralf Ehrlich

#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "Offline/DataProducts/inc/CRVId.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/CRVConditions/inc/CRVOrdinal.hh"
#include "Offline/CRVConditions/inc/CRVCalib.hh"

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
#include <TCanvas.h>
#include <TVirtualPad.h>
#include <TText.h>
#include <TROOT.h>

#include <TError.h>

namespace
{
  ErrorHandlerFunc_t oldHandler = nullptr;

  void QuietPrintHandler(int level, Bool_t abort, const char* location, const char* msg)
  {
    // Suppress only TCanvas::Print informational messages.
    if (level == kInfo && strcmp(location, "TCanvas::Print") == 0) return;
    oldHandler(level, abort, location, msg);
  }
}

namespace mu2e
{
  class CrvTimeOffsets : public art::EDAnalyzer
  {
    public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct ModuleConfig
    {
      fhicl::Atom<std::string> name{Name("moduleName")};
      fhicl::Atom<int>         firstFEB{Name("firstFEB")};
      fhicl::Atom<int>         sides{Name("sides")};
      fhicl::Atom<bool>        muonTaggers{Name("muonTaggers")};
      fhicl::Atom<std::string> neighborName{Name("neighborName")};
    };
    struct TimeOffset
    {
      fhicl::Atom<int>         feb{Name("feb")};
      fhicl::Atom<double>      timeOffset{Name("timeOffset")};
    };
    struct Config
    {
      fhicl::Atom<std::string> crvRecoPulsesModuleLabel{Name("crvRecoPulseModuleLabel"), Comment("CrvRecoPulse Label")};
      fhicl::Atom<double>      PEthreshold{Name("PEthreshold"), Comment("PE threshold")};
      fhicl::Atom<bool>        removeTimeOffsets{Name("removeTimeOffsets"), Comment("remove time offsets added by reco")};
      fhicl::Atom<std::string> calibFileName{Name("calibFileName"), Comment("name of the DB file name for the time offsets")};
      fhicl::Atom<std::string> pdfFileName{Name("pdfFileName"), Comment("name of the pdf file name with the time difference plots")};
      fhicl::Atom<int>         diagLevel{fhicl::Name("diagLevel"), fhicl::Comment("diagnostic Level"), 0};

      //specifically for a geometry
      fhicl::Sequence<int>                         specialFEBsType1{Name("specialFEBsType1"), Comment("FEB numbers where the FPGAs are compared 0/2, 1/3, 1/2")};  //extracted: none
      fhicl::Sequence<int>                         specialFEBsType2{Name("specialFEBsType2"), Comment("FEB numbers where the FPGAs are compared 0/1, 1/2, 2/3")};  //extracted: 28,29
      fhicl::Sequence<fhicl::Table<ModuleConfig> > moduleConfigs{Name("moduleConfigs"), Comment("module-specific settings")};
      fhicl::Sequence<std::string>                 checkOppositeSides{Name("checkOppositeSides"), Comment("Check time diffs between opposite sides")};  //extracted: "EX_0", "T2_0"
      fhicl::Sequence<fhicl::Table<TimeOffset> >   specialTimeOffsets{Name("specialTimeOffsets"), Comment("special time offsets")};
    };

    using Parameters = art::EDAnalyzer::Table<Config>;
    explicit CrvTimeOffsets(const Parameters& conf);
    ~CrvTimeOffsets() override;
    void analyze(const art::Event& e) override;
    void beginJob() override;
    void endJob() override;
    void PreparePlots();

    private:
    std::string   _crvRecoPulsesModuleLabel;
    double        _PEthreshold;
    bool          _removeTimeOffsets;
    bool          _firstEvent;
    std::string   _calibFileName;
    std::string   _pdfFileName;
    int           _diagLevel;

    //specifically for a geometry
    std::vector<int>          _specialFEBsType1;
    std::vector<int>          _specialFEBsType2;
    std::vector<ModuleConfig> _moduleConfigs;
    std::vector<std::string>  _checkOppositeSides;
    std::vector<TimeOffset>   _specialTimeOffsets;

    std::map<std::pair<int,int>,TH1F*>  _histTimeDiffs; //between FPGAs
    std::map<int,CRVROC>                _channels;
    std::set<int>                       _febs;

    ProditionsHandle<CRVOrdinal> _crvChannelMap_h;
    ProditionsHandle<CRVCalib>   _calib_h;

    struct Plot
    {
      int   _pad;
      int   _feb1, _feb2;
      int   _fpga1, _fpga2;
      TH1F* _h;
      Plot(int pad, int feb1, int feb2, int fpga1, int fpga2, std::map<std::pair<int,int>,TH1F*> &histTimeDiffs) :
           _pad(pad), _feb1(feb1), _feb2(feb2), _fpga1(fpga1), _fpga2(fpga2)
      {
        std::pair<int,int> histIndex(feb1*CRVId::nFPGAPerFEB+fpga1,feb2*CRVId::nFPGAPerFEB+fpga2);
        art::ServiceHandle<art::TFileService> tfs;
        _h = tfs->make<TH1F>(Form("TimeDiffs__FEB_FPGA__%i_%i__%i_%i",feb1+1,fpga1,feb2+1,fpga2),
                             Form("FEB/FPGA %i/%i, FEB/FPGA %i/%i;time difference [ns];counts",feb1+1,fpga1,feb2+1,fpga2),
                             300,-150,150);
        histTimeDiffs[histIndex] = _h;
      }
      void FindMean(TVirtualPad *c, std::map<std::pair<int,int>,float> &measuredTimeDiffs)
      {
        c->cd(_pad);
        _h->Draw();
        float mean=_h->GetMean();  //TODO: get the mean from a Gaussian fit (where possible)
        if(_h->Integral(1,300)<25) //TODO: use a better cutoff
        {
          std::cerr<<"Couldn't find the mean for histogram "<<_h->GetName()<<"."<<std::endl;
        }
        else
        {
          std::pair<int,int> histIndex(_feb1*CRVId::nFPGAPerFEB+_fpga1,_feb2*CRVId::nFPGAPerFEB+_fpga2);
          measuredTimeDiffs[histIndex]=mean;
        }
      }
    };
    struct PlotPage
    {
      std::string       _name;
      std::vector<Plot> _plots;
      PlotPage(const std::string &name) : _name(name) {}
    };
    std::vector<PlotPage> _plotPages;
  };

  CrvTimeOffsets::CrvTimeOffsets(const Parameters& conf) :
    art::EDAnalyzer{conf},
    _crvRecoPulsesModuleLabel(conf().crvRecoPulsesModuleLabel()),
    _PEthreshold(conf().PEthreshold()),
    _removeTimeOffsets(conf().removeTimeOffsets()),
    _firstEvent(true),
    _calibFileName(conf().calibFileName()),
    _pdfFileName(conf().pdfFileName()),
    _diagLevel(conf().diagLevel()),
    _specialFEBsType1(conf().specialFEBsType1()),
    _specialFEBsType2(conf().specialFEBsType2()),
    _moduleConfigs(conf().moduleConfigs()),
    _checkOppositeSides(conf().checkOppositeSides()),
    _specialTimeOffsets(conf().specialTimeOffsets())
  {}

  void CrvTimeOffsets::PreparePlots() //create all required histograms of time differences between FPGAs
  {
    //(1) compare FPGAs of one FEB
    for(int feb : _febs)
    {
      std::string pageName(Form("FEB %i",feb+1));
      _plotPages.emplace_back(pageName);
      if(std::find(_specialFEBsType1.begin(),_specialFEBsType1.end(),feb)!=_specialFEBsType1.end()) //reverse offset
      {
        _plotPages.back()._plots.emplace_back(1,feb,feb,0,1,_histTimeDiffs);
        _plotPages.back()._plots.emplace_back(2,feb,feb,1,3,_histTimeDiffs);
        _plotPages.back()._plots.emplace_back(3,feb,feb,1,2,_histTimeDiffs);
      }
      else if(std::find(_specialFEBsType2.begin(),_specialFEBsType2.end(),feb)!=_specialFEBsType2.end()) //muon taggers
      {
        _plotPages.back()._plots.emplace_back(1,feb,feb,0,1,_histTimeDiffs);
        _plotPages.back()._plots.emplace_back(2,feb,feb,1,2,_histTimeDiffs);
        _plotPages.back()._plots.emplace_back(3,feb,feb,2,3,_histTimeDiffs);
      }
      else //regular modules
      {
        _plotPages.back()._plots.emplace_back(1,feb,feb,0,2,_histTimeDiffs);
        _plotPages.back()._plots.emplace_back(2,feb,feb,1,3,_histTimeDiffs);
        _plotPages.back()._plots.emplace_back(3,feb,feb,0,3,_histTimeDiffs);
      }
    }

    //(2) compare both FEBs at a module's readout side
    //module 8 is for the muon taggers
    for(auto moduleConfig : _moduleConfigs)
    {
      std::string pageName(Form("Module %s (FEBs on same side)",moduleConfig.name().c_str()));
      _plotPages.emplace_back(pageName);
      if(moduleConfig.muonTaggers())
      {
        _plotPages.back()._plots.emplace_back(1,moduleConfig.firstFEB()+0,moduleConfig.firstFEB()+1,3,0,_histTimeDiffs);
      }
      else
      {
        if(moduleConfig.sides()==2)
        {
          _plotPages.back()._plots.emplace_back(1,moduleConfig.firstFEB()+0,moduleConfig.firstFEB()+1,3,1,_histTimeDiffs); //one module side
          _plotPages.back()._plots.emplace_back(2,moduleConfig.firstFEB()+2,moduleConfig.firstFEB()+3,3,1,_histTimeDiffs); //other module side
        }
        else
        {
          _plotPages.back()._plots.emplace_back(1,moduleConfig.firstFEB()+0,moduleConfig.firstFEB()+1,3,1,_histTimeDiffs); //one module side
        }
      }
    }

    //(3) compare neighboring modules
    for(auto moduleConfig : _moduleConfigs)
    {
      std::string neighborName = moduleConfig.neighborName();
      if(neighborName=="") continue;  //last module in sector type

      auto neighborModuleConfig = std::find_if(_moduleConfigs.begin(),_moduleConfigs.end(),[&neighborName](const ModuleConfig &m){return m.name()==neighborName;});
      if(neighborModuleConfig==_moduleConfigs.end())
      {
        std::cerr<<"Neighbor module "<<neighborName<<" for module "<<moduleConfig.name()<<" does not exist."<<std::endl;
        continue;
      }

      std::string pageName(Form("Module %s / Module %s",moduleConfig.name().c_str(),neighborName.c_str()));
      _plotPages.emplace_back(pageName);

      if(!moduleConfig.muonTaggers())
      {
        if(moduleConfig.sides()==2)
        {
          _plotPages.back()._plots.emplace_back(1,moduleConfig.firstFEB()+0,neighborModuleConfig->firstFEB()+1,1,2,_histTimeDiffs); //one module side
          _plotPages.back()._plots.emplace_back(2,moduleConfig.firstFEB()+2,neighborModuleConfig->firstFEB()+3,1,2,_histTimeDiffs); //other module side
        }
        else
        {
          _plotPages.back()._plots.emplace_back(1,moduleConfig.firstFEB()+0,neighborModuleConfig->firstFEB()+1,1,2,_histTimeDiffs); //one module side
        }
      }
    }

    //(4) compare FEBs of opposite readout sides
    for(auto moduleConfig : _moduleConfigs)
    {
      std::string moduleName = moduleConfig.name();
      if(std::find(_checkOppositeSides.begin(),_checkOppositeSides.end(),moduleName)!=_checkOppositeSides.end())
      {
        std::string pageName(Form("Module %s (FEBs on opposite sides)",moduleName.c_str()));
        _plotPages.emplace_back(pageName);
        _plotPages.back()._plots.emplace_back(1,moduleConfig.firstFEB()+0,moduleConfig.firstFEB()+2,0,0,_histTimeDiffs);
      }
    }
  }

  CrvTimeOffsets::~CrvTimeOffsets()
  {
  }

  void CrvTimeOffsets::beginJob()
  {
  }

  void CrvTimeOffsets::endJob()
  {
    std::map<std::pair<int,int>,float> measuredTimeDiffs;
    std::map<int,float> timeOffsets;

    for(const auto& [feb,specialTimeOffset] : _specialTimeOffsets)
    {
      timeOffsets[feb()*CRVId::nFPGAPerFEB]=specialTimeOffset();  //reference time offsets to start chains of FPGAs with time differences
    }

    oldHandler = SetErrorHandler(QuietPrintHandler);  //suppress the messages caused by TCanvas::Print()
    TCanvas c0;
    c0.Print(Form("%s[", _pdfFileName.c_str()), "pdf");

    //Find all time differences
    for(size_t i=0; i<_plotPages.size(); ++i)
    {
      TCanvas *c = new TCanvas(_plotPages.at(i)._name.c_str(),_plotPages.at(i)._name.c_str(),800,800);
      TText *t = new TText(.1,.8,_plotPages.at(i)._name.c_str());
      gROOT->Add(c);
      gROOT->Add(t);
      c->Divide(1,2);

      c->cd(1);
      t->SetTextSize(0.08);
      t->Draw();

      c->cd(2);
      TVirtualPad *pad=gPad;
      pad->Divide(_plotPages.at(i)._plots.size(),1);

      for(size_t j=0; j<_plotPages.at(i)._plots.size(); ++j)
      {
        _plotPages.at(i)._plots.at(j).FindMean(pad, measuredTimeDiffs);
      }

      c->Print(_pdfFileName.c_str(),"pdf");
    }

    c0.Print(Form("%s]", _pdfFileName.c_str()), "pdf");
    SetErrorHandler(oldHandler);

    //Find a chain of FPGAs starting from a reference FPGA to find time offsets all FPGAs relative to the reference FPGA
    size_t prevSize=measuredTimeDiffs.size();
    for(auto measuredTimeDiff=measuredTimeDiffs.cbegin(); ; )
    {
      if(measuredTimeDiff==measuredTimeDiffs.cend())
      {
        if(measuredTimeDiffs.size()==0) break;  //all measured time differences used.
        if(prevSize==measuredTimeDiffs.size()) break;  //loops don't seem to make a difference anymore.

        prevSize=measuredTimeDiffs.size();
        measuredTimeDiff=measuredTimeDiffs.cbegin();  //there are still some measured time differences left. start the loop again.
      }

      //check if this particular measured time difference can be connected to a point where the time offset is known already.
      bool erased=false;
      for(const auto& timeOffset: timeOffsets)
      {
        //Longer cables result in an earlier reco time, because the start t0 arrives later at the FEBs, so that the time difference between t0 and tHit is shorter (=tReco).
        //Therefore, the dt plots result in a negative mean value.
        //The time offset needs to be positive to counter act it. That's why a negative sign is used below.
        if(measuredTimeDiff->first.first==timeOffset.first)
        {
          timeOffsets[measuredTimeDiff->first.second]=timeOffset.second+measuredTimeDiff->second;
          measuredTimeDiff=measuredTimeDiffs.erase(measuredTimeDiff);
          erased=true;
          break;
        }
        if(measuredTimeDiff->first.second==timeOffset.first)
        {
          timeOffsets[measuredTimeDiff->first.first]=timeOffset.second-measuredTimeDiff->second;
          measuredTimeDiff=measuredTimeDiffs.erase(measuredTimeDiff);
          erased=true;
          break;
        }
      }
      if(erased) continue;

      ++measuredTimeDiff;
    } //done with connecting all FPGAs

    if(measuredTimeDiffs.size()>0) std::cerr<<"There are still some unused measured time diffs!"<<std::endl;
    for(const auto& measuredTimeDiff: measuredTimeDiffs)
    {
      std::cerr<<"Unused time diff for "<<measuredTimeDiff.first.first<<"/"<<measuredTimeDiff.first.second<<": "<<measuredTimeDiff.second<<std::endl;
    }

    if(_diagLevel>0)
    {
      for(const auto& timeOffset: timeOffsets)
      {
        int fpga=timeOffset.first%CRVId::nFPGAPerFEB;
        int globalfeb=timeOffset.first/CRVId::nFPGAPerFEB;
        int feb=globalfeb%CRVId::nFEBPerROC+1;
        int roc=globalfeb/CRVId::nFEBPerROC+1;
        std::cout<<"ROC "<<roc<<"    FEB "<<feb<<"   fpga "<<fpga<<"     timeOffset "<<timeOffset.second<<std::endl;
      }
    }

    //write out time calibration file
    std::ofstream calibFile(_calibFileName);
    if(!calibFile.is_open()) throw cet::exception("CRVCALIB") << "Couldn't create output txt file " << _calibFileName << ".";
    calibFile<<"TABLE CRVTime"<<std::endl;
    calibFile<<"#channel,timeOffset"<<std::endl;

    //find time offsets for each offline channel
    int nChanPerFPGA = CRVId::nChanPerFEB/CRVId::nFPGAPerFEB;
    GeomHandle<CosmicRayShield> CRS;
    const std::vector<std::shared_ptr<CRSScintillatorBar> > &counters = CRS->getAllCRSScintillatorBars();
    for(size_t offlineChannel=0; offlineChannel<counters.size()*CRVId::nChanPerBar; ++offlineChannel)
    {
      if(_channels.find(offlineChannel)==_channels.end())
      {
        calibFile<<offlineChannel<<",0"<<std::endl;
      }
      else
      {
        const CRVROC &onlineChannel = _channels.find(offlineChannel)->second;
        int globalfeb    = (onlineChannel.ROC()-1)*CRVId::nFEBPerROC + (onlineChannel.FEB()-1);
        int globalfpga   = globalfeb*CRVId::nFPGAPerFEB + onlineChannel.FEBchannel()/nChanPerFPGA;
        if(timeOffsets.find(globalfpga)==timeOffsets.end())
        {
          std::cerr<<"Couldn't find a time offset for roc "<<onlineChannel.ROC()<<"  feb "<<onlineChannel.FEB()<<"  febChannel "<<onlineChannel.FEBchannel()<<std::endl;
          calibFile<<offlineChannel<<",0"<<std::endl;
        }
        else calibFile<<offlineChannel<<","<<timeOffsets.at(globalfpga)<<std::endl;
      }
    }

    calibFile.close();
  } //endJob

  void CrvTimeOffsets::analyze(const art::Event& event)
  {
    art::Handle<CrvRecoPulseCollection> crvRecoPulseCollection;
    if(!event.getByLabel(_crvRecoPulsesModuleLabel,"",crvRecoPulseCollection)) return;

    auto const& crvChannelMap = _crvChannelMap_h.get(event.id());
    auto const& calib = _calib_h.get(event.id());

    if(_firstEvent)
    {
      _firstEvent=false;

      GeomHandle<CosmicRayShield> CRS;
      const std::vector<std::shared_ptr<CRSScintillatorBar> > &counters = CRS->getAllCRSScintillatorBars();
      for(size_t offlineChannel=0; offlineChannel<counters.size()*CRVId::nChanPerBar; ++offlineChannel)
      {
        if(!crvChannelMap.onlineExists(offlineChannel)) continue;
        CRVROC onlineChannel = crvChannelMap.online(offlineChannel);
        _channels.emplace(offlineChannel,onlineChannel);
        int feb = (onlineChannel.ROC()-1)*CRVId::nFEBPerROC+(onlineChannel.FEB()-1);
        _febs.insert(feb);
      }

      PreparePlots();
    }


    int previousOfflineChannel=-1;
    std::map<uint16_t,std::vector<double> > fpgaTimes;
    for(auto crvRecoPulse=crvRecoPulseCollection->begin(); crvRecoPulse!=crvRecoPulseCollection->end(); ++crvRecoPulse)
    {
      if(crvRecoPulse->GetPEs()<_PEthreshold) continue;  //ignore pulses below a threshold
      if(crvRecoPulse->GetRecoPulseFlags().test(CrvRecoPulseFlagEnums::failedFit)) continue;  //ignore pulses with a failed fit

      auto barIndex         = crvRecoPulse->GetScintillatorBarIndex();
      int  SiPM             = crvRecoPulse->GetSiPMNumber();
      double recoPulseTime  = crvRecoPulse->GetPulseTime();

      uint16_t offlineChannel = barIndex.asUint()*CRVId::nChanPerBar + SiPM;
      if(offlineChannel==previousOfflineChannel) continue; //only look at first pulse of a channel (pulses arrive time ordered)
      previousOfflineChannel=offlineChannel;
      CRVROC   onlineChannel  = crvChannelMap.online(offlineChannel);
      uint16_t ROC            = onlineChannel.ROC();
      uint16_t feb            = onlineChannel.FEB();
      uint16_t febChannel     = onlineChannel.FEBchannel();

      int      nChanPerFPGA   = CRVId::nChanPerFEB/CRVId::nFPGAPerFEB;
      uint16_t fpgaIndex      = ((ROC-1)*CRVId::nFEBPerROC*CRVId::nChanPerFEB+(feb-1)*CRVId::nChanPerFEB+febChannel)/nChanPerFPGA;

      if(_removeTimeOffsets) //remove time offsets introduced during reconstruction to see the "pure" time differences, e.g. to generate new time calibration tables
      {
        double timeOffset = calib.timeOffset(offlineChannel);
        recoPulseTime-=timeOffset;
      }

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
    for(const auto& [histIndex, hist] : _histTimeDiffs)
    {
      int fpga1=histIndex.first;
      int fpga2=histIndex.second;
      auto avgTime1 = fpgaAverageTimes.find(fpga1);
      auto avgTime2 = fpgaAverageTimes.find(fpga2);
      if(avgTime1==fpgaAverageTimes.end() || avgTime2==fpgaAverageTimes.end()) continue;  //plot cannot be created
      double timeDiff=avgTime1->second-avgTime2->second;
      hist->Fill(timeDiff);
    }

  } // end analyze

} // end namespace mu2e

using mu2e::CrvTimeOffsets;
DEFINE_ART_MODULE(CrvTimeOffsets)
