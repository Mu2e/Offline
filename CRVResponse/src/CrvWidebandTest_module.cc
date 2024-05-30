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

namespace
{
double LandauGaussFunction(double *x, double *par)
{
    //From $ROOTSYS/tutorials/fit/langaus.C
    //Fit parameters:
    //par[0]=Width (scale) parameter of Landau density
    //par[1]=Most Probable (MP, location) parameter of Landau density
    //par[2]=Total area (integral -inf to inf, normalization constant)
    //par[3]=Width (sigma) of convoluted Gaussian function
    //
    //In the Landau distribution (represented by the CERNLIB approximation),
    //the maximum is located at x=-0.22278298 with the location parameter=0.
    //This shift is corrected within this function, so that the actual
    //maximum is identical to the MP parameter.

    // Numeric constants
    constexpr Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
    constexpr Double_t mpshift  = -0.22278298;       // Landau maximum location

    // Control constants
    constexpr Double_t np = 100.0;      // number of convolution steps
    constexpr Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

    // Variables
    Double_t xx = 0.0;
    Double_t mpc = 0.0;
    Double_t fland = 0.0;
    Double_t sum = 0.0;
    Double_t xlow = 0.0, xupp = 0.0;
    Double_t step = 0.0;
    Int_t    i = 0.0;

    // MP shift correction
    mpc = par[1] - mpshift * par[0];

    // Range of convolution integral
    xlow = x[0] - sc * par[3];
    xupp = x[0] + sc * par[3];
    step = (xupp-xlow) / np;

    // Convolution integral of Landau and Gaussian by sum
    for(i=1.0; i<=np/2; i++)
    {
      xx = xlow + (i-.5) * step;
      fland = TMath::Landau(xx,mpc,par[0]) / par[0];
      sum += fland * TMath::Gaus(x[0],xx,par[3]);

      xx = xupp - (i-.5) * step;
      fland = TMath::Landau(xx,mpc,par[0]) / par[0];
      sum += fland * TMath::Gaus(x[0],xx,par[3]);
    }

    return (par[2] * step * sum * invsq2pi / par[3]);
}
void LandauGauss(TH1F &h, float &mpv, float &fwhm, float &signals, float &chi2)
{
    std::multimap<float,float> bins;  //binContent,binCenter
    for(int i=8; i<=h.GetNbinsX(); i++) bins.emplace(h.GetBinContent(i),h.GetBinCenter(i));  //ordered from smallest to largest bin entries
    if(bins.size()<4) return;
    if(bins.rbegin()->first<20) return;  //low statistics
    int nBins=0;
    float binSum=0;
    for(auto bin=bins.rbegin(); bin!=bins.rend(); ++bin)
    {
      nBins++;
      binSum+=bin->second;
      if(nBins==4) break;
    }
    float maxX=binSum/4;
    float fitRangeStart=0.7*maxX;  //0.6 @ 24
    float fitRangeEnd  =2.0*maxX;
    if(fitRangeStart<15.0) fitRangeStart=15.0;

    //Parameters
    Double_t startValues[4], parLimitsLow[4], parLimitsHigh[4];
    //Most probable value
    startValues[1]=maxX;
    parLimitsLow[1]=fitRangeStart;
    parLimitsHigh[1]=fitRangeEnd;
    //Area
    startValues[2]=h.Integral(h.FindBin(fitRangeStart),h.FindBin(fitRangeEnd));
    parLimitsLow[2]=0.01*startValues[2];
    parLimitsHigh[2]=100*startValues[2];
    //Other parameters
    startValues[0]=5.0;   startValues[3]=10.0;
    parLimitsLow[0]=2.0;  parLimitsLow[3]=2.0;
    parLimitsHigh[0]=15.0; parLimitsHigh[3]=20.0; //7 and 15 @ 21  //6 and 13 @ 23

    TF1 fit("LandauGauss",LandauGaussFunction,fitRangeStart,fitRangeEnd,4);
    fit.SetParameters(startValues);
    fit.SetLineColor(kRed);
    fit.SetParNames("Width","MP","Area","GSigma");
    for(int i=0; i<4; i++) fit.SetParLimits(i, parLimitsLow[i], parLimitsHigh[i]);
    TFitResultPtr fr = h.Fit(&fit,"LQRS");
    fit.Draw("same");

    mpv = fit.GetMaximumX();
    chi2 = (fr->Ndf()>0?fr->Chi2()/fr->Ndf():NAN);
    if(mpv==fitRangeStart) {mpv=0; return;}
    float halfMaximum = fit.Eval(mpv)/2.0;
    float leftX = fit.GetX(halfMaximum,0.0,mpv);
    float rightX = fit.GetX(halfMaximum,mpv,10.0*mpv);
    fwhm = rightX-leftX;
    signals = fit.Integral(0,150);
}

} //end anonymous namespace for LandauGauss function


namespace mu2e
{
  class CrvWidebandTest : public art::EDAnalyzer
  {
    public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct Config
    {
      fhicl::Atom<std::string> crvStepsModuleLabel{ Name("crvStepModuleLabel"), Comment("CrvStep label")};
      fhicl::Atom<std::string> crvRecoPulsesModuleLabel{ Name("crvRecoPulseModuleLabel"), Comment("CrvRecoPulse Label")};
      fhicl::Atom<std::string> crvCoincidenceClusterModuleLabel{ Name("crvCoincidenceClusterModuleLabel"), Comment("CrvCoincidenceCluster Label")};
      fhicl::Atom<std::string> crvCoincidenceClusterMCModuleLabel{ Name("crvCoincidenceClusterMCModuleLabel"), Comment("CrvCoincidenceClusterMC Label")};
      fhicl::Atom<float> minTrackFitPEs{ Name("minTrackFitPEs"), Comment("minimum number of PEs to be considered for track fit")};
      fhicl::Atom<art::InputTag> protonBunchTimeTag{ Name("protonBunchTimeTag"), Comment("ProtonBunchTime producer"),"EWMProducer" };
    };

    using Parameters = art::EDAnalyzer::Table<Config>;
    explicit CrvWidebandTest(const Parameters& conf);
    ~CrvWidebandTest() override;
    void analyze(const art::Event& e) override;
    void beginJob() override;
    void endJob() override;

    private:
    std::string   _crvStepsModuleLabel;
    std::string   _crvRecoPulsesModuleLabel;
    std::string   _crvCoincidenceClusterModuleLabel;
    std::string   _crvCoincidenceClusterMCModuleLabel;
    const float   _minTrackFitPEs{5};
    art::InputTag _protonBunchTimeTag;

    int     _nFEBs{0};
    int     _nSectorTypes{0};

    int     _runNumber{0};
    int     _subrunNumber{0};
    int     _spillNumber{1};  //to have the same variables as in the Wideband files
    int     _spillIndex{1};   //to have the same variables as in the Wideband files
    int     _eventNumber{0};

    float  *_recoPEs{NULL};
    float  *_recoTime{NULL};
    int    *_fitStatus{NULL};
    float  *_depositedEnergy{NULL};
    int    *_coincidencePDGid{NULL};
    float  *_coincidenceTime{NULL};
    float  *_coincidencePosX{NULL};
    float  *_coincidencePosY{NULL};
    float  *_coincidencePosZ{NULL};

    float  *_trackSlope{NULL};      //using slope=dx/dy to avoid inf for vertical tracks
    float  *_trackIntercept{NULL};  //x value, where y=0
    int    *_trackPoints{NULL};
    float  *_trackPEs{NULL};
    float  *_trackChi2{NULL};

    int     _eventsRecorded{0};
    float  *_summaryPEs{NULL};
    float  *_summaryFWHMs{NULL};
    float  *_summarySignals{NULL};
    float  *_summaryChi2s{NULL};

    std::vector<TH1F*> _histPEs;

    TTree  *_tree;
    TTree  *_treeSummary;

    ProditionsHandle<CRVOrdinal> _crvChannelMap_h;
  };

  CrvWidebandTest::CrvWidebandTest(const Parameters& conf) :
    art::EDAnalyzer{conf},
    _crvStepsModuleLabel(conf().crvStepsModuleLabel()),
    _crvRecoPulsesModuleLabel(conf().crvRecoPulsesModuleLabel()),
    _crvCoincidenceClusterModuleLabel(conf().crvCoincidenceClusterModuleLabel()),
    _crvCoincidenceClusterMCModuleLabel(conf().crvCoincidenceClusterMCModuleLabel()),
    _protonBunchTimeTag(conf().protonBunchTimeTag())
  {
    art::ServiceHandle<art::TFileService> tfs;
    _tree = tfs->make<TTree>("run", "run");
    _treeSummary = tfs->make<TTree>("runSummary", "runSummary");
  }

  CrvWidebandTest::~CrvWidebandTest()
  {
    delete[] _recoPEs;
    delete[] _depositedEnergy;
    delete[] _coincidencePDGid;
  }

  void CrvWidebandTest::beginJob()
  {
  }

  void CrvWidebandTest::endJob()
  {
    _summaryPEs     = new float[_nFEBs*CRVId::nChanPerFEB];
    _summaryFWHMs   = new float[_nFEBs*CRVId::nChanPerFEB];
    _summarySignals = new float[_nFEBs*CRVId::nChanPerFEB];
    _summaryChi2s   = new float[_nFEBs*CRVId::nChanPerFEB];

    _treeSummary->Branch("runNumber",&_runNumber,"runNumber/I");
    _treeSummary->Branch("subrunNumber",&_subrunNumber,"subrunNumber/I");
    _treeSummary->Branch("eventsRecorded",&_eventsRecorded,"eventsRecorded/I");
    _treeSummary->Branch("PEs",_summaryPEs,Form("PEs[%i][%i]/F",_nFEBs,(unsigned int)CRVId::nChanPerFEB));
    _treeSummary->Branch("FWHWs",_summaryFWHMs,Form("FWHMs[%i][%i]/F",_nFEBs,(unsigned int)CRVId::nChanPerFEB));
    _treeSummary->Branch("signals",_summarySignals,Form("signals[%i][%i]/F",_nFEBs,(unsigned int)CRVId::nChanPerFEB));
    _treeSummary->Branch("chi2s",_summaryChi2s,Form("chi2s[%i][%i]/F",_nFEBs,(unsigned int)CRVId::nChanPerFEB));

    for(int iFEB=0; iFEB<_nFEBs; ++iFEB)
    for(int iChannel=0; iChannel<(int)CRVId::nChanPerFEB; ++iChannel)
    {
      int index=iFEB*CRVId::nChanPerFEB+iChannel;
      _summaryPEs[index]=0;
      _summaryFWHMs[index]=0;
      _summarySignals[index]=0;
      _summaryChi2s[index]=0;
      LandauGauss(*_histPEs[index], _summaryPEs[index], _summaryFWHMs[index], _summarySignals[index], _summaryChi2s[index]);
    }

    _treeSummary->Fill();
  }

  void CrvWidebandTest::analyze(const art::Event& event)
  {
    art::Handle<CrvStepCollection> crvStepsCollection;
    if(!event.getByLabel(_crvStepsModuleLabel,"",crvStepsCollection)) return;

    art::Handle<CrvRecoPulseCollection> crvRecoPulseCollection;
    if(!event.getByLabel(_crvRecoPulsesModuleLabel,"",crvRecoPulseCollection)) return;

    art::Handle<CrvCoincidenceClusterCollection> crvCoincidenceClusterCollection;
    if(!event.getByLabel(_crvCoincidenceClusterModuleLabel,"",crvCoincidenceClusterCollection)) return;

    art::Handle<CrvCoincidenceClusterMCCollection> crvCoincidenceClusterMCCollection;
    if(!event.getByLabel(_crvCoincidenceClusterMCModuleLabel,"",crvCoincidenceClusterMCCollection)) return;

    art::Handle<ProtonBunchTime> protonBunchTime;
    event.getByLabel(_protonBunchTimeTag, protonBunchTime);
//    double TDC0time = -protonBunchTime->pbtime_;

    _runNumber = event.run();
    _subrunNumber = event.subRun();
    _eventNumber = event.event();
    ++_eventsRecorded;

    auto const& crvChannelMap = _crvChannelMap_h.get(event.id());

    GeomHandle<CosmicRayShield> CRS;
    const std::vector<std::shared_ptr<CRSScintillatorBar> > &counters = CRS->getAllCRSScintillatorBars();
    std::vector<std::shared_ptr<CRSScintillatorBar> >::const_iterator iter;

    //setup ROOT tree at the first event
    if(_recoPEs==NULL)
    {
      //find out how many FEBs are present and what the top/bottom locations are
      //find out how many sector types are present for which coincidences can be found
      _nFEBs = 0;
      _nSectorTypes = 0;
      for(iter=counters.begin(); iter!=counters.end(); iter++)
      {
        //get counter properties
        const CRSScintillatorBar &crvCounter = **iter;
        const CRSScintillatorBarIndex &barIndex = crvCounter.index();
        int sectorNumber = crvCounter.id().getShieldNumber();
        int sectorType = CRS->getCRSScintillatorShields().at(sectorNumber).getSectorType();
        if(sectorType>=_nSectorTypes) _nSectorTypes=sectorType+1;

        for(int SiPM=0; SiPM<4; SiPM++)
        {
          if(!crvCounter.getBarDetail().hasCMB(SiPM%CRVId::nSidesPerBar)) continue;  //don't check non-existing SiPMs
          uint16_t offlineChannel = barIndex.asUint()*4 + SiPM;
          CRVROC   onlineChannel  = crvChannelMap.online(offlineChannel);
          uint16_t feb            = onlineChannel.FEB();
          if(feb>=_nFEBs) _nFEBs=feb+1;
        }
      }

      _recoPEs          = new float[_nFEBs*CRVId::nChanPerFEB];
      _recoTime         = new float[_nFEBs*CRVId::nChanPerFEB];
      _fitStatus        = new int[_nFEBs*CRVId::nChanPerFEB];
      _depositedEnergy  = new float[_nFEBs*CRVId::nChanPerFEB];
      _coincidencePDGid = new int[_nSectorTypes];
      _coincidenceTime  = new float[_nSectorTypes];
      _coincidencePosX  = new float[_nSectorTypes];
      _coincidencePosY  = new float[_nSectorTypes];
      _coincidencePosZ  = new float[_nSectorTypes];

      _trackSlope       = new float[_nSectorTypes];
      _trackIntercept   = new float[_nSectorTypes];
      _trackPoints      = new int[_nSectorTypes];
      _trackPEs         = new float[_nSectorTypes];
      _trackChi2        = new float[_nSectorTypes];

      _tree->Branch("runNumber",&_runNumber,"runNumber/I");
      _tree->Branch("subrunNumber",&_subrunNumber,"subrunNumber/I");
      _tree->Branch("spillNumber",&_spillNumber,"spillNumber/I");
      _tree->Branch("spillIndex",&_spillIndex,"spillIndex/I");
      _tree->Branch("eventNumber",&_eventNumber,"eventNumber/I");
      _tree->Branch("PEs",_recoPEs,Form("PEs[%i][%i]/F",_nFEBs,(unsigned int)CRVId::nChanPerFEB));
      _tree->Branch("time",_recoTime,Form("time[%i][%i]/F",_nFEBs,(unsigned int)CRVId::nChanPerFEB));
      _tree->Branch("fitStatus",_fitStatus,Form("fitStatus[%i][%i]/I",_nFEBs,(unsigned int)CRVId::nChanPerFEB));
      _tree->Branch("depositedEnergy",_depositedEnergy,Form("depositedEnergy[%i][%i]/F",_nFEBs,(unsigned int)CRVId::nChanPerFEB));
      _tree->Branch("coincidencePDGid",_coincidencePDGid,Form("coincidencePDGid[%i]/I",_nSectorTypes));
      _tree->Branch("coincidenceTime",_coincidenceTime,Form("coincidenceTime[%i]/F",_nSectorTypes));
      _tree->Branch("coincidencePosX",_coincidencePosX,Form("coincidencePosX[%i]/F",_nSectorTypes));
      _tree->Branch("coincidencePosY",_coincidencePosY,Form("coincidencePosY[%i]/F",_nSectorTypes));
      _tree->Branch("coincidencePosZ",_coincidencePosZ,Form("coincidencePosZ[%i]/F",_nSectorTypes));
      _tree->Branch("trackSlope", _trackSlope, Form("trackSlope[%i]/F",_nSectorTypes));
      _tree->Branch("trackIntercept", _trackIntercept, Form("trackIntercept[%i]/F",_nSectorTypes));
      _tree->Branch("trackPoints", _trackPoints, Form("trackPoints[%i]/I",_nSectorTypes));
      _tree->Branch("trackPEs", _trackPEs, Form("trackPEs[%i]/F",_nSectorTypes));
      _tree->Branch("trackChi2", _trackChi2, Form("trackChi2[%i]/F",_nSectorTypes));

      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir = tfs->mkdir("plots");

      _histPEs.resize(_nFEBs*CRVId::nChanPerFEB);
      for(int iFEB=0; iFEB<_nFEBs; ++iFEB)
      for(int iChannel=0; iChannel<(int)CRVId::nChanPerFEB; ++iChannel)
      {
        int index=iFEB*CRVId::nChanPerFEB+iChannel;
        _histPEs[index] = tfdir.make<TH1F>(Form("PEs_%i_%i",iFEB,iChannel),
                                           Form("PE Distribution FEB %i Channel %i;PE;Counts",iFEB,iChannel),
                                           75,0,150);
      }
    }//setup at first event

    //fits for the entire stack of modules/sectors (sectorType=0) and for individual modules/sectors (sectorType=1,...)
    for(int iSectorType=0; iSectorType<_nSectorTypes; ++iSectorType)
    {
     //initialize track variables
     float sumX     =0;
     float sumY     =0;
     float sumXY    =0;
     float sumYY    =0;
     _trackSlope[iSectorType]    =0;
     _trackIntercept[iSectorType]=0;
     _trackPEs[iSectorType]      =0;
     _trackPoints[iSectorType]   =0;
     _trackChi2[iSectorType]     =-1;

      int widthDirection = counters.at(0)->getBarDetail().getWidthDirection();  //assumes that all counters are oriented in the same way
      int thicknessDirection = counters.at(0)->getBarDetail().getThicknessDirection();

      //loop through all counters
      for(iter=counters.begin(); iter!=counters.end(); iter++)
      {
        //get counter properties
        const CRSScintillatorBar &crvCounter = **iter;
        const CRSScintillatorBarIndex &barIndex = crvCounter.index();

        int sectorNumber = crvCounter.id().getShieldNumber();
        int sectorType = CRS->getCRSScintillatorShields().at(sectorNumber).getSectorType();
        if(iSectorType>0 && iSectorType!=sectorType) continue;


        CLHEP::Hep3Vector crvCounterPos = crvCounter.getPosition();
        if(widthDirection!=crvCounter.getBarDetail().getWidthDirection() ||
           thicknessDirection!=crvCounter.getBarDetail().getThicknessDirection())
        throw std::logic_error("crv modules are oriented in different directions.");

        double x=crvCounterPos[widthDirection];
        double y=crvCounterPos[thicknessDirection];

        //get visible deposited energy for each counter
        double depositedEnergy=0;
        if(crvStepsCollection.isValid())
        {
          for(size_t istep=0; istep<crvStepsCollection->size(); istep++)
          {
            CrvStep const& step(crvStepsCollection->at(istep));
            if(step.barIndex()==barIndex) depositedEnergy+=step.visibleEDep();
          }
        }

        //get reco PEs of each SiPM of this counter
        float counterPEs=0;
        for(int SiPM=0; SiPM<4; SiPM++)
        {
          if(!(*iter)->getBarDetail().hasCMB(SiPM%CRVId::nSidesPerBar)) continue;  //don't check non-existing SiPMs

          float recoPEs=0;
          float recoTime=-1;
          int   fitStatus=0;
          for(size_t recoPulseIndex=0; recoPulseIndex<crvRecoPulseCollection->size(); recoPulseIndex++)
          {
            const CrvRecoPulse &crvRecoPulse = crvRecoPulseCollection->at(recoPulseIndex);
            if(crvRecoPulse.GetScintillatorBarIndex()==barIndex && crvRecoPulse.GetSiPMNumber()==SiPM)
            {
              if(recoPEs<crvRecoPulse.GetPEs())  //record the largest pulse to remove noise hits, after pulses, ...
              {
                recoPEs   = crvRecoPulse.GetPEs();
                recoTime  = crvRecoPulse.GetPulseTime();
                fitStatus = 1;
                if(crvRecoPulse.GetRecoPulseFlags().test(CrvRecoPulseFlagEnums::failedFit)) fitStatus=2;
              }
            }
          }

          if(iSectorType==0)
          {
            uint16_t offlineChannel = barIndex.asUint()*4 + SiPM;
            CRVROC   onlineChannel  = crvChannelMap.online(offlineChannel);
            uint16_t feb            = onlineChannel.FEB();
            uint16_t febChannel     = onlineChannel.FEBchannel();

            int index=feb*CRVId::nChanPerFEB+febChannel;
            _recoPEs[index]   = recoPEs;
            _recoTime[index]  = recoTime;
            _fitStatus[index] = fitStatus;
            _depositedEnergy[index] = depositedEnergy;
            _histPEs[index]->Fill(recoPEs);
          }

          counterPEs+=recoPEs;
        }

        //update track variables for track fit
        if(counterPEs<_minTrackFitPEs) continue;

        sumX +=x*counterPEs;
        sumY +=y*counterPEs;
        sumXY+=x*y*counterPEs;
        sumYY+=y*y*counterPEs;
        _trackPEs[iSectorType]+=counterPEs;
        ++_trackPoints[iSectorType];
      }

      //track fit
      if(_trackPEs[iSectorType]>=2*_minTrackFitPEs && _trackPoints[iSectorType]>1)
      {
        if(_trackPEs[iSectorType]*sumYY-sumY*sumY!=0)
        {
          _trackSlope[iSectorType]=(_trackPEs[iSectorType]*sumXY-sumX*sumY)/(_trackPEs[iSectorType]*sumYY-sumY*sumY);
          _trackIntercept[iSectorType]=(sumX-_trackSlope[iSectorType]*sumY)/_trackPEs[iSectorType];
        }
      }

      //calculate chi2
      for(iter=counters.begin(); iter!=counters.end(); iter++)
      {
        //get counter properties
        const CRSScintillatorBar &crvCounter = **iter;
        const CRSScintillatorBarIndex &barIndex = crvCounter.index();

        CLHEP::Hep3Vector crvCounterPos = crvCounter.getPosition();
        if(widthDirection!=crvCounter.getBarDetail().getWidthDirection() ||
           thicknessDirection!=crvCounter.getBarDetail().getThicknessDirection())
        throw std::logic_error("crv modules are oriented in different directions.");

        double x=crvCounterPos[widthDirection];
        double y=crvCounterPos[thicknessDirection];

        //get reco PEs of each SiPM of this counter
        float counterPEs=0;
        for(int SiPM=0; SiPM<4; SiPM++)
        {
          if(!(*iter)->getBarDetail().hasCMB(SiPM%CRVId::nSidesPerBar)) continue;  //don't check non-existing SiPMs

          float recoPEs=0;
          for(size_t recoPulseIndex=0; recoPulseIndex<crvRecoPulseCollection->size(); recoPulseIndex++)
          {
            const CrvRecoPulse &crvRecoPulse = crvRecoPulseCollection->at(recoPulseIndex);
            if(crvRecoPulse.GetScintillatorBarIndex()==barIndex && crvRecoPulse.GetSiPMNumber()==SiPM)
            {
              if(recoPEs<crvRecoPulse.GetPEs())  //record the largest pulse to remove noise hits, after pulses, ...
              {
                recoPEs = crvRecoPulse.GetPEs();
              }
            }
          }
          counterPEs+=recoPEs;
        }
        if(counterPEs<_minTrackFitPEs) continue;

        float xFit = _trackSlope[iSectorType]*y + _trackIntercept[iSectorType];
        _trackChi2[iSectorType]+=(xFit-x)*(xFit-x)*counterPEs;  //PE-weighted chi2
      }
      _trackChi2[iSectorType]/=_trackPEs[iSectorType];
    }//track fits

    //get PDG ID of all sectors (if available)
    std::vector<float> totalEnergyDeposited(_nSectorTypes,0);
    for(int i=0; i<_nSectorTypes; ++i)
    {
      _coincidencePDGid[i]=0;
      _coincidenceTime[i]=0;
      _coincidencePosX[i]=0;
      _coincidencePosY[i]=0;
      _coincidencePosZ[i]=0;
    }
    for(size_t iCluster=0; iCluster<crvCoincidenceClusterCollection->size(); ++iCluster)
    {
      const CrvCoincidenceCluster   &cluster   = crvCoincidenceClusterCollection->at(iCluster);
      const CrvCoincidenceClusterMC &clusterMC = crvCoincidenceClusterMCCollection->at(iCluster);
      int sectorType = cluster.GetCrvSectorType();

      if(totalEnergyDeposited[sectorType]<clusterMC.GetTotalEnergyDeposited())
      {
        totalEnergyDeposited[sectorType]=clusterMC.GetTotalEnergyDeposited();
        _coincidencePDGid[sectorType]=clusterMC.GetMostLikelySimParticle()->pdgId();
        _coincidenceTime[sectorType] =clusterMC.GetAvgHitTime();
        _coincidencePosX[sectorType] =clusterMC.GetAvgHitPos().x();
        _coincidencePosY[sectorType] =clusterMC.GetAvgHitPos().y();
        _coincidencePosZ[sectorType] =clusterMC.GetAvgHitPos().z();
      }

    }

    _tree->Fill();

  } // end analyze

} // end namespace mu2e

using mu2e::CrvWidebandTest;
DEFINE_ART_MODULE(CrvWidebandTest)
