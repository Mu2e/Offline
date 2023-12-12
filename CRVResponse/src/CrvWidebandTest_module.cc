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
    ~CrvWidebandTest();
    void analyze(const art::Event& e);
    void beginJob();
    void endJob();

    private:
    std::string   _crvStepsModuleLabel;
    std::string   _crvRecoPulsesModuleLabel;
    std::string   _crvCoincidenceClusterModuleLabel;
    std::string   _crvCoincidenceClusterMCModuleLabel;
    const float   _minTrackFitPEs{5};
    art::InputTag _protonBunchTimeTag;

    int     _nFEBs;
    int     _nSectorTypes;

    int     _runNumber;
    int     _subrunNumber;
    int     _spillNumber{1};  //to have the same variables as in the Wideband files
    int     _spillIndex{1};   //to have the same variables as in the Wideband files
    int     _eventNumber;

    float  *_recoPEs{NULL};
    float  *_recoTime{NULL};
    int    *_fitStatus{NULL};
    float  *_depositedEnergy{NULL};
    int    *_coincidencePDGid{NULL};
    float  *_coincidenceTime{NULL};
    float  *_coincidencePosX{NULL};
    float  *_coincidencePosY{NULL};
    float  *_coincidencePosZ{NULL};

    float   _trackSlope;      //using slope=dx/dy to avoid inf for vertical tracks
    float   _trackIntercept;  //x value, where y=0
    int     _trackPoints;
    float   _trackPEs;
    float   _trackChi2;

    TTree  *_tree;

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
    _tree = tfs->make<TTree>("WidebandTree", "WidebandTree");
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

      _tree->Branch("runNumber",&_runNumber,"runNumber/I");
      _tree->Branch("subrunNumber",&_subrunNumber,"subrunNumber/I");
      _tree->Branch("spillNumber",&_spillNumber,"spillNumber/I");
      _tree->Branch("spillIndex",&_spillNumber,"spillIndex/I");
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
      _tree->Branch("trackSlope", &_trackSlope, "trackSlope/F");
      _tree->Branch("trackIntercept", &_trackIntercept, "trackIntercept/F");
      _tree->Branch("trackPoints", &_trackPoints, "trackPoints/I");
      _tree->Branch("trackPEs", &_trackPEs, "trackPEs/F");
      _tree->Branch("trackChi2", &_trackChi2, "trackChi2/F");
    }

    //initialize track variables
    float sumX     =0;
    float sumY     =0;
    float sumXY    =0;
    float sumYY    =0;
    _trackSlope    =0;
    _trackIntercept=0;
    _trackPEs      =0;
    _trackPoints   =0;
    _trackChi2     =-1;

    int widthDirection = counters.at(0)->getBarDetail().getWidthDirection();  //assumes that all counters are oriented in the same way
    int thicknessDirection = counters.at(0)->getBarDetail().getThicknessDirection();

    //loop through all counters
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

        uint16_t offlineChannel = barIndex.asUint()*4 + SiPM;
        CRVROC   onlineChannel  = crvChannelMap.online(offlineChannel);
        uint16_t feb            = onlineChannel.FEB();
        uint16_t febChannel     = onlineChannel.FEBchannel();

        _recoPEs[feb*CRVId::nChanPerFEB+febChannel]   = recoPEs;
        _recoTime[feb*CRVId::nChanPerFEB+febChannel]  = recoTime;
        _fitStatus[feb*CRVId::nChanPerFEB+febChannel] = fitStatus;
        _depositedEnergy[feb*CRVId::nChanPerFEB+febChannel] = depositedEnergy;

        counterPEs+=recoPEs;
      }

      //update track variables for track fit
      if(counterPEs<_minTrackFitPEs) continue;

      sumX +=x*counterPEs;
      sumY +=y*counterPEs;
      sumXY+=x*y*counterPEs;
      sumYY+=y*y*counterPEs;
      _trackPEs+=counterPEs;
      ++_trackPoints;
    }

    //track fit
    if(_trackPEs>=2*_minTrackFitPEs && _trackPoints>1)
    {
      if(_trackPEs*sumYY-sumY*sumY!=0)
      {
        _trackSlope=(_trackPEs*sumXY-sumX*sumY)/(_trackPEs*sumYY-sumY*sumY);
        _trackIntercept=(sumX-_trackSlope*sumY)/_trackPEs;
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

      float xFit = _trackSlope*y + _trackIntercept;
      _trackChi2+=(xFit-x)*(xFit-x)*counterPEs;  //PE-weighted chi2
    }
    _trackChi2/=_trackPEs;

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
