//
// A module to extract number of PEs, arrival times, hit positions, etc. from the CRV waveforms
//
//
// Original Author: Ralf Ehrlich

#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"

#include "Offline/CRVConditions/inc/CRVDigitizationPeriod.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/MCDataProducts/inc/CrvPhotons.hh"
#include "Offline/MCDataProducts/inc/CrvSiPMCharges.hh"
#include "Offline/RecoDataProducts/inc/CrvDigi.hh"
#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"
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
  class CRVTest : public art::EDAnalyzer
  {
    public:
    explicit CRVTest(fhicl::ParameterSet const& pset);
    void analyze(const art::Event& e);
    void beginJob();
    void endJob();

    private:
    std::string _crvStepsModuleLabel;
    std::string _crvPhotonsModuleLabel;
    std::string _crvSiPMChargesModuleLabel;
    std::string _crvDigiModuleLabel;
    std::string _crvRecoPulsesModuleLabel;
    art::InputTag _protonBunchTimeTag;

    TNtuple  *_recoPulses;
    TH1F     *_sipmTimes0;
    TH1F     *_sipmTimes1;
    TH1F     *_sipmTimeDiffs;
    TH1F     *_photonTimeDiffs;
    TH1F     *_adcs0;
    TH1F     *_adcs1;
  };

  CRVTest::CRVTest(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    _crvStepsModuleLabel(pset.get<std::string>("crvStepsModuleLabel")),
    _crvPhotonsModuleLabel(pset.get<std::string>("crvPhotonsModuleLabel")),
    _crvSiPMChargesModuleLabel(pset.get<std::string>("crvSiPMChargesModuleLabel")),
    _crvDigiModuleLabel(pset.get<std::string>("crvDigiModuleLabel")),
    _crvRecoPulsesModuleLabel(pset.get<std::string>("crvRecoPulsesModuleLabel")),
    _protonBunchTimeTag(pset.get<art::InputTag>("protonBunchTimeTag","EWMProducer"))
  {
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = tfs->mkdir("CrvSingleCounter");
    _recoPulses = tfdir.make<TNtuple>("RecoPulses", "RecoPulses", "event:barIndex:SiPM:nRecoPulses:recoPEs:recoPulseHeight:recoPulseWidth:recoPulseTime:recoLEtime:chi2:MCPEs:visibleEnergyDeposited");
    _sipmTimes0 = tfdir.make<TH1F>("SiPMtimes0","SiPMtimes0", 100,800,900);
    _sipmTimes1 = tfdir.make<TH1F>("SiPMtimes1","SiPMtimes1", 100,800,900);
    _sipmTimeDiffs = tfdir.make<TH1F>("SiPMtimeDiffs","SiPMtimeDiffs", 400,0,100);
    _photonTimeDiffs = tfdir.make<TH1F>("photonTimeDiffs","photonTimeDiffs", 400,0,100);
    _adcs0 = tfdir.make<TH1F>("adcs0","adcs0", 100,800,900);
    _adcs1 = tfdir.make<TH1F>("adcs1","adcs1", 100,800,900);
  }

  void CRVTest::beginJob()
  {
  }

  void CRVTest::endJob()
  {
  }

  void CRVTest::analyze(const art::Event& event)
  {
    art::Handle<CrvStepCollection> crvStepsCollection;
    event.getByLabel(_crvStepsModuleLabel,"",crvStepsCollection);

    art::Handle<CrvPhotonsCollection> crvPhotonsCollection;
    event.getByLabel(_crvPhotonsModuleLabel,"",crvPhotonsCollection);

    art::Handle<CrvSiPMChargesCollection> crvSiPMChargesCollection;
    event.getByLabel(_crvSiPMChargesModuleLabel,"",crvSiPMChargesCollection);

    art::Handle<CrvDigiCollection> crvDigiCollection;
    event.getByLabel(_crvDigiModuleLabel,"",crvDigiCollection);

    art::Handle<CrvRecoPulseCollection> crvRecoPulseCollection;
    event.getByLabel(_crvRecoPulsesModuleLabel,"",crvRecoPulseCollection);

    art::Handle<ProtonBunchTime> protonBunchTime;
    event.getByLabel(_protonBunchTimeTag, protonBunchTime);
    double TDC0time = -protonBunchTime->pbtime_;

    int eventID = event.id().event();

    GeomHandle<CosmicRayShield> CRS;
    const std::vector<std::shared_ptr<CRSScintillatorBar> > &counters = CRS->getAllCRSScintillatorBars();
    std::vector<std::shared_ptr<CRSScintillatorBar> >::const_iterator iter;
    for(iter=counters.begin(); iter!=counters.end(); iter++)
    {
      const CRSScintillatorBarIndex &barIndex = (*iter)->index();

      double visibleEnergyDeposited=0;
      if(crvStepsCollection.isValid())
      {
        for(size_t istep=0; istep<crvStepsCollection->size(); istep++)
        {
          CrvStep const& step(crvStepsCollection->at(istep));
          if(step.barIndex()==barIndex) visibleEnergyDeposited+=step.visibleEDep();
        }
      }
      double avgPhotonTime0=0;
      double avgPhotonTime1=0;
      int nPhotons0=0;
      int nPhotons1=0;

      double avgTime0=0;
      double avgTime1=0;
      double charge0=0;
      double charge1=0;
      for(int SiPM=0; SiPM<4; SiPM++)
      {

        int    nRecoPulses=0;
        int    recoPEs=0;
        double recoPulseHeight=0;
        double recoPulseBeta=0;
        double recoPulseTime=0;
        double recoLEtime=0;
        double MCPEs=0;
        double chi2=0;

        for(size_t recoPulseIndex=0; recoPulseIndex<crvRecoPulseCollection->size(); recoPulseIndex++)
        {
          const CrvRecoPulse &crvRecoPulse = crvRecoPulseCollection->at(recoPulseIndex);
          if(crvRecoPulse.GetScintillatorBarIndex()==barIndex && crvRecoPulse.GetSiPMNumber()==SiPM)
          {
            nRecoPulses++;
            if(recoPEs<crvRecoPulse.GetPEs()  //record the largest pulse to remove noise hits, after pulses, ...
               && crvRecoPulse.GetRecoPulseFlags().none())
            {
              recoPEs            = crvRecoPulse.GetPEs();
              recoPulseHeight    = crvRecoPulse.GetPulseHeight();
              recoPulseBeta      = crvRecoPulse.GetPulseBeta();
              recoPulseTime      = crvRecoPulse.GetPulseTime();
              recoLEtime         = crvRecoPulse.GetLEtime();
              chi2               = crvRecoPulse.GetPulseFitChi2();
            }
          }
        }

        for(size_t photonIndex=0; photonIndex<crvPhotonsCollection->size(); photonIndex++)
        {
          const CrvPhotons &crvPhotons = crvPhotonsCollection->at(photonIndex);
          if(crvPhotons.GetScintillatorBarIndex()==barIndex && crvPhotons.GetSiPMNumber()==SiPM)
          {
            const std::vector<CrvPhotons::SinglePhoton> &singlePhotons = crvPhotons.GetPhotons();
            for(size_t i=0; i<singlePhotons.size(); i++)
            {
              double time = singlePhotons[i]._time;
              if(SiPM==0 || SiPM==2) {avgPhotonTime0+=time; ++nPhotons0;}
              if(SiPM==1 || SiPM==3) {avgPhotonTime1+=time; ++nPhotons1;}
            }
          }
        }

        for(size_t chargeIndex=0; chargeIndex<crvSiPMChargesCollection->size(); chargeIndex++)
        {
          const CrvSiPMCharges &crvSiPMCharges = crvSiPMChargesCollection->at(chargeIndex);
          if(crvSiPMCharges.GetScintillatorBarIndex()==barIndex && crvSiPMCharges.GetSiPMNumber()==SiPM)
          {
            const std::vector<CrvSiPMCharges::SingleCharge> &singleCharges = crvSiPMCharges.GetCharges();
            for(size_t i=0; i<singleCharges.size(); i++)
            {
              double charge = singleCharges[i]._chargeInPEs;
              double time = singleCharges[i]._time;
              MCPEs+=charge;
              if(SiPM==0 || SiPM==2) _sipmTimes0->Fill(time);
              if(SiPM==1 || SiPM==3) _sipmTimes1->Fill(time);
              if(SiPM==0 || SiPM==2) {avgTime0+=time*charge; charge0+=charge;}
              if(SiPM==1 || SiPM==3) {avgTime1+=time*charge; charge1+=charge;}
            }
          }
        }

        for(size_t digiIndex=0; digiIndex<crvDigiCollection->size(); digiIndex++)
        {
          const CrvDigi &crvDigi = crvDigiCollection->at(digiIndex);
          if(crvDigi.GetScintillatorBarIndex()==barIndex && crvDigi.GetSiPMNumber()==SiPM)
          {
            uint16_t startTDC=crvDigi.GetStartTDC();
            const std::vector<int16_t> &adcs = crvDigi.GetADCs();
            for(size_t ii=0; ii<adcs.size(); ++ii)
            {
              int16_t  adc=adcs.at(ii);
              if(SiPM==0 || SiPM==2) _adcs0->Fill(TDC0time+(startTDC+ii)*CRVDigitizationPeriod,adc);
              if(SiPM==1 || SiPM==3) _adcs1->Fill(TDC0time+(startTDC+ii)*CRVDigitizationPeriod,adc);
            }
          }
        }

        _recoPulses->Fill(eventID,barIndex.asInt(),SiPM,nRecoPulses,recoPEs,recoPulseHeight,recoPulseBeta,recoPulseTime,recoLEtime,chi2,MCPEs,visibleEnergyDeposited);
      }

      avgTime0/=charge0;
      avgTime1/=charge1;
      _sipmTimeDiffs->Fill(avgTime1-avgTime0);

      avgPhotonTime0/=nPhotons0;
      avgPhotonTime1/=nPhotons1;
      _photonTimeDiffs->Fill(avgPhotonTime1-avgPhotonTime0);
    }

  } // end analyze

} // end namespace mu2e

using mu2e::CRVTest;
DEFINE_ART_MODULE(CRVTest)
