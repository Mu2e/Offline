//
// A module to extract number of PEs, arrival times, hit positions, etc. from the CRV waveforms
//
// 
// Original Author: Ralf Ehrlich

#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"

#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/CrvParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "MCDataProducts/inc/CrvPhotons.hh"
#include "MCDataProducts/inc/CrvSiPMCharges.hh"
#include "RecoDataProducts/inc/CrvDigiCollection.hh"
#include "RecoDataProducts/inc/CrvRecoPulse.hh"

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"

#include <string>

#include <TCanvas.h>
#include <TGraph.h>
#include <TGaxis.h>
#include <TH1.h>
#include <TMath.h>
#include <TMarker.h>
#include <TStyle.h>

namespace mu2e 
{
  class CrvPlot : public art::EDAnalyzer
  {
    public:
    explicit CrvPlot(fhicl::ParameterSet const& pset);
    void analyze(const art::Event& e);
    void beginRun(art::Run const &run);
    void beginJob();

    private:
    std::vector<int> _events;
    std::vector<int> _crvBarIndices;
    std::string _crvPhotonsModuleLabel;
    std::string _crvSiPMChargesModuleLabel;
    std::string _crvDigiModuleLabel;
    std::string _crvRecoPulsesModuleLabel;
    double      _timeStart, _timeEnd;
    double      _recoPulsePedestal;

    double      _microBunchPeriod;
    double      _digitizationPeriod;

    std::map<std::pair<int,int>, TCanvas*> _canvas;
  };

  CrvPlot::CrvPlot(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    _events(pset.get<std::vector<int> >("events")),
    _crvBarIndices(pset.get<std::vector<int> >("crvBarIndices")),
    _crvPhotonsModuleLabel(pset.get<std::string>("crvPhotonsModuleLabel")),
    _crvSiPMChargesModuleLabel(pset.get<std::string>("crvSiPMChargesModuleLabel")),
    _crvDigiModuleLabel(pset.get<std::string>("crvDigiModuleLabel")),
    _crvRecoPulsesModuleLabel(pset.get<std::string>("crvRecoPulsesModuleLabel")),
    _timeStart(pset.get<double>("timeStart")),
    _timeEnd(pset.get<double>("timeEnd"))
  {
  }

  void CrvPlot::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;
    for(size_t i=0; i<_events.size(); ++i)
    for(size_t j=0; j<_crvBarIndices.size(); ++j)
    {
      std::stringstream canvasTitle;
      canvasTitle<<"CrvPlot_Event_"<<_events[i]<<"_BarIndex_"<<_crvBarIndices[j];
      std::pair<int,int> eventBarIndex(_events[i],_crvBarIndices[j]);
      _canvas.emplace(eventBarIndex, tfs->makeAndRegister<TCanvas>(canvasTitle.str().c_str(),canvasTitle.str().c_str()));
      TCanvas *c = _canvas[eventBarIndex];
      c->SetCanvasSize(1000,1000);
      c->cd();
      c->Draw();
      c->Divide(2,2);
    }
  }

  void CrvPlot::beginRun(art::Run const &run)
  {
    mu2e::ConditionsHandle<mu2e::AcceleratorParams> accPar("ignored");
    _microBunchPeriod = accPar->deBuncherPeriod;

    mu2e::ConditionsHandle<mu2e::CrvParams> crvPar("ignored");
    _digitizationPeriod  = crvPar->digitizationPeriod;
    _recoPulsePedestal   = crvPar->pedestal;
  }

  void CrvPlot::analyze(const art::Event& event) 
  {
    int eventID = event.id().event();

    if(std::find(_events.begin(),_events.end(),eventID)==_events.end()) return;

    gStyle->SetOptStat(0);

    art::Handle<CrvPhotonsCollection> crvPhotonsCollection;
    event.getByLabel(_crvPhotonsModuleLabel,"",crvPhotonsCollection);

    art::Handle<CrvSiPMChargesCollection> crvSiPMChargesCollection;
    event.getByLabel(_crvSiPMChargesModuleLabel,"",crvSiPMChargesCollection);

    art::Handle<CrvDigiCollection> crvDigiCollection;
    event.getByLabel(_crvDigiModuleLabel,"",crvDigiCollection);

    art::Handle<CrvRecoPulseCollection> crvRecoPulseCollection;
    event.getByLabel(_crvRecoPulsesModuleLabel,"",crvRecoPulseCollection);

    GeomHandle<CosmicRayShield> CRS;
    for(size_t bi=0; bi<_crvBarIndices.size(); ++bi)
    {
      const CRSScintillatorBarIndex barIndex(_crvBarIndices[bi]);
      std::pair<int,int> eventBarIndex(eventID,_crvBarIndices[bi]);

      std::vector<TH1D*> histVector;
      std::vector<TGraph*> graphVector;
      std::vector<TGaxis*> axisVector;

//Photon Times
      std::vector<double> photonTimesAdjusted[4];
      double histMax[4]={0};
      CrvPhotonsCollection::const_iterator photons;
      for(photons=crvPhotonsCollection->begin(); photons!=crvPhotonsCollection->end(); ++photons)
      {
        if(photons->GetScintillatorBarIndex()==barIndex)
        {
          int SiPM=photons->GetSiPMNumber();
          if(SiPM<0 || SiPM>3) throw cet::exception("SIM")<<"mu2e::CrvPlot: Invalid SiPM# at event "<<eventID<<std::endl;

          const std::vector<CrvPhotons::SinglePhoton> &photonTimes = photons->GetPhotons();
          std::vector<CrvPhotons::SinglePhoton>::const_iterator timeIter;
          for(timeIter=photonTimes.begin(); timeIter!=photonTimes.end(); ++timeIter)
          {
            double time = timeIter->_time;
            time = fmod(time,_microBunchPeriod); 
            photonTimesAdjusted[SiPM].push_back(time);
          }
        }
      }

      for(int SiPM=0; SiPM<4; ++SiPM)
      {
        _canvas[eventBarIndex]->cd(SiPM+1);

        std::ostringstream s2, s3;
        s2<<"Photons_"<<eventBarIndex.first<<"_"<<eventBarIndex.second<<"_"<<SiPM;
        s3<<"Fiber: "<<SiPM/2<<",  Side: "<<SiPM%2;
        TH1D *hist=new TH1D(s2.str().c_str(),s3.str().c_str(),100,_timeStart,_timeEnd);
        histVector.push_back(hist);

        for(unsigned int i=0; i<photonTimesAdjusted[SiPM].size(); ++i)
        {
          hist->Fill(photonTimesAdjusted[SiPM][i]);
        }

        hist->SetLineColor(kBlue);
        hist->GetXaxis()->SetTitle("t [ns]");
        hist->GetYaxis()->SetTitle("Photons");
        hist->GetYaxis()->SetTitleOffset(1.0);
        hist->GetYaxis()->SetAxisColor(kBlue);
        hist->GetYaxis()->SetTitleColor(kBlue);
        hist->GetYaxis()->SetLabelColor(kBlue);
        hist->Draw();

        if(histMax[SiPM]<hist->GetMaximum()) histMax[SiPM]=hist->GetMaximum();
      }

//SiPM response
      std::vector<double> siPMtimes[4], siPMcharges[4];
      CrvSiPMChargesCollection::const_iterator charges;
      for(charges=crvSiPMChargesCollection->begin(); charges!=crvSiPMChargesCollection->end(); ++charges)
      {
        if(charges->GetScintillatorBarIndex()==barIndex)
        {
          int SiPM=charges->GetSiPMNumber();
          if(SiPM<0 || SiPM>3) throw cet::exception("SIM")<<"mu2e::CrvPlot: Invalid SiPM# at event "<<eventID<<std::endl;

          const std::vector<CrvSiPMCharges::SingleCharge> &timesAndCharges = charges->GetCharges();
          for(size_t i=0; i<timesAndCharges.size(); ++i)
          {
            siPMtimes[SiPM].push_back(timesAndCharges[i]._time);
            siPMcharges[SiPM].push_back(timesAndCharges[i]._chargeInPEs);
          }
        }
      }

      for(int SiPM=0; SiPM<4; ++SiPM)
      {
        _canvas[eventBarIndex]->cd(SiPM+1);

        double scaleSiPMResponse = 0.5;
        std::ostringstream s4;
        s4<<"SiPMcharges_"<<eventBarIndex.first<<"_"<<eventBarIndex.second<<"_"<<SiPM;
        TH1D *histSiPMResponse=new TH1D(s4.str().c_str(),"",100,_timeStart,_timeEnd);
        histVector.push_back(histSiPMResponse);
        for(unsigned int i=0; i<siPMtimes[SiPM].size(); ++i)
        {
          histSiPMResponse->Fill(siPMtimes[SiPM][i], siPMcharges[SiPM][i]*scaleSiPMResponse);
        }
        histSiPMResponse->SetLineColor(kOrange);
        histSiPMResponse->Draw("same");

        if(histMax[SiPM]==0) histMax[SiPM]=1;

        TGaxis *axisSiPMResponse = new TGaxis(_timeEnd,0,_timeEnd,histMax[SiPM],0,histMax[SiPM]/scaleSiPMResponse,10,"+L");
        axisVector.push_back(axisSiPMResponse);
        axisSiPMResponse->SetTitle("SiPM charges [PE]");
        axisSiPMResponse->SetTitleOffset(1.0);
        axisSiPMResponse->SetTitleColor(kOrange);
        axisSiPMResponse->SetLineColor(kOrange);
        axisSiPMResponse->SetLabelColor(kOrange);
        axisSiPMResponse->Draw("same");
      }

//waveforms
      std::vector<std::vector<double> > ADCs[4];  //there can be multiple disconnected ADC waveforms per SiPM
      std::vector<std::vector<double> > times[4];
      double scale[4]={NAN};
      double maxADC[4]={NAN};
      CrvDigiCollection::const_iterator digis;
      for(digis=crvDigiCollection->begin(); digis!=crvDigiCollection->end(); ++digis)
      {
        if(digis->GetScintillatorBarIndex()==barIndex) 
        {
          int SiPM=digis->GetSiPMNumber();
          if(SiPM<0 || SiPM>3) throw cet::exception("SIM")<<"mu2e::CrvPlot: Invalid SiPM# at event "<<eventID<<std::endl;

          std::vector<double> ADCsTmp;
          std::vector<double> timesTmp;
          for(size_t j=0; j<digis->GetADCs().size(); ++j)
          {    
            double ADC=digis->GetADCs()[j];
            double time=(digis->GetStartTDC()+j)*_digitizationPeriod;
            ADCsTmp.push_back(ADC);
            timesTmp.push_back(time); 
            if(maxADC[SiPM]<ADC || isnan(maxADC[SiPM])) maxADC[SiPM]=ADC;
          }
          ADCs[SiPM].push_back(ADCsTmp);
          times[SiPM].push_back(timesTmp); 
        }
      }

      for(int SiPM=0; SiPM<4; ++SiPM)
      {
        if(maxADC[SiPM]<=0) continue;  //TODO: What should be done for cases, where the baseline is below 0?
        scale[SiPM] = histMax[SiPM]/maxADC[SiPM];

        _canvas[eventBarIndex]->cd(SiPM+1);
        for(size_t i=0; i<ADCs[SiPM].size(); ++i)
        {
          std::vector<double> &ADCsTmp=ADCs[SiPM][i];
          std::vector<double> &timesTmp=times[SiPM][i];

          size_t n = ADCsTmp.size();
          if(n>0)
          {
            for(size_t j=0; j<n; ++j) ADCsTmp[j]*=scale[SiPM];

            TGraph *graphW=new TGraph(n,timesTmp.data(),ADCsTmp.data());
            graphVector.push_back(graphW);
            graphW->SetTitle("");
            graphW->SetMarkerStyle(20);
            graphW->SetMarkerSize(1.5);
            graphW->SetMarkerColor(kRed);
            graphW->Draw("sameP");
          }
        }

        TGaxis *axis = new TGaxis(_timeEnd,0,_timeEnd,histMax[SiPM],0,maxADC[SiPM],10,"-R");
        axisVector.push_back(axis);
        axis->SetTitle("ADC");
        axis->SetTitleOffset(1.5);
        axis->SetTitleColor(kRed);
        axis->SetLineColor(kRed);
        axis->SetLabelColor(kRed);
        axis->Draw("same");
      }

//fit
      CrvRecoPulseCollection::const_iterator recoPulse;
      for(recoPulse=crvRecoPulseCollection->begin(); recoPulse!=crvRecoPulseCollection->end(); ++recoPulse)
      {
        if(recoPulse->GetScintillatorBarIndex()==barIndex) 
        {
          int SiPM=recoPulse->GetSiPMNumber();
          if(SiPM<0 || SiPM>3) throw cet::exception("SIM")<<"mu2e::CrvPlot: Invalid SiPM# at event "<<eventID<<std::endl;
          if(isnan(scale[SiPM])) continue; //don't draw anything, if there are no digis
          _canvas[eventBarIndex]->cd(SiPM+1);

          double fitParam0 = recoPulse->GetPulseHeight()*TMath::E();
          double fitParam1 = recoPulse->GetPulseTime();
          double fitParam2 = recoPulse->GetPulseBeta();

          double tF1=fitParam1 - 2.0*fitParam2;
          double tF2=fitParam1 + 8.0*fitParam2;
          int nF=(tF2-tF1)/1.0 + 1;
          double *tF = new double[nF];
          double *vF = new double[nF];
          for(int iF=0; iF<nF; ++iF)
          {
            tF[iF] = tF1 + iF*1.0;
            vF[iF] = fitParam0*TMath::Exp(-(tF[iF]-fitParam1)/fitParam2-TMath::Exp(-(tF[iF]-fitParam1)/fitParam2));
            vF[iF]+=_recoPulsePedestal;
            vF[iF]*=scale[SiPM];
            if(isnan(vF[iF])) nF=0;
          }
          if(nF>0)
          {
            TGraph *graphF=new TGraph(nF,tF,vF);
            graphVector.push_back(graphF);
            graphF->SetTitle("");
            graphF->SetLineWidth(2);
            graphF->SetLineColor(kGreen);
            graphF->Draw("same");
          }

          delete[] tF;
          delete[] vF;
        }
      }
    } // end bar index
  } // end analyze

} // end namespace mu2e

using mu2e::CrvPlot;
DEFINE_ART_MODULE(CrvPlot)
