//
// A module to extract number of PEs, arrival times, hit positions, etc. from the CRV waveforms
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:40 $
// 
// Original Author: Ralf Ehrlich

#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"

#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "MCDataProducts/inc/CrvPhotonArrivalsCollection.hh"
#include "MCDataProducts/inc/CrvSiPMResponsesCollection.hh"
#include "MCDataProducts/inc/CrvDigiMCCollection.hh"
#include "RecoDataProducts/inc/CrvRecoPulsesCollection.hh"

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
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

    private:
    std::string _crvPhotonArrivalsModuleLabel;
    std::string _crvSiPMResponsesModuleLabel;
    std::string _crvWaveformsModuleLabel;
    std::string _crvRecoPulsesModuleLabel;
    int         _crvBarIndex;
    double      _timeStart, _timeEnd;

    double      _microBunchPeriod;

    std::vector<TCanvas*>  _canvas;
    int                    _icanvas;
  };

  CrvPlot::CrvPlot(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    _crvPhotonArrivalsModuleLabel(pset.get<std::string>("crvPhotonArrivalsModuleLabel")),
    _crvSiPMResponsesModuleLabel(pset.get<std::string>("crvSiPMResponsesModuleLabel")),
    _crvWaveformsModuleLabel(pset.get<std::string>("crvWaveformsModuleLabel")),
    _crvRecoPulsesModuleLabel(pset.get<std::string>("crvRecoPulsesModuleLabel")),
    _crvBarIndex(pset.get<int>("crvBarIndex")),
    _timeStart(pset.get<double>("timeStart")),
    _timeEnd(pset.get<double>("timeEnd"))
  {
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = tfs->mkdir("Plots");
    for(int i=0; i<20; i++)
    {
      std::stringstream canvasTitle;
      canvasTitle<<"CrvPlot_"<<i;
      _canvas.push_back(tfdir.make<TCanvas>(canvasTitle.str().c_str(),canvasTitle.str().c_str(),1000,1000));
    }
    _icanvas=-1;
  }

  void CrvPlot::beginRun(art::Run const &run)
  {
    mu2e::ConditionsHandle<mu2e::AcceleratorParams> accPar("ignored");
    _microBunchPeriod = accPar->deBuncherPeriod;
  }

  void CrvPlot::analyze(const art::Event& event) 
  {
    _icanvas++;
    if(_icanvas>=20) return;

    art::Handle<CrvPhotonArrivalsCollection> crvPhotonArrivalsCollection;
    event.getByLabel(_crvPhotonArrivalsModuleLabel,"",crvPhotonArrivalsCollection);

    art::Handle<CrvSiPMResponsesCollection> crvSiPMResponsesCollection;
    event.getByLabel(_crvSiPMResponsesModuleLabel,"",crvSiPMResponsesCollection);

    art::Handle<CrvDigiMCCollection> crvDigiMCCollection;
    event.getByLabel(_crvWaveformsModuleLabel,"",crvDigiMCCollection);

    art::Handle<CrvRecoPulsesCollection> crvRecoPulsesCollection;
    event.getByLabel(_crvRecoPulsesModuleLabel,"",crvRecoPulsesCollection);

    GeomHandle<CosmicRayShield> CRS;
    const CRSScintillatorBarIndex barIndex(_crvBarIndex);

    CrvPhotonArrivalsCollection::const_iterator  photonArrivals = crvPhotonArrivalsCollection->find(barIndex);
    CrvSiPMResponsesCollection::const_iterator   siPMResponses  = crvSiPMResponsesCollection->find(barIndex);
    CrvDigiMCCollection::const_iterator          digiMC         = crvDigiMCCollection->find(barIndex);
    CrvRecoPulsesCollection::const_iterator      recoPulses     = crvRecoPulsesCollection->find(barIndex);

    gStyle->SetOptStat(0);
    int runID = event.id().run();
    int subRunID = event.id().subRun();
    int eventID = event.id().event();
    std::stringstream canvasTitle;
    canvasTitle<<"Run "<<runID<<"  SubRun "<<subRunID<<"  Event "<<eventID;
    _canvas[_icanvas]->SetTitle(canvasTitle.str().c_str());
    _canvas[_icanvas]->cd();
    _canvas[_icanvas]->Draw();
    _canvas[_icanvas]->Divide(2,2);

    std::vector<TH1D*> histVector;
    std::vector<TGraph*> graphVector;
    std::vector<TGaxis*> axisVector;

    for(int SiPM=0; SiPM<4; SiPM++)
    {
      _canvas[_icanvas]->cd(SiPM+1);

//Photon Arrival Times
      std::vector<double> photonArrivalTimesAdjusted;
      if(photonArrivals!=crvPhotonArrivalsCollection->end())
      {
        const std::vector<CrvPhotonArrivals::SinglePhoton> &photonArrivalTimes = photonArrivals->second.GetPhotonArrivalTimes(SiPM);
        std::vector<CrvPhotonArrivals::SinglePhoton>::const_iterator timeIter;
        for(timeIter=photonArrivalTimes.begin(); timeIter!=photonArrivalTimes.end(); timeIter++)
        {
          double time = timeIter->_time;
          time = fmod(time,_microBunchPeriod); 
          photonArrivalTimesAdjusted.push_back(time);
        }
      }

      std::ostringstream s2, s3;
      s2<<"Photons_"<<_icanvas<<"__"<<SiPM;
      s3<<"Fiber: "<<SiPM/2<<",  Side: "<<SiPM%2;
      TH1D *hist=new TH1D(s2.str().c_str(),s3.str().c_str(),100,_timeStart,_timeEnd);
      histVector.push_back(hist);

      for(unsigned int i=0; i<photonArrivalTimesAdjusted.size(); i++)
      {
        hist->Fill(photonArrivalTimesAdjusted[i]);
      }

      hist->SetLineColor(kBlue);
      hist->GetXaxis()->SetTitle("t [ns]");
      hist->GetYaxis()->SetTitle("Photons");
      hist->GetYaxis()->SetTitleOffset(0.5);
      hist->GetYaxis()->SetAxisColor(kBlue);
      hist->GetYaxis()->SetTitleColor(kBlue);
      hist->GetYaxis()->SetLabelColor(kBlue);
      hist->Draw();

//SiPM response
      std::vector<double> siPMtimes, siPMcharges;
      if(siPMResponses!=crvSiPMResponsesCollection->end())
      {
        const std::vector<CrvSiPMResponses::CrvSingleSiPMResponse> &timesAndCharges = siPMResponses->second.GetSiPMResponses(SiPM);
        for(size_t i=0; i<timesAndCharges.size(); i++)
        {
          siPMtimes.push_back(timesAndCharges[i]._time);
          siPMcharges.push_back(timesAndCharges[i]._chargeInPEs);
        }
      }

      double scaleSiPMResponse = 0.5;
      TH1D *histSiPMResponse=new TH1D((s2.str()+"SiPMResponse").c_str(),"",100,_timeStart,_timeEnd);
      histVector.push_back(histSiPMResponse);
      for(unsigned int i=0; i<siPMtimes.size(); i++)
      {
        histSiPMResponse->Fill(siPMtimes[i], siPMcharges[i]*scaleSiPMResponse);
      }
      histSiPMResponse->SetLineColor(kOrange);
      histSiPMResponse->Draw("same");

//waveforms
      std::vector<double> voltages;
      std::vector<double> voltagesTimes;
      double histMax = hist->GetMaximum();
      double scale=NAN;
      if(digiMC!=crvDigiMCCollection->end())
      {
        double digitizationPrecision = digiMC->second.GetDigitizationPrecision(); //ns
        const std::vector<CrvDigiMC::CrvSingleWaveform> &singleWaveforms = digiMC->second.GetSingleWaveforms(SiPM);
        for(size_t i=0; i<singleWaveforms.size(); i++)
        {
          for(size_t j=0; j<singleWaveforms[i]._voltages.size(); j++)
          {    
            voltages.push_back(singleWaveforms[i]._voltages[j]);
            voltagesTimes.push_back(singleWaveforms[i]._startTime+digitizationPrecision*j); 
          }
        }

        unsigned int n = voltages.size();

        if(n>0)
        {
          double maxVoltage=*std::max_element(voltages.begin(),voltages.end());
          scale = histMax/maxVoltage;

          for(unsigned int i=0; i<n; i++) voltages[i]*=scale;
          TGraph *graphW=new TGraph(n,voltagesTimes.data(),voltages.data());
          graphVector.push_back(graphW);
          graphW->SetTitle("");
          graphW->SetMarkerStyle(20);
          graphW->SetMarkerSize(1.5);
          graphW->SetMarkerColor(kRed);
          graphW->Draw("sameP");
        }
      }

//fit
      if(recoPulses!=crvRecoPulsesCollection->end() && !isnan(scale))
      {
        const std::vector<CrvRecoPulses::CrvSingleRecoPulse> &pulseVector = recoPulses->second.GetRecoPulses(SiPM);
        for(unsigned int i = 0; i<pulseVector.size(); i++) 
        {
          const CrvRecoPulses::CrvSingleRecoPulse &pulse = pulseVector[i];
          double fitParam0 = pulse._pulseHeight*2.718;
          double fitParam1 = pulse._pulseTime;
          double fitParam2 = pulse._pulseWidth/1.283;

          double tF1=fitParam1 - 2.0*fitParam2;
          double tF2=fitParam1 + 8.0*fitParam2;
          int nF=(tF2-tF1)/1.0 + 1;
          double *tF = new double[nF];
          double *vF = new double[nF];
          for(int iF=0; iF<nF; iF++)
          {
            tF[iF] = tF1 + iF*1.0;
            vF[iF] = fitParam0*TMath::Exp(-(tF[iF]-fitParam1)/fitParam2-TMath::Exp(-(tF[iF]-fitParam1)/fitParam2));
            vF[iF]*=scale;
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

      if(histMax>0 && !isnan(scale) && !isnan(scaleSiPMResponse))
      {
        TGaxis *axis = new TGaxis(_timeEnd*0.93,0,_timeEnd*0.93,histMax,0,histMax/scale,10,"+L");
        axisVector.push_back(axis);
        axis->SetTitle("voltage [V]");
        axis->SetTitleOffset(-0.5);
        axis->SetTitleColor(kRed);
        axis->SetLineColor(kRed);
        axis->SetLabelColor(kRed);
        axis->Draw("same");

        TGaxis *axisSiPMResponse = new TGaxis(_timeEnd,0,_timeEnd,histMax,0,histMax/scaleSiPMResponse,10,"+L");
        axisVector.push_back(axisSiPMResponse);
        axisSiPMResponse->SetTitle("SiPM output [PE]");
        axisSiPMResponse->SetTitleOffset(1.0);
        axisSiPMResponse->SetTitleColor(kOrange);
        axisSiPMResponse->SetLineColor(kOrange);
        axisSiPMResponse->SetLabelColor(kOrange);
        axisSiPMResponse->Draw("same");
      }

    } //SiPM

    _canvas[_icanvas]->Write();

//    delete _canvas[_icanvas];

    for(size_t i=0; i<histVector.size(); i++) delete histVector[i];
    for(size_t i=0; i<graphVector.size(); i++) delete graphVector[i];
    for(size_t i=0; i<axisVector.size(); i++) delete axisVector[i];

    histVector.clear();
    graphVector.clear();
    axisVector.clear();

  } // end analyze

} // end namespace mu2e

using mu2e::CrvPlot;
DEFINE_ART_MODULE(CrvPlot)
