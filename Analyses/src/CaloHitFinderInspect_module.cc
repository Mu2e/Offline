//
// Recall, the caloFigiColl structure is the following:
// nTotWords - nWords_roID - roiID - nWord_roID_ihit - time_roID_ihit - Data_roID_ihit - ...
//

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Provenance.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/ConditionsService/inc/ConditionsHandle.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include "Offline/RecoDataProducts/inc/CaloRecoDigi.hh"
#include "Offline/MCDataProducts/inc/CaloShowerSim.hh"

#include "TH2F.h"
#include "TFile.h"

#include <iostream>
#include <string>


namespace mu2e {



  class CaloHitFinderInspect : public art::EDAnalyzer {

    public:

       explicit CaloHitFinderInspect(fhicl::ParameterSet const& pset);
       virtual ~CaloHitFinderInspect() {}

       virtual void beginJob();
       virtual void analyze(const art::Event& e);



    private:
       typedef std::vector<const CaloShowerSim*> CaloShowerSimVec;

       std::string  _caloDigiModuleLabel;
       std::string  _caloShowerSimModuleLabel;
       double       _digiSampling;
       int          _minDigiHitLength;

       TH1F*  _firstAmpOk;
       TH1F*  _firstAmpNo;
       TH1F*  _secondAmpOk;
       TH1F*  _secondAmpNo;
       TH1F*  _secondPSDOk;
       TH1F*  _secondPSDNo;
       TH2F*  _second2dOk;
       TH2F*  _second2dNo;
       TH1F*  _secondDTNo;
       TH1F*  _secondDTOk;

       std::map<int, std::vector<const CaloShowerSim*>> _caloShowerSimsMap;

       void   extractAmplitude(int roId, int crystalId, std::vector<int> &waveform, double adc2MeV, CaloShowerSimVec& caloShowerSims);
       void   findPeak(const std::vector<double>& x, const std::vector<double>& y, CaloShowerSimVec& caloShowerSims);
       double fitfunction(double x, double *par, int nparTot, int nparFcn);
       double logn(double x, double *par);
       double meanParabol(double x1, double x2, double x3, double y1, double y2, double y3);



  };


     CaloHitFinderInspect::CaloHitFinderInspect(fhicl::ParameterSet const& pset) :
      art::EDAnalyzer(pset),
      _caloDigiModuleLabel     (pset.get<std::string> ("caloDigiModuleLabel")),
      _caloShowerSimModuleLabel(pset.get<std::string> ("caloShowerSimModuleLabel")),
      _digiSampling            (pset.get<double>      ("digiSampling")),
      _minDigiHitLength        (pset.get<int>         ("minDigiHitLength"))
   {}




  //------------------------------------------
  void CaloHitFinderInspect::beginJob()
  {
       art::ServiceHandle<art::TFileService> tfs;
       _firstAmpOk   = tfs->make<TH1F>("firstAmpOk", "Amplitude primary ok",     200,0.,2000);
       _firstAmpNo   = tfs->make<TH1F>("firstAmpNo", "Amplitude primary not ok", 200,0.,2000);
       _secondAmpOk  = tfs->make<TH1F>("secondAmpOk","Amplitude second ok",      200,0.,200);
       _secondAmpNo  = tfs->make<TH1F>("secondAmpNo","Amplitude second not ok",  200,0.,200);
       _secondPSDOk  = tfs->make<TH1F>("secondPSDOk","PSD second ok",            100,0.,1);
       _secondPSDNo  = tfs->make<TH1F>("secondPSDNo","PSD second not ok",        100,0.,1);
       _second2dOk   = tfs->make<TH2F>("second2dOk", "AMp vs PSD second ok",     200,0., 100, 100,0.,1);
       _second2dNo   = tfs->make<TH2F>("second2dNo", "AMp vs PSD second not ok", 100,0., 100, 100,0.,1);
       _secondDTOk   = tfs->make<TH1F>("secondDTOk", "Second diff Time ok",       20,0.,20);
       _secondDTNo   = tfs->make<TH1F>("secondDTNo", "Second diff Time not ok",   20,0.,20);
  }


  //-------------------------------------------------------
  void CaloHitFinderInspect::analyze(const art::Event& event)
  {

      /*
      ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");

      art::ServiceHandle<GeometryService> geom;
      if( ! geom->hasElement<Calorimeter>() ) return;
      Calorimeter const & cal = *(GeomHandle<Calorimeter>());

      art::Handle<CaloDigiCollection> caloDigisHandle;
      event.getByLabel(_caloDigiModuleLabel, caloDigisHandle);
      CaloDigiCollection const& caloDigisColl(*caloDigisHandle);

      art::Handle<CaloShowerSimCollection> caloShowerSimHandle;
      event.getByLabel(_caloShowerSimModuleLabel, caloShowerSimHandle);
      CaloShowerSimCollection const& caloShowerSims(*caloShowerSimHandle);

      std::map<int, std::vector<const CaloShowerSim*>>  caloShowerSimsMap;
      for (auto const& caloShowerSim: caloShowerSims) caloShowerSimsMap[caloShowerSim.crystalID()].push_back(&caloShowerSim);


      for (const auto& caloDigis : caloDigisColl)
      {
         std::vector<int> caloFromDigi = caloDigis.output();
         int caloFromDigiSize = caloFromDigi.size();

         int index(1);
         while ( index < caloFromDigiSize )
         {
              int    digitizedHitLength     = caloFromDigi.at(index);
              int    roId                   = caloFromDigi.at(index+1);
              int    crystalId              = cal.caloIDMapper().crystalIDFromSiPMID(roId);
              double adc2MeV                = calorimeterCalibrations->ADC2MeV(roId);
              CaloShowerSimVec& caloShowers = caloShowerSimsMap[crystalId];

              if (digitizedHitLength > _minDigiHitLength)
              {
                  int wfSize             = digitizedHitLength - 2;
                  int wfOffset           = index+2;

                  std::vector<int> waveform(&caloFromDigi[wfOffset],&caloFromDigi[wfOffset+wfSize]);
                  extractAmplitude(roId, crystalId, waveform, adc2MeV, caloShowers);
              }

              index += digitizedHitLength;

         }
      }
      */

      return;
  }




  void CaloHitFinderInspect::extractAmplitude(int roId, int crystalId, std::vector<int>& waveform,
                                              double adc2MeV, CaloShowerSimVec& caloShowerSims)
  {


      std::vector<double> x,y;
      unsigned int index(0);

      while(index < waveform.size())
      {
          int    size    = waveform.at(index);
          double t0      = waveform.at(index+1);
          int    wfsize  = size-2;
          int    wfindex = index+2;

          x.clear();
          y.clear();
          for (int i=0;i<wfsize;++i)
          {
              x.push_back(t0 + (i+0.5)*_digiSampling);
              y.push_back(waveform.at(wfindex+i));
          }

          findPeak(x,y,caloShowerSims);

          index += size;
      }

   }


  void CaloHitFinderInspect::findPeak(const std::vector<double>& x, const std::vector<double>& y, CaloShowerSimVec& caloShowerSims)
  {

          int nparTot(0);
          double parInit[99]={0};
          std::vector<int> peakLocationInit,peakLocationRes;

          int windowPeak(2);
          int nparFcn(4);
          int minPeakAmplitude(5);

          double dummy[4]={1.0,1.0,-0.384473,11.1171};
          double peakFactor = 1.0/logn(dummy[1],dummy);

          for (unsigned int i=windowPeak; i<x.size()-windowPeak;++i)
          {
               auto maxp = std::max_element(&y[i-windowPeak],&y[i+windowPeak+1]);
               if (maxp != &y[i]) continue;
               if (y[i] > minPeakAmplitude) peakLocationInit.push_back(i);


               bool peakFound(false);
               for (const auto & showerSim : caloShowerSims) if (std::abs(showerSim->time()+20.0-x[i]) < 5.0) {peakFound=true; break;}

               if (peakFound) _firstAmpOk->Fill(y[i]);
               else           _firstAmpNo->Fill(y[i]);
          }


          for (unsigned int ipeak : peakLocationInit)
          {
               double currentAmplitudeX = fitfunction(x[ipeak],parInit,nparTot,nparFcn);

               double loc(x[ipeak]);
               if (ipeak>0 && ipeak<x.size()-1) loc = meanParabol(x[ipeak],x[ipeak-1],x[ipeak+1],y[ipeak],y[ipeak-1],y[ipeak+1]);

               parInit[nparTot++] = peakFactor*(y[ipeak] - currentAmplitudeX);
               parInit[nparTot++] = loc;
               parInit[nparTot++] = -0.384473;
               parInit[nparTot++] = 11.1171;
          }

          std::vector<double> residual;
          for (unsigned int i=0;i<x.size();++i)
          {
               double val = (y[i] > 0 ) ? y[i] - fitfunction(x[i],parInit,nparTot,nparFcn) : 0 ;
               residual.push_back(val);
          }
            for (auto ipeak : peakLocationInit) std::fill(&residual[ipeak-windowPeak],&residual[ipeak+windowPeak],0);


          for (unsigned int i=windowPeak; i<x.size()-windowPeak; ++i)
          {
               auto maxp = std::max_element(&residual[i-windowPeak], &residual[i+windowPeak+1]);
               if (maxp != &residual[i]) continue;

               bool peakFound(false);
               for (const auto & showerSim : caloShowerSims) if (std::abs(showerSim->time()+20-x[i]) < 20) {peakFound=true; break;}

               if (peakFound) _secondAmpOk->Fill(residual[i]);
               else           _secondAmpNo->Fill(residual[i]);

               if (peakFound) _secondPSDOk->Fill(residual[i]/y[i]);
               else           _secondPSDNo->Fill(residual[i]/y[i]);

               if (peakFound) _second2dOk->Fill(residual[i],residual[i]/y[i]);
               else           _second2dNo->Fill(residual[i],residual[i]/y[i]);
          }

  }





   //------------------------------------------------------
   double CaloHitFinderInspect::fitfunction(double x, double *par, int nparTot, int nparFcn)
   {
       double result(0);
       for (int i=0;i<nparTot; i+=nparFcn) result += logn(x,&par[i]);
       return result;
   }

   //-----------------------------------------------
   double CaloHitFinderInspect::logn(double x, double *par)
   {
        double logterms0 = 1.175*par[2]+sqrt(1.0+1.357225*par[2]*par[2]);
        double s0        = log(logterms0)/1.175;
        double Aterm     = par[2]/(2.506628274631*par[3]*s0);
        double logterm   = 1.0-(par[2]/par[3])*(x-par[1]);

        if (logterm<0) logterm = 1e-6;

        double expterm = log(logterm)/s0;
        return par[0]*Aterm*exp(-0.5*expterm*expterm)/1.07755263; //1.07755263 for normalization
   }

   //-----------------------------------
   double CaloHitFinderInspect::meanParabol(double x1, double x2, double x3, double y1, double y2, double y3)
   {
       double a = ((y1-y2)/(x1-x2)-(y1-y3)/(x1-x3))/(x2-x3);
       double b = (y1-y2)/(x1-x2) - a*(x1+x2);
       if (std::abs(a) < 1e-6) return (x1+x2+x3)/3.0;
       return -b/2.0/a;
   }
}

DEFINE_ART_MODULE(mu2e::CaloHitFinderInspect)
