//
// Simulate the readout waveform for each sensors from CaloshowerStepROs.
// Individual photo-electrons are generated for each readout, including photo-statistic fluctuations
// Simulate digitization procedure and produce CaloDigis. 
//
//
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include "Mu2eUtilities/inc/CaloPulseShape.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/CalorimeterCalibrations.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "MCDataProducts/inc/CaloShowerStepRO.hh"
#include "RecoDataProducts/inc/CaloDigi.hh"
#include "SeedService/inc/SeedService.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandFlat.h"

#include "TH2.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TStyle.h"
#include "TGraph.h"

#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <map>
#include <vector>
#include <numeric>


namespace mu2e {


  class CaloDigiFromShower : public art::EDProducer 
  {
     public:
         struct Config 
         {
             using Name    = fhicl::Name;
             using Comment = fhicl::Comment;
             
             fhicl::Atom<art::InputTag> caloShowerCollection { Name("caloShowerROCollection"), Comment("CaloShowerRO collection name") }; 
             fhicl::Atom<double>        blindTime            { Name("blindTime"),              Comment("Microbunch blind time") }; 
             fhicl::Atom<bool>          addNoise             { Name("addNoise"),               Comment("Add noise to waveform") }; 
             fhicl::Atom<double>        noiseElecRMS         { Name("noiseElecADC"),           Comment("Electronics noise level - ADC equivalent") }; 
             fhicl::Atom<double>        rinBufferTime        { Name("rinBufferTime"),          Comment("RIN buffer time to simulate noise") }; 
             fhicl::Atom<double>        nPhotPerDigi         { Name("nPhotPerDigi"),           Comment("number of PE / digitized bin for RIN noise") }; 
             fhicl::Atom<double>        digiSampling         { Name("digiSampling"),           Comment("Digitization time sampling") }; 
             fhicl::Atom<int>           nBits                { Name("nBits"),                  Comment("ADC Number of bits") }; 
             fhicl::Atom<int>           nBinsPeak            { Name("nBinsPeak"),              Comment("Window size for finding local maximum to digitize wf") }; 
             fhicl::Atom<double>        MeVToADC             { Name("MeVToADC"),               Comment("MeV to ADC conversion factor") }; 
             fhicl::Atom<int>           minPeakADC           { Name("minPeakADC"),             Comment("Minimum ADC hits of local peak to digitize") }; 
             fhicl::Atom<double>        endTimeBuffer        { Name("endTimeBuffer"),          Comment("Number of extra timestamps after end of pulse") }; 
             fhicl::Atom<int>           bufferDigi           { Name("bufferDigi"),             Comment("Number of timeStamps for the buffer digi") }; 
             fhicl::Atom<int>           pulseIntegralSteps   { Name("pulseIntegralSteps"),     Comment("Numer of time sub-division for CaloPulseChape") }; 
             fhicl::Atom<int>           diagLevel            { Name("diagLevel"),              Comment("Diag Level"),0 };
         };
         
         explicit CaloDigiFromShower(const art::EDProducer::Table<Config>& config) :
            EDProducer{config},
            caloShowerToken_{consumes<CaloShowerStepROCollection>(config().caloShowerCollection())},
            blindTime_         (config().blindTime()),
            addNoise_          (config().addNoise()),
            noiseElecRMS_      (config().noiseElecRMS()),
            rinBufferTime_     (config().rinBufferTime()),
            nPhotPerDigi_      (config().nPhotPerDigi()),
            digiSampling_      (config().digiSampling()),
            bufferDigi_        (config().bufferDigi()),
            nBinsPeak_         (config().nBinsPeak()),
            minPeakADC_        (config().minPeakADC()),
            MeVToADC_          (config().MeVToADC()),
            maxADCCounts_      (1 << config().nBits()),
            endTimeBuffer_     (config().endTimeBuffer()),
            pulseIntegralSteps_(config().pulseIntegralSteps()),
            diagLevel_         (config().diagLevel()),
            engine_            (createEngine(art::ServiceHandle<SeedService>()->getSeed())),
            randPoisson_       (engine_),
            randGauss_         (engine_),
            randFlat_          (engine_),
            pulseShape_        (CaloPulseShape(config().digiSampling()))
         {
             produces<CaloDigiCollection>();
         }
         
         void produce(art::Event& e)   override;
         void beginRun(art::Run& aRun) override;

    private:       
       void generateNoise       (std::vector<std::vector<double>>&, std::vector<int>&, const ConditionsHandle<CalorimeterCalibrations>&);
       void generateNoiseApprox (std::vector<std::vector<double>>&, std::vector<int>&, const ConditionsHandle<CalorimeterCalibrations>&);
       void makeDigitization    (const CaloShowerStepROCollection&, CaloDigiCollection&);
       void fillWaveforms       (std::vector<std::vector<double>>&, const CaloShowerStepROCollection&, const ConditionsHandle<CalorimeterCalibrations>&);
       void buildOutputDigi     (std::vector<std::vector<double>>&, std::vector<int>&, CaloDigiCollection&);
       void diag0               (int, const std::vector<int>&);
       void diag1               (int, double, size_t, const std::vector<int>&, int);
       void diag2               (const CaloDigiCollection&);
       void plotWF              (const std::vector<int>& waveform,    const std::string& pname, int pedestal);
       void plotWF              (const std::vector<double>& waveform, const std::string& pname, int pedestal);
       void plotNoise           (const std::vector<std::vector<double>>& waveforms, const std::vector<int>& pedestals);

       
       const art::ProductToken<CaloShowerStepROCollection> caloShowerToken_;
       double                  blindTime_;
       double                  mbtime_;
       bool                    addNoise_;
       double                  noiseElecRMS_;
       double                  rinBufferTime_;
       double                  nPhotPerDigi_;
       double                  digiSampling_;
       int                     bufferDigi_;
       int  		       nBinsPeak_;
       int  		       minPeakADC_;
       double                  MeVToADC_;
       int                     maxADCCounts_;
       double                  endTimeBuffer_;
       int                     pulseIntegralSteps_;
       int                     diagLevel_;
       CLHEP::HepRandomEngine& engine_;
       CLHEP::RandPoissonQ     randPoisson_;
       CLHEP::RandGaussQ       randGauss_;
       CLHEP::RandFlat         randFlat_;
       CaloPulseShape          pulseShape_;
       const Calorimeter*      calorimeter_;
  };


  //-----------------------------------------------------------------------------
  void CaloDigiFromShower::beginRun(art::Run& aRun)
  {
      pulseShape_.buildShapes();
  }



  //---------------------------------------------------------
  void CaloDigiFromShower::produce(art::Event& event)
  {
      if ( diagLevel_ > 0 ) std::cout<<"[CaloDigiFromShower::produce] begin" << std::endl;

      //update condition cache
      ConditionsHandle<AcceleratorParams> accPar("ignored");
      mbtime_ = accPar->deBuncherPeriod;
      
      auto caloShowerStepHandle     = event.getValidHandle(caloShowerToken_);
      const auto& caloShowerStepROs = *caloShowerStepHandle;
      auto caloDigiColl             = std::make_unique<CaloDigiCollection>();
      
      makeDigitization(caloShowerStepROs, *caloDigiColl);

      event.put(std::move(caloDigiColl));

      if ( diagLevel_ > 0 ) std::cout<<"[CaloDigiFromShower::produce] end" << std::endl;
  }

  
  //-----------------------------------------------------------------------------------------------------------------------------
  void CaloDigiFromShower::makeDigitization(const CaloShowerStepROCollection& caloShowerStepROs,CaloDigiCollection& caloDigiColl)
  {
      mu2e::GeomHandle<mu2e::Calorimeter> ch;
      calorimeter_ = ch.get();

      ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");
      if (calorimeter_->nCrystal()<1 || calorimeter_->caloInfo().nROPerCrystal()<1) return;
      
      unsigned nWaveforms   = calorimeter_->nCrystal()*calorimeter_->caloInfo().nROPerCrystal();
      unsigned waveformSize = (mbtime_ - blindTime_ + endTimeBuffer_) / digiSampling_;     
      std::vector<std::vector<double>> waveforms(nWaveforms,std::vector<double>(waveformSize,0.0));
      std::vector<int>                 pedestals(nWaveforms,0);
      
      
      if (addNoise_) generateNoise(waveforms,pedestals,calorimeterCalibrations);
      fillWaveforms(waveforms,caloShowerStepROs,calorimeterCalibrations);
      buildOutputDigi(waveforms, pedestals, caloDigiColl);
      
      if (diagLevel_ > 1) diag2(caloDigiColl);
  }


  //------------------------------------------------------------------------------------------------------------------
  void CaloDigiFromShower::generateNoise(std::vector<std::vector<double>>& waveforms, std::vector<int>& pedestals, 
                                         const ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations)
  {
     double              totalTime   = mbtime_ - blindTime_ + endTimeBuffer_ + rinBufferTime_;
     int                 nPh         = int(totalTime/digiSampling_*nPhotPerDigi_);
     int                 ibuffer     = int(rinBufferTime_/digiSampling_);
     int                 nbins       = waveforms[0].size();
     std::vector<double> wf          = pulseShape_.digitizedPulse(0.0);
     int                 pulseSize   = wf.size();
     
     for (unsigned i=0; i<waveforms.size(); ++i)
     {
         auto& waveform      = waveforms[i];
         double scaleFactor  = MeVToADC_/calorimeterCalibrations->peMeV(i);

         //Generate the radiation induced noise (RIN)
         for (int i=0;i<nPh;++i)
         {
	     double t0 = randFlat_.fire(0.0,totalTime);       
	     wf = pulseShape_.digitizedPulse(t0);

	     int i0 = int(t0/digiSampling_) - ibuffer;
	     int l0 = (i0<0) ? -i0 : 0;      
	     int l1 = std::min(pulseSize,nbins-i0);
	     for (int l=l0;l<l1;++l) waveform[i0+l] += wf[l]*scaleFactor;
         }
         
         //add electronics noise        
         for (auto& val : waveform) val += randGauss_.fire(0.0,noiseElecRMS_);
         
         //estimate pedestal for this waveform - set it to theoretical value for the time being to allow fluctuations
         //pedestals[i] = int(std::accumulate(waveform.begin(),waveform.end(),0.0)/waveform.size());
         pedestals[i] = int(digiSampling_/totalTime*nPh*std::accumulate(wf.begin(),wf.end(),0.0))*scaleFactor;
     }
     
     if (diagLevel_==99) plotNoise(waveforms,pedestals);
  }


  //-------------------------------------------------------------------------------------------------------------------------
  void CaloDigiFromShower::generateNoiseApprox(std::vector<std::vector<double>>& waveforms, std::vector<int>& pedestals, 
                                               const ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations)
  {
     
     int                 totalTime  = int((mbtime_ - blindTime_ + endTimeBuffer_ + rinBufferTime_)/digiSampling_);
     int                 nPh        = int(totalTime/digiSampling_*nPhotPerDigi_);
     int                 ibuffer    = int(rinBufferTime_/digiSampling_);
     int                 nbins      = waveforms[0].size();
     std::vector<double> wf0        = pulseShape_.digitizedPulse(0.0);
     std::vector<double> wf1        = pulseShape_.digitizedPulse(1.0);
     std::vector<double> wf2        = pulseShape_.digitizedPulse(2.0);
     std::vector<double> wf3        = pulseShape_.digitizedPulse(3.0);
     std::vector<double> wf4        = pulseShape_.digitizedPulse(4.0);
     int                 pulseSize  = wf0.size();
     

     for (unsigned i=0; i<waveforms.size(); ++i)
     {
         auto& waveform      = waveforms[i];
         double scaleFactor  = MeVToADC_/calorimeterCalibrations->peMeV(i);

         //Generate approximation of the radiation induced noise (RIN)         
         for (int i=0;i<nPh/5;++i)
         {
	     int i0 = int(randFlat_.fire(0.0,totalTime)) - ibuffer;
	     int l0 = (i0<0) ? -i0 : 0;      
	     int l1 = std::min(pulseSize,nbins-i0);
	     for (int l=l0;l<l1;++l) waveform[i0+l] += wf0[l]*scaleFactor;
	     
             i0 = int(randFlat_.fire(0.0,totalTime)) - ibuffer;
	     l0 = (i0<0) ? -i0 : 0;      
	     l1 = std::min(pulseSize,nbins-i0);
	     for (int l=l0;l<l1;++l) waveform[i0+l] += wf1[l]*scaleFactor;
	     
             i0 = int(randFlat_.fire(0.0,totalTime)) - ibuffer;
	     l0 = (i0<0) ? -i0 : 0;      
	     l1 = std::min(pulseSize,nbins-i0);
	     for (int l=l0;l<l1;++l) waveform[i0+l] += wf2[l]*scaleFactor;
	     
             i0 = int(randFlat_.fire(0.0,totalTime)) - ibuffer;
	     l0 = (i0<0) ? -i0 : 0;      
	     l1 = std::min(pulseSize,nbins-i0);
	     for (int l=l0;l<l1;++l) waveform[i0+l] += wf3[l]*scaleFactor;
	     
             i0 = int(randFlat_.fire(0.0,totalTime)) - ibuffer;
	     l0 = (i0<0) ? -i0 : 0;      
	     l1 = std::min(pulseSize,nbins-i0);
	     for (int l=l0;l<l1;++l) waveform[i0+l] += wf4[l]*scaleFactor;
         }

         //add electronics noise        
         for (auto& val : waveform) val += randGauss_.fire(0.0,noiseElecRMS_);
         
         //estimate pedestal for this waveform - set it to theoretical value for the time being to allow fluctuations
         //pedestals[i] = int(std::accumulate(waveform.begin(),waveform.end(),0.0)/waveform.size());
         pedestals[i] = int(digiSampling_/totalTime*nPh*std::accumulate(wf0.begin(),wf0.end(),0.0))*scaleFactor;
     }
     
     if (diagLevel_==99) plotNoise(waveforms,pedestals);
  }










  //----------------------------------------------------------------------------------------------------------
  void CaloDigiFromShower::fillWaveforms(std::vector<std::vector<double>>& waveforms, const CaloShowerStepROCollection& caloShowerStepROs,
                                         const ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations)
  {
      unsigned totalPE(0);
      for (const auto& caloShowerStepRO : caloShowerStepROs)
      {
          int ROID = caloShowerStepRO.ROID();
          auto& waveform = waveforms.at(ROID);
          for (const float PEtime : caloShowerStepRO.PETime())
          {        
              float       time           = PEtime - blindTime_;         
              int         startSample    = time/digiSampling_;
              const auto& pulse          = pulseShape_.digitizedPulse(time);
              int         stopSample     = std::min(startSample+pulse.size(), waveform.size());
              
              for (int timeSample = startSample; timeSample < stopSample; ++timeSample) 
                 waveform.at(timeSample) += pulse.at(timeSample - startSample)/calorimeterCalibrations->peMeV(ROID)*MeVToADC_;
              
          }
          totalPE += caloShowerStepRO.NPE();
      }
            
      if (diagLevel_ > 1) std::cout<<"[CaloDigiFromShower::fillWaveforms] total PE processed "<<totalPE<<std::endl;
  }



  //-------------------------------------------------------------------------------------------------------------------
  void CaloDigiFromShower::buildOutputDigi(std::vector<std::vector<double>>& waveforms, std::vector<int>& pedestals, 
                                           CaloDigiCollection& caloDigiColl)
  {
      float totEdepEq(0);
      std::vector<int> hitStarts, hitStops, wf(waveforms[0].size(),0);
      for (unsigned int iRO=0; iRO<waveforms.size(); ++iRO)
      {    
           // round the waveform into integers and apply maxADC cut
           int wfSize = waveforms[iRO].size();
           for (int i=0;i<wfSize;++i) wf[i] = std::min(maxADCCounts_,int(waveforms[iRO].at(i))-pedestals[iRO]);

           if (diagLevel_ > 2) diag0(iRO, wf);

	   hitStarts.clear();
	   hitStops.clear();
           int timeSample(nBinsPeak_);
           while (timeSample < wfSize-nBinsPeak_) // chech if bound is correct.
	   {
	       // find the local maximum over a window of size peakWindow_ and above threshold value
	       if (wf[timeSample] < minPeakADC_) {++timeSample; continue;}
	       if (std::max_element(&wf[timeSample-nBinsPeak_],&wf[timeSample+nBinsPeak_+1]) != &wf[timeSample]) {++timeSample; continue;}

	       // find the starting / stopping point of the peak
	       // the stopping point is the first value below the threshold

	       int sampleStart = std::max(timeSample - bufferDigi_, 0);
	       int sampleStop  = timeSample;
	       while (sampleStop < wfSize)
	       {
		      if (wf[sampleStop]< minPeakADC_) break;
		      ++sampleStop;
	       }

	       //if (sampleStop-sampleStart<2) continue;
               if (diagLevel_ > 2) std::cout<<"[CaloDigiFromShower] found peak with startSample="<<sampleStart<<"  stopSample="<<sampleStop<<"  timePeak="<<timeSample<<std::endl;

	       hitStarts.push_back(sampleStart);
	       hitStops.push_back(sampleStop);
               
	       //fast forward to end of waveform to search for next one 
	       timeSample = sampleStop+1;   
	   }
           
           
           // Concatenate peaks
           unsigned iprev(0), icurrent(1);
           while (icurrent < hitStarts.size())
           {
               if (hitStops[iprev] > hitStarts[icurrent]) {hitStops[iprev]=std::max(hitStops[icurrent],hitStops[iprev]); hitStarts[icurrent]=hitStops[icurrent]=-1;}
               else                                       {iprev = icurrent;}
               ++icurrent;
           }

	   
           // Build digi for concatenated hits   
	   for (unsigned ihit=0;ihit<hitStarts.size();ihit++)
           {
	        if (hitStarts[ihit]<0) continue;
	        int sampleStart = hitStarts[ihit];
	        int sampleStop  = hitStops[ihit];
	        int t0          = int(sampleStart*digiSampling_+ blindTime_);

	        auto peakPosition = std::max_element(&wf[sampleStart],&wf[sampleStop+1])-&wf[sampleStart];
	        if (diagLevel_ >2) std::cout<<"[CaloDigiFromShower] Start=" << sampleStart << " Stop=" << sampleStop << " peak in position " << peakPosition << std::endl; 

	        std::vector<int> wfsample;
	        for (int i=sampleStart; i<=sampleStop; ++i) wfsample.push_back(std::min(int(waveforms[iRO][i]),maxADCCounts_));

	        // make the CaloDigi
	        caloDigiColl.emplace_back(CaloDigi(iRO,t0,wfsample,peakPosition) );

	        if (diagLevel_ > 2) diag1(iRO, t0, peakPosition, wfsample, pedestals[iRO]);
	        if (iRO%2==0) totEdepEq += float(wfsample[peakPosition])/MeVToADC_;
	   }
      }
      if (diagLevel_ >0) std::cout<<"[CaloDigiFromShower] Total energy equivalent digitized "<<totEdepEq<<std::endl;
    }




  //-------------------------------------------------------------------------------------------------------------------
  void CaloDigiFromShower::diag0(int iRO, const std::vector<int>& wf)
  {
     if (*std::max_element(wf.begin(),wf.end())<1) return;
     std::cout<<"CaloDigiFromShower::fillOutoutRO] Waveform content for readout "<<iRO<<std::endl;
     for (unsigned i=0;i<wf.size();++i) {if (i%10==0 && i>0) std::cout<<"- "; std::cout<<wf[i]<<" ";}
     std::cout<<std::endl;
  }
  void CaloDigiFromShower::diag1(int iRO, double time, size_t peakP, const std::vector<int>& wf, int ped)
  {
     std::cout<<"Created caloDigi with roID = "<<iRO<<"  t0="<<time<<" peak="<<peakP<<" and content ";
     for (const auto  &v : wf) {if (v-ped >=minPeakADC_ ) std::cout<< "**"; std::cout<<v-ped<<" ";}
     std::cout<<std::endl;
  }
  void CaloDigiFromShower::diag2(const CaloDigiCollection& caloDigiColl)
  {      
     std::map<int,double> enerMap;
     for (const auto& digi : caloDigiColl) enerMap[digi.roId()] += digi.waveform().at(digi.peakpos())/MeVToADC_;
     std::cout<<"[CaloDigiFromShower] energy equivalent per RoID"<<std::endl;
     for (auto& kv : enerMap) std::cout<<" roID: "<<kv.first<<"   Ener: "<<kv.second<<std::endl;
  }



  //-------------------------------------------------------------------------------------------------------------------
  void CaloDigiFromShower::plotNoise(const std::vector<std::vector<double>>& waveforms, const std::vector<int>& pedestals)
  {      
      const unsigned nGr = waveforms.size();
      TGraph gr[nGr];
      
      double maxe(0);
      for (unsigned i=0;i<waveforms.size();++i)
      {
         const auto& wf = waveforms[i];
         double xv[9999],yv[9999];
         for (unsigned j=0;j<wf.size();++j) {xv[j]=j*digiSampling_+0.5*digiSampling_;yv[j]=wf[j]-pedestals[i]; maxe=std::max(maxe,std::abs(yv[j]));}
         gr[i] = TGraph(wf.size(),xv,yv);
         gr[i].SetLineColor(i%20+1);
      }
      
      static int nplots(0);
      std::stringstream ss;ss<<"noise_"; ss<<nplots;ss<<".pdf";
      ++nplots;

      gStyle->SetOptStat(0);
      TCanvas c1("c1","c1");
      TH1F empty("e","waveform",10,0,waveforms[0].size()*digiSampling_);
      empty.GetYaxis()->SetRangeUser(-maxe,maxe);
      empty.GetXaxis()->SetTitle("time [ns]");
      empty.GetYaxis()->SetTitle("ADC");
      empty.Draw();
      for (unsigned i=0;i<std::min(100u,nGr);++i) gr[i].Draw("L");
      c1.SaveAs(ss.str().c_str());   
   }


  //-------------------------------------------------------------------------------------------------------------------
  void CaloDigiFromShower::plotWF(const std::vector<int>& waveform, const std::string& pname,int pedestal)
  {      
     TH1F h("h","waaveform",waveform.size(),blindTime_,waveform.size()*digiSampling_+blindTime_);
     for (unsigned i=0;i<waveform.size();++i) h.SetBinContent(i,waveform[i]);
     TLine line;
     line.SetLineStyle(2);
     TLine line2;
     line2.SetLineStyle(3);
          
     gStyle->SetOptStat(0);
     TCanvas c1("c1","c1");
     h.Draw();
     line.DrawLine(blindTime_,pedestal,waveform.size()*digiSampling_+blindTime_,pedestal);
     line2.DrawLine(blindTime_,pedestal+16,waveform.size()*digiSampling_+blindTime_,pedestal+16);
     c1.SaveAs(pname.c_str());
  }
  
  void CaloDigiFromShower::plotWF(const std::vector<double>& waveform, const std::string& pname, int pedestal)
  {      
     std::vector<int> v;
     for (unsigned i=0;i<waveform.size();++i) v.push_back(int(waveform[i]));
     plotWF(v,pname, pedestal);          
  }


}

DEFINE_ART_MODULE(mu2e::CaloDigiFromShower);
