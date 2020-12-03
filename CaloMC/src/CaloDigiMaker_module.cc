//
// Simulate the readout waveform for each sensors from CaloShowerROs.
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
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "Mu2eUtilities/inc/CaloPulseShape.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/CalorimeterCalibrations.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "MCDataProducts/inc/CaloShowerRO.hh"
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


  class CaloDigiMaker : public art::EDProducer 
  {
     public:
         struct Config 
         {
             using Name    = fhicl::Name;
             using Comment = fhicl::Comment;
             fhicl::Atom<art::InputTag> caloShowerCollection { Name("caloShowerROCollection"), Comment("CaloShowerRO collection name") }; 
             fhicl::Atom<double>        blindTime            { Name("blindTime"),              Comment("Microbunch blind time") }; 
             fhicl::Atom<bool>          addNoise             { Name("addNoise"),               Comment("Add noise to waveform") }; 
             fhicl::Atom<double>        elecNphotPerNs       { Name("elecNphotPerNs"),         Comment("Electronics noise number of PE / ns ") }; 
             fhicl::Atom<double>        rinNphotPerNs        { Name("rinNphotPerNs"),          Comment("RIN noise number of PE / ns ") }; 
             fhicl::Atom<double>        darkNphotPerNs       { Name("darkNphotPerNs"),         Comment("SiPM Dark noise number of PE / ns ") }; 
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
         
         explicit CaloDigiMaker(const art::EDProducer::Table<Config>& config) :
            EDProducer{config},
            caloShowerToken_{consumes<CaloShowerROCollection>(config().caloShowerCollection())},
            blindTime_         (config().blindTime()),
            addNoise_          (config().addNoise()),
            noiseRinDark_      (config().rinNphotPerNs() + config().darkNphotPerNs()),
            noiseElec_         (config().elecNphotPerNs()),
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
         void beginJob() override;
         void beginRun(art::Run& aRun) override;

    private:       
       void generateNoise     (std::vector<std::vector<double>>&, std::vector<int>&, const ConditionsHandle<CalorimeterCalibrations>&);
       void generateNoiseFast (std::vector<std::vector<double>>&, std::vector<int>&, const ConditionsHandle<CalorimeterCalibrations>&);
       void makeDigitization  (const CaloShowerROCollection&, CaloDigiCollection&);
       void fillWaveforms     (std::vector<std::vector<double>>&, const CaloShowerROCollection&, const ConditionsHandle<CalorimeterCalibrations>&);
       void buildOutputDigi   (std::vector<std::vector<double>>&, std::vector<int>&, CaloDigiCollection&);
       void diag0             (int, const std::vector<int>&);
       void diag1             (int, double, size_t, const std::vector<int>&, int);
       void diag2             (const CaloDigiCollection&);
       void plotWF            (const std::vector<int>& waveform,    const std::string& pname, int pedestal);
       void plotWF            (const std::vector<double>& waveform, const std::string& pname, int pedestal);
       void plotNoise         (const std::vector<std::vector<double>>& waveforms, const std::vector<int>& pedestals);

       
       const art::ProductToken<CaloShowerROCollection> caloShowerToken_;
       double                  blindTime_;
       double                  mbtime_;
       bool                    addNoise_;
       double                  noiseRinDark_;
       double                  noiseElec_;
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
       
       TH1F*                   hEtot_;
       TH1F*                   hPEtot_;
  };


  //-----------------------------------------------
  void CaloDigiMaker::beginJob()
  {      
      if (diagLevel_ > 1)
      {
          art::ServiceHandle<art::TFileService> tfs;
          hEtot_   = tfs->make<TH1F>("hEtot",      "Total E dep",   150,     0,   150);
          hPEtot_  = tfs->make<TH1F>("hPEtot",     "Total E dep",   150,     0,   15000);
      }
  }

  //-----------------------------------------------------------------------------
  void CaloDigiMaker::beginRun(art::Run& aRun)
  {
      pulseShape_.buildShapes();
  }



  //---------------------------------------------------------
  void CaloDigiMaker::produce(art::Event& event)
  {

      if ( diagLevel_ > 0 ) std::cout<<"[CaloDigiMaker::produce] begin" << std::endl;

      ConditionsHandle<AcceleratorParams> accPar("ignored");
      mbtime_ = accPar->deBuncherPeriod;
      
      auto caloShowerStepHandle = event.getValidHandle(caloShowerToken_);
      const auto& CaloShowerROs = *caloShowerStepHandle;
      
      auto caloDigiColl             = std::make_unique<CaloDigiCollection>();
     
      makeDigitization(CaloShowerROs, *caloDigiColl);

      event.put(std::move(caloDigiColl));

      if ( diagLevel_ > 0 ) std::cout<<"[CaloDigiMaker::produce] end" << std::endl;    
  }

  
  //-----------------------------------------------------------------------------------------------------------------------------
  void CaloDigiMaker::makeDigitization(const CaloShowerROCollection& CaloShowerROs,CaloDigiCollection& caloDigiColl)
  {
      mu2e::GeomHandle<mu2e::Calorimeter> ch;
      calorimeter_ = ch.get();

      ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");
      if (calorimeter_->nCrystal()<1 || calorimeter_->caloInfo().getInt("nSiPMPerCrystal")<1) return;
      
      unsigned nWaveforms   = calorimeter_->nCrystal()*calorimeter_->caloInfo().getInt("nSiPMPerCrystal");
      unsigned waveformSize = (mbtime_ - blindTime_ + endTimeBuffer_) / digiSampling_;     
      std::vector<std::vector<double>> waveforms(nWaveforms,std::vector<double>(waveformSize,0.0));
      std::vector<int>                 pedestals(nWaveforms,0);
      
      if (addNoise_)
      {
          if (noiseRinDark_<1.0) generateNoise(waveforms,pedestals,calorimeterCalibrations);
          else                   generateNoiseFast(waveforms,pedestals,calorimeterCalibrations);
      }

      fillWaveforms(waveforms,CaloShowerROs,calorimeterCalibrations);
      buildOutputDigi(waveforms, pedestals, caloDigiColl);
      
      if (diagLevel_ > 1) diag2(caloDigiColl);
  }


  //------------------------------------------------------------------------------------------------------------------
  void CaloDigiMaker::generateNoise(std::vector<std::vector<double>>& waveforms, std::vector<int>& pedestals, 
                                         const ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations)
  {
     
      std::vector<double> wf            = pulseShape_.digitizedPulse(0.0);
      int                 pulseSize     = wf.size();     
      double              rinBufferTime = 0.75*pulseSize*digiSampling_;
      double              totalTime     = mbtime_ - blindTime_ + endTimeBuffer_ + rinBufferTime;
      int                 nPhMean       = int(totalTime*noiseRinDark_);
      int                 nbuffer       = int(rinBufferTime/digiSampling_);
      int                 nbins         = waveforms[0].size();

      for (unsigned i=0; i<waveforms.size(); ++i)
      {
          auto& waveform     = waveforms[i];
          double scaleFactor = MeVToADC_/calorimeterCalibrations->peMeV(i);

          //Generate the radiation induced noise (RIN)
          int nPh = randPoisson_(nPhMean);
          for (int i=0;i<nPh;++i)
          {
	      double t0 = randFlat_.fire(0.0,totalTime);       
	      wf = pulseShape_.digitizedPulse(t0);

	      int i0 = int(t0/digiSampling_) - nbuffer;
	      int l0 = (i0<0) ? -i0 : 0;      
	      int l1 = std::min(pulseSize,nbins-i0);
	      for (int l=l0;l<l1;++l) waveform[i0+l] += wf[l]*scaleFactor;
          }

          //add electronics noise        
          double noiseADC = noiseElec_*digiSampling_*scaleFactor;
          for (auto& val : waveform) val += randGauss_.fire(0.0,noiseADC);
 
          //estimate pedestal for this waveform - set it to theoretical value for the time
          pedestals[i] = int(noiseRinDark_*digiSampling_*std::accumulate(wf.begin(),wf.end(),0.0)*scaleFactor);      
      }

      if (diagLevel_==99) plotNoise(waveforms,pedestals);
  }


  //-------------------------------------------------------------------------------------------------------------------------
  // This approximation  divides the digitization period into 5 fixed intervals to generate each waveform.
  void CaloDigiMaker::generateNoiseFast(std::vector<std::vector<double>>& waveforms, std::vector<int>& pedestals, 
                                        const ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations)
  {
     
      std::vector<double> wf0 = pulseShape_.digitizedPulse(0.0);
      std::vector<double> wf1 = pulseShape_.digitizedPulse(0.2*digiSampling_);
      std::vector<double> wf2 = pulseShape_.digitizedPulse(0.4*digiSampling_);
      std::vector<double> wf3 = pulseShape_.digitizedPulse(0.6*digiSampling_);
      std::vector<double> wf4 = pulseShape_.digitizedPulse(0.8*digiSampling_);

      unsigned   pulseSize     = wf0.size();
      double     rinBufferTime = 0.75*pulseSize*digiSampling_;
      unsigned   nbuffer       = unsigned(rinBufferTime/digiSampling_);
      unsigned   nbins         = waveforms[0].size();

      for (unsigned i=0; i<waveforms.size(); ++i)
      {
          auto& waveform     = waveforms[i];
          double scaleFactor = MeVToADC_/calorimeterCalibrations->peMeV(i);

          //Generate approximation of the radiation induced noise (RIN)         
          for (unsigned i=0;i<nbuffer+nbins; ++i)
          {
	      unsigned pulseStart = (i>nbuffer) ? 0 : nbuffer-i;
	      unsigned pulseEnd   = std::min(pulseSize,nbins+nbuffer-i);

	      int nP0 = randPoisson_(noiseRinDark_);
	      int nP1 = randPoisson_(noiseRinDark_);
	      int nP2 = randPoisson_(noiseRinDark_);
	      int nP3 = randPoisson_(noiseRinDark_);
	      int nP4 = randPoisson_(noiseRinDark_);

	      for (unsigned l=pulseStart;l<pulseEnd;++l)
              {
                 waveform[i-nbuffer+l] += (nP0*wf0[l]+nP1*wf1[l]+nP2*wf2[l]+nP3*wf3[l]+nP4*wf4[l])*scaleFactor;
	      }
          }
          //add electronics noise        
          double noiseADC = noiseElec_*digiSampling_*scaleFactor;
          for (auto& val : waveform) val += randGauss_.fire(0.0,noiseADC);

          //estimate pedestal for this waveform - set it to theoretical value for the time being
          pedestals[i] = int(noiseRinDark_*digiSampling_*std::accumulate(wf0.begin(),wf0.end(),0.0)*scaleFactor);
      }

      if (diagLevel_==99) plotNoise(waveforms,pedestals);
  }



  //----------------------------------------------------------------------------------------------------------
  void CaloDigiMaker::fillWaveforms(std::vector<std::vector<double>>& waveforms, const CaloShowerROCollection& CaloShowerROs,
                                    const ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations)
  {
      unsigned totalPE(0);
      for (const auto& CaloShowerRO : CaloShowerROs)
      {
          int SiPMID = CaloShowerRO.SiPMID();
          auto& waveform = waveforms.at(SiPMID);
          for (const float PEtime : CaloShowerRO.PETime())
          {        
              float       time           = PEtime - blindTime_;         
              int         startSample    = time/digiSampling_;
              const auto& pulse          = pulseShape_.digitizedPulse(time);
              int         stopSample     = std::min(startSample+pulse.size(), waveform.size());
              
              for (int timeSample = startSample; timeSample < stopSample; ++timeSample) 
                 waveform.at(timeSample) += pulse.at(timeSample - startSample)/calorimeterCalibrations->peMeV(SiPMID)*MeVToADC_;
              
          }
          totalPE += CaloShowerRO.NPE();
      }
            
      if (diagLevel_ > 1) 
      {
         hEtot_->Fill(totalPE/calorimeterCalibrations->peMeV(0)/2.0);
         hPEtot_->Fill(totalPE);         
         std::cout<<"[CaloDigiMaker::fillWaveforms] total PE processed "<<totalPE<<std::endl;
      }
  }



  //-------------------------------------------------------------------------------------------------------------------
  void CaloDigiMaker::buildOutputDigi(std::vector<std::vector<double>>& waveforms, std::vector<int>& pedestals, 
                                      CaloDigiCollection& caloDigiColl)
  {
      float totEdepEq(0);
      std::vector<int> hitStarts{}, hitStops{}, wf(waveforms[0].size(),0);
      for (unsigned int iRO=0; iRO<waveforms.size(); ++iRO)
      {    
           // round the waveform into integers and apply maxADC cut
           int wfSize = waveforms[iRO].size();
           for (int i=0;i<wfSize;++i) wf[i] = std::min(maxADCCounts_,int(waveforms[iRO][i]) - pedestals[iRO]);

           if (diagLevel_ > 2) diag0(iRO, wf);

	   hitStarts.clear();
	   hitStops.clear();
           int timeSample(nBinsPeak_);
           while (timeSample < wfSize-nBinsPeak_) // chech if bound is correct.
	   {
	       // find the local maximum over a window of size peakWindow_ and above threshold value
	       if (wf[timeSample] < minPeakADC_) {++timeSample; continue;}
               
               auto it1 =  wf.begin()+timeSample-nBinsPeak_, it2=wf.begin()+timeSample+nBinsPeak_+1;
               if (std::max_element(it1,it2) != wf.begin()+timeSample) {++timeSample; continue;}

	       // find the starting / stopping point of the peak
	       // the stopping point is the first value below the threshold
	       int sampleStart = std::max(timeSample - bufferDigi_, 0);
	       int sampleStop  = timeSample;
	       while (sampleStop < wfSize-1)
	       {
		   if (wf[sampleStop]< minPeakADC_) break;
		   ++sampleStop;
	       }

	       //if (sampleStop-sampleStart<2) continue;
               if (diagLevel_ > 2) std::cout<<"[CaloDigiMaker] found peak with startSample="<<sampleStart<<"  stopSample="<<sampleStop
                                            <<"  timePeak="<<timeSample<<std::endl;

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

	        auto it1 = wf.begin()+sampleStart, it2 = wf.begin()+sampleStop+1;
                auto peakPosition = std::max_element(it1,it2) - it1;
	        if (diagLevel_ >2) std::cout<<"[CaloDigiMaker] Start=" << sampleStart << " Stop=" << sampleStop << " peak in position " << peakPosition << std::endl; 

	        std::vector<int> wfsample{};
	        for (int i=sampleStart; i<=sampleStop; ++i) wfsample.push_back(std::min(int(waveforms[iRO][i]),maxADCCounts_));

	        // make the CaloDigi
	        caloDigiColl.emplace_back(CaloDigi(iRO,t0,wfsample,peakPosition) );

	        if (diagLevel_ > 2) diag1(iRO, t0, peakPosition, wfsample, pedestals[iRO]);
	        if (iRO%2==0) totEdepEq += float(wfsample[peakPosition])/MeVToADC_;
	   }
      }
      if (diagLevel_ >0) std::cout<<"[CaloDigiMaker] Total energy equivalent digitized "<<totEdepEq<<std::endl;
    }




  //-------------------------------------------------------------------------------------------------------------------
  void CaloDigiMaker::diag0(int iSiPM, const std::vector<int>& wf)
  {
      if (*std::max_element(wf.begin(),wf.end())<1) return;
      std::cout<<"CaloDigiMaker::fillOutoutRO] Waveform content for readout "<<iSiPM<<std::endl;
      for (unsigned i=0;i<wf.size();++i) {if (i%10==0 && i>0) std::cout<<"- "; std::cout<<wf[i]<<" ";}
      std::cout<<std::endl;
  }
  void CaloDigiMaker::diag1(int iSiPM, double time, size_t peakP, const std::vector<int>& wf, int ped)
  {
      std::cout<<"Created caloDigi with SiPM = "<<iSiPM<<"  t0="<<time<<" peak="<<peakP<<" and content ";
      for (const auto  &v : wf) {if (v-ped >=minPeakADC_ ) std::cout<< "**"; std::cout<<v-ped<<" ";}
      std::cout<<std::endl;
  }
  void CaloDigiMaker::diag2(const CaloDigiCollection& caloDigiColl)
  {      
      std::map<int,double> enerMap;
      for (const auto& digi : caloDigiColl) enerMap[digi.SiPMID()] += digi.waveform().at(digi.peakpos())/MeVToADC_;
      std::cout<<"[CaloDigiMaker] energy equivalent per SiPMID"<<std::endl;
      for (auto& kv : enerMap) std::cout<<" SiPMID: "<<kv.first<<"   Ener: "<<kv.second<<std::endl;
  }



  //-------------------------------------------------------------------------------------------------------------------
  void CaloDigiMaker::plotNoise(const std::vector<std::vector<double>>& waveforms, const std::vector<int>& pedestals)
  {      
       const unsigned nMaxPlot(100u);
       const unsigned nGr = waveforms.size();
       TGraph gr[nGr];

       double maxe(0);
       for (unsigned i=0;i<waveforms.size();++i)
       {
          const auto& wf = waveforms[i];
          std::vector<double> xv,yv;
          for (unsigned j=0;j<wf.size();++j) {xv.push_back(j*digiSampling_+0.5*digiSampling_);yv.push_back(wf[j]-pedestals[i]); maxe=std::max(maxe,std::abs(yv[j]));}
          gr[i] = TGraph(xv.size(),xv.data(),yv.data());
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
       for (unsigned i=0;i<std::min(nMaxPlot,nGr);++i) gr[i].Draw("L");
       c1.SaveAs(ss.str().c_str());   
   }


  //-------------------------------------------------------------------------------------------------------------------
  void CaloDigiMaker::plotWF(const std::vector<int>& waveform, const std::string& pname,int pedestal)
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
  
  void CaloDigiMaker::plotWF(const std::vector<double>& waveform, const std::string& pname, int pedestal)
  {      
      std::vector<int> v;
      for (unsigned i=0;i<waveform.size();++i) v.push_back(int(waveform[i]));
      plotWF(v,pname, pedestal);  
  }


}

DEFINE_ART_MODULE(mu2e::CaloDigiMaker);
