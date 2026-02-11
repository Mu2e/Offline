#include "Offline/CaloMC/inc/CaloNoiseSimGenerator.hh"
#include "Offline/ConditionsService/inc/ConditionsHandle.hh"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "Offline/SeedService/inc/SeedService.hh"

#include "TFile.h"
#include "TH2.h"
#include "TGraph.h"
#include "TCanvas.h"

#include <algorithm>
#include <string>
#include <iostream>
#include <vector>
#include <numeric>




namespace mu2e {

   CaloNoiseSimGenerator::CaloNoiseSimGenerator(const Config& config, CLHEP::HepRandomEngine& engine, int iRO) :
     iRO_           (iRO),
     waveform_      (config.noiseWFSize(),0.0),
     pedestal_      (0.0),
     digiNoise_     (),
     digiNoiseProb_ (),
     digiSampling_  (config.digiSampling()),
     noiseRinDark_  (config.rinNphotPerNs() + config.darkNphotPerNs()),
     noiseElec_     (config.elecNphotPerNs()),
     minPeakADC_    (config.minPeakADC()),
     pePerMeV_      (config.pePerMeV()),
     MeVToADC_      (1.0/config.ADCToMeV()),
     randPoisson_   (engine),
     randGauss_     (engine),
     randFlat_      (engine),
     nMaxFragment_  (config.nMaxFragment()),
     pulseShape_    (config.pulseFileName(),config.pulseHistName(),digiSampling_),
     diagLevel_     (config.diagLevel())
   {}


   //------------------------------------------------------------------------------------------------------------------
   void CaloNoiseSimGenerator::initialize(const CaloWFExtractor& wfExtractor)
   {
       pulseShape_.buildShapes();
       generateWF(waveform_);
       generateFragments(wfExtractor);
   }

   //------------------------------------------------------------------------------------------------------------------
   void CaloNoiseSimGenerator::generateWF(std::vector<double>& wfVector)
   {
       float scaleFactor(MeVToADC_/pePerMeV_);

       std::fill(wfVector.begin(),wfVector.end(),0);

       const auto&      pulse         = pulseShape_.digitizedPulse(0.0);
       const unsigned   pulseSize     = pulse.size();
       const unsigned   bufferSize    = int(0.75*pulseSize);
       const unsigned   noiseSize     = wfVector.size();
       const double     totalTime     = (noiseSize+bufferSize)*digiSampling_;
       const int        noiseLevelPE  = int(totalTime*noiseRinDark_);

       //Generate the radiation induced noise (RIN)
       const int nPh = randPoisson_(noiseLevelPE);
       for (int i=0;i<nPh;++i)
       {
           double t0 = randFlat_.fire(0.0,totalTime);
           const auto& wf = pulseShape_.digitizedPulse(t0);

           int i0 = int(t0/digiSampling_) - bufferSize;
           int l0 = (i0<0) ? -i0 : 0;
           int l1 = std::min(pulseSize,noiseSize-i0);
           for (int l=l0;l<l1;++l) wfVector[i0+l] += wf[l]*scaleFactor;
       }

       //add electronics noise
       double noiseADC = noiseElec_*digiSampling_*scaleFactor;
       for (auto& val : wfVector) val += randGauss_.fire(0.0,noiseADC);

       //estimate pedestal for this waveform - set it to theoretical value for the time being
       pedestal_ = std::floor(noiseRinDark_*digiSampling_*std::accumulate(pulse.begin(),pulse.end(),0.0)*scaleFactor);
   }

   //------------------------------------------------------------------------------------------------------------------
   void CaloNoiseSimGenerator::generateFragments(const CaloWFExtractor& wfExtractor)
   {
       constexpr unsigned enoughFragments(20);
       unsigned nwf(0), nfound(0), length(255);
       for (;nwf<nMaxFragment_;++nwf)
       {
          std::vector<double> temp(length,0.0);

          std::vector<int> wf;
          wf.reserve(temp.size());
          for (const auto& val : temp) wf.emplace_back(val - pedestal_);

          std::vector<size_t> starts, stops;
          starts.reserve(16); stops.reserve(16);
          wfExtractor.extract(wf,starts,stops);
          if (starts.empty()) continue;

          std::vector<double> fragment;
          fragment.reserve(stops[0]-starts[0]);
          std::copy(temp.begin()+starts[0], temp.begin()+stops[0]+1, std::back_inserter(fragment));
          digiNoise_.push_back(fragment);

          ++nfound;
          if (nfound==enoughFragments) break;
       }

       digiNoiseProb_= float(nfound)/float(nwf)/float(length);
   }

   //------------------------------------------------------------------------------------------------------------------
   void CaloNoiseSimGenerator::refresh() {generateWF(waveform_);}



   //------------------------------------------------------------------------------------------------------------------
   void CaloNoiseSimGenerator::addSampleNoise(std::vector<double>& wfVector, unsigned istart, unsigned ilength)
   {
       if (ilength >=waveform_.size())
          throw cet::exception("CATEGORY")<<"[CaloNoiseSimGenerator] noise length request too long";

       unsigned irandom = unsigned(randFlat_.fire(0.,waveform_.size()-ilength));
       for (unsigned i=0;i<ilength;++i) wfVector[istart+i] += waveform_[irandom+i];

   }


   //------------------------------------------------------------------------------------------------------------------
   void CaloNoiseSimGenerator::addSaltAndPepper(std::vector<double>& wfVector)
   {
       double muNoise = waveform_.size()*digiNoiseProb_;
       int    nNoise  = randPoisson_(muNoise);
       for (int in=0;in<nNoise;++in)
       {
           unsigned idigi = unsigned(randFlat_.fire(0.,digiNoise_.size()));
           const std::vector<double>& digi = digiNoise_[idigi];
           if (wfVector.size() < digi.size()) continue;

           unsigned istart = unsigned(randFlat_.fire(0.,wfVector.size()-digi.size()));
           for (unsigned i=0;i<digi.size();++i)
           {
                if (wfVector[istart+i] < minPeakADC_) wfVector[istart+i] += digi[i];
           }
       }
   }




   //------------------------------------------------------------------------------------------------------------------
   void CaloNoiseSimGenerator::plotNoise(const std::string& name)
   {
       std::vector<double> x(waveform_.size()),y;
       std::iota(x.begin(),x.end(),0);
       TGraph gr(x.size(),x.data(),waveform_.data());
       gr.SetTitle("Original waveform");

       TH1F h1("h1","Projection waveform",100,-50,50);
       for (const auto& val: waveform_) h1.Fill(val-pedestal_);

       TCanvas c1("c1","c1");
       c1.Divide(2,2);
       c1.cd(1);
       gr.Draw("AL");
       c1.cd(2);
       h1.Draw();
       c1.SaveAs(name.c_str());
   }

}


