#include "cetlib_except/exception.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/Mu2eUtilities/inc/CaloNoiseUtil.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"

#include "TFile.h"
#include "TH1F.h"
#include "TKey.h"
#include "TDirectory.h"

#include <string>
#include <memory>
#include <vector>
#include <iostream>
#include <numeric>


namespace mu2e {

   CaloNoiseUtil::CaloNoiseUtil(const Config& config, CLHEP::HepRandomEngine& engine) :
      generate_      {config.generate()},
      fileName_      {config.histoFileName()},
      prefix_        {config.histoPrefix()},
      digiSampling_  {config.digiSampling()},
      noiseRinDark_  {config.rinNphotPerNs() + config.darkNphotPerNs()},
      noiseElec_     {config.elecNphotPerNs()},
      randPoisson_   {engine},
      randGauss_     {engine},
      randFlat_      {engine},
      pulseCache_    {config.pulseCache()},
      dumpGenerated_ {config.dumpGenerated()},
      histoBaseID_   {-1},
      pedestal_      {0},
      noiseMap_      {}
   {}


   //----------------------------------------------------------------------------------------------------------------------
   void CaloNoiseUtil::prepare(int histoID, double peToADC)
   {
       int baseID = histoID / base;
       if (baseID == histoBaseID_) return;
       if (generate_) generateCache(baseID, peToADC);
       else           fillCache(baseID);
   }

   //----------------------------------------------------------------------------------------------------------------------
   void CaloNoiseUtil::fillCache(int histoBaseID)
   {
      // Cache is for all baseID, clear it
      noiseMap_.clear();
      histoBaseID_ = histoBaseID;

      // Refill the cache with all histos sharing the same baseID
      ConfigFileLookupPolicy resolveFullPath;
      std::string fullFileName = resolveFullPath(fileName_);

      TFile file(fullFileName.c_str());
      if (!file.IsOpen()) throw cet::exception("NOISEREADER")<<"Filename"<<fullFileName.c_str()<<" does not exist\n";

      TIter next(file.GetListOfKeys());
      TKey* key;

      while ((key = (TKey*)next())) {
        TObject* obj = key->ReadObj();
        if (!obj->InheritsFrom(TH1F::Class()))
          continue;

        TH1F* histo = static_cast<TH1F*>(obj);
        std::string name = histo->GetName();

        // Parse the integer after the prefix
        std::string suffix = name.substr(prefix_.size());
        std::istringstream iss(suffix);

        int hid;
        if (!(iss >> hid) || !iss.eof())
          throw cet::exception("NOISEREADER")<<"Hitsogram "<<name.c_str()<<" is invalid\n";

        // Only take histos in the right batch to reduce memory footprint
        // the batch is defined as histogrm id/base

        if (hid/base != histoBaseID/base) continue;

        const float* array = histo->GetArray();
        noiseMap_[hid].assign(array + 1, array + histo->GetNbinsX() + 1);

        // estimate pedestal, take a single value for every waveform
        // (pedestals will be stored somewhere else later)

        float sum = 0.0;
        for (int i = 1; i <= histo->GetNbinsX(); ++i) sum += histo->GetBinContent(i);
        pedestal_ = std::trunc(sum /histo->GetNbinsX() );
      }
   }


   //----------------------------------------------------------------------------------------------------------------------
   void CaloNoiseUtil::generateCache(int histoBaseID, double peToADC)
   {
       // Cache is for all baseID, clear it
       noiseMap_.clear();
       histoBaseID_ = histoBaseID;
       constexpr unsigned noiseSize{10000};
       std::vector<float> waveform(noiseSize,0.0);

       pulseCache_.buildCache();
       const auto&        pulse         = pulseCache_.digitizedPulse(0.0);
       const unsigned     pulseSize     = pulse.size();
       const unsigned     bufferSize    = int(0.75*pulseSize);
       const double       totalTime     = (noiseSize+bufferSize)*digiSampling_;
       const int          noiseLevelPE  = int(totalTime*noiseRinDark_);

       //Generate the radiation induced noise (RIN)
       const int nPh = randPoisson_(noiseLevelPE);
       for (int i=0;i<nPh;++i) {
           double t0 = randFlat_.fire(0.0,totalTime);
           const auto& wf = pulseCache_.digitizedPulse(t0);
           int i0 = int(t0/digiSampling_) - bufferSize;
           int l0 = (i0<0) ? -i0 : 0;
           int l1 = std::min(pulseSize,noiseSize-i0);
           for (int l=l0;l<l1;++l) waveform[i0+l] += wf[l];
       }
       for (auto& wf : waveform) wf *= peToADC;

       //add electronics noise
       double noiseADC = noiseElec_*digiSampling_*peToADC;
       for (auto& val : waveform) val += randGauss_.fire(0.0,noiseADC);

       noiseMap_.emplace(histoBaseID, waveform);

       //estimate pedestal for this waveform - set it to theoretical value for the time being
       pedestal_ = std::trunc(noiseRinDark_*digiSampling_*std::accumulate(pulse.begin(),pulse.end(),0.0)*peToADC);

       if (dumpGenerated_){
         dumpNoise("caloNoise.root",waveform);
         dumpGenerated_ = false;
       }
   }


   //----------------------------------------------------------------------------------------------------------------------
   std::span<float> CaloNoiseUtil::noiseSegment(int histoID, size_t istart, size_t ilength)
   {
     int baseID = histoID/base;

     if (baseID != histoBaseID_) {
       if (generate_)  throw cet::exception("CaloNoiseUtil")
            << "noiseSegment called before prepare() for baseID " << histoID/base << "\n";
       else fillCache(baseID);
     }

     auto iter = noiseMap_.find(histoID);
     if (iter == noiseMap_.end())
       throw cet::exception("CALONOISEUTIL")<<"histoID "<<histoID<<" is invalid\n";
     auto& vec = iter->second;

     if (ilength >= vec.size())
       throw cet::exception("CALONOISEUTIL")<<"noise length request too long\n";

     size_t irandom = size_t(randFlat_.fire(0.,vec.size()-ilength));
     return std::span<float>(vec.data() + irandom, ilength);
   }

   //----------------------------------------------------------------------------------------------------------------------
   int CaloNoiseUtil::pedestal() {return pedestal_;}

   //----------------------------------------------------------------------------------------------------------------------
   void CaloNoiseUtil::printCache()
   {
     std::cout<<"CaloNoiseUtil cache\n";
     for (const auto& kv : noiseMap_) std::cout<<"Histo id "<<kv.first<<"  noise waveform length "<<kv.second.size()<<"\n";
   }

   //------------------------------------------------------------------------------------------------------------------
   void CaloNoiseUtil::dumpNoise(const std::string& fname, const std::vector<float>& wave)
   {
      TFile outfile(fname.c_str(), "RECREATE");

      TH1F h("histo_0","histo_0", wave.size(), 0, wave.size());
      for (size_t i = 0; i < wave.size(); ++i) h.SetBinContent(i+1, wave[i]);

      h.Write();
      outfile.Close();
      std::cout<<"CaloNoiseUtil written waveform in "<<fname<<"\n";
   }
}
