// Purpose: Filter events based on RecoNPOT to create a flat distribution for even statistics per bin when sampling
// author: H Applegate 2025

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

#include "fhiclcpp/ParameterSet.h"

#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Random/RandFlat.h"

// Mu2e includes.
#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"

// MC DataProduct
#include "Offline/MCDataProducts/inc/ProtonBunchIntensity.hh"

// ROOT
#include "TFile.h"
#include "TH1.h"

#include <iostream>
#include <string>

//using namespace std;
namespace mu2e {

  class RecoNPOTFilter : public art::EDFilter {
  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<art::InputTag>     recoNPOTTag      {Name("recoNPOTTag"         ) ,Comment("The reconstructed N(POT) estimate")};
      fhicl::Atom<float>             minWeight        {Name("minWeight"           ) ,Comment("The smallest bin height in the POT distribution to invert")};
      fhicl::Atom<std::string>       fileName         {Name("fileName"            ) ,Comment("Data file with the nominal POT distribution to assume")};
      fhicl::Atom<int>               debugLevel       {Name("debugLevel"          ) ,Comment("Debug level"), 0};
    };
    explicit RecoNPOTFilter(const art::EDFilter::Table<Config>& config);
    virtual bool filter(art::Event& event) override;


    //-----------------------------------------------------------------------------
    // helper functions
    //-----------------------------------------------------------------------------
    void findData(const art::Event& evt);

  private:

    art::RandomNumberGenerator::base_engine_t& _eng;
    CLHEP::RandFlat _randomFlat;

    // inputs
    art::InputTag   _recoNPOTTag ;
    float           _minWeight;
    std::string     _fileName;
    int             _debugLevel;


    const art::Event*                _event;
    const ProtonBunchIntensity*      _recoNPOT;
    TH1*                             _hPOT; // input POT distribution
    double                           _binWidth; // nPOT histogram bin width
    int                              _nBins; // number of POT histogram bins

  };


  RecoNPOTFilter::RecoNPOTFilter(const art::EDFilter::Table<Config>& config) :
    EDFilter{config}                                                                       ,
    _eng                   {createEngine(art::ServiceHandle<SeedService>()->getSeed())}    ,
    _randomFlat            {_eng}                                                          ,
    _recoNPOTTag           {config().recoNPOTTag()}                                        ,
    _minWeight             {config().minWeight()}                                          ,
    _fileName              {config().fileName()}                                           ,
    _debugLevel            {config().debugLevel()}
  {
    if(_minWeight <= 0.) {
      throw cet::exception("BADCONFIG") << " Minimum weight must be greater than 0 to invert!";
    }
    consumes<ProtonBunchIntensity>(_recoNPOTTag);

    // Read in the POT distribution histogram
    ConfigFileLookupPolicy configFile;
    std::string filePath = configFile(_fileName);
    TFile* f = TFile::Open(filePath.c_str(), "READ");
    if(!f) throw cet::exception("BADCONFIG") << " POT distribution file " << filePath << " not found!";
    _hPOT = (TH1*) f->Get("POT");
    if(!_hPOT) throw cet::exception("BADCONFIG") << " POT distribution in " << filePath << " file not found!";
    _hPOT->SetDirectory(0);
    f->Close();

    // assume the distribution has fixed-width binning
    _binWidth = _hPOT->GetBinWidth(1);
    _nBins = _hPOT->GetNbinsX();

    // scale the histogram to have a maximum of 1
    const double max_weight = _hPOT->GetMaximum();
    if(max_weight <= 0.) throw cet::exception("BADCONFIG") << " Maximum POT bin must be greater than 0!";
    _hPOT->Scale(1./max_weight);
  }


  //-----------------------------------------------------------------------------
  // find data
  //-----------------------------------------------------------------------------
  void RecoNPOTFilter::findData(const art::Event& evt) {

    //Reco N(POT)
    auto recoNPOTH = evt.getValidHandle<ProtonBunchIntensity>(_recoNPOTTag);
    _recoNPOT = recoNPOTH.product();
  }

  //-----------------------------------------------------------------------------
  // main filter function
  //-----------------------------------------------------------------------------
  bool RecoNPOTFilter::filter(art::Event& event) {
    findData(event);
    if (_recoNPOT) {
      const auto npotreco = _recoNPOT->intensity();
      //if recoNPOT exceeds our weights array, default its bin index to the last index
      const int bin = 1 + std::min(int(npotreco / _binWidth), _nBins-1);
      //If a bin's weight < minimum binWeight, set it to the minWeight (and therefore 100% sampling probability)
      const float weight = std::max(_minWeight, (float) _hPOT->GetBinContent(bin));
      const float inverseWeight = 1./weight;

      const float inverseWeightScaled = inverseWeight*_minWeight; //ensure between 0 and 1
      const float randomVal = _randomFlat.fire(); //random integer between 0 and 1

      //Filtering Step
      const bool passed = randomVal < inverseWeightScaled;

      if(_debugLevel > 1) printf("[RecoNPOTFilter::%s: N(POT) = %8llu, index = %2i, weight = %-5.3g, scaled inv = %-6.3g, rand = %-5.3g --> passed = %o\n",
                                 __func__, npotreco, bin, weight, inverseWeightScaled, randomVal, passed);
      return passed;
    } else {
      std::cout << "[RecoNPOTFilter] Reco N(POT) not found, returning false" << std::endl;
      return false;
    }
    return false;
  }
}

using mu2e::RecoNPOTFilter;
DEFINE_ART_MODULE(RecoNPOTFilter)
