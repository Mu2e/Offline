
// Purpose: Filter events based on RecoNPOT to create a flat distribution for even statistics per bin when sampling
// author: H Applegate 2025

//For use getting results:    processName: "S5Stn"      filterName: "lumiStream"
//Filter fcl parameters must be updated in proloc.fcl<CalPatRec> and /mu2e-trig-config/core/filters/trigSupFilters.fcl
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

#include "fhiclcpp/ParameterSet.h"
#include "art_root_io/TFileService.h"

#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Random/RandFlat.h"

// Mu2e includes.
#include "Offline/SeedService/inc/SeedService.hh"

//MC DataProduct
#include "Offline/MCDataProducts/inc/ProtonBunchIntensity.hh"


#include <iostream>
#include <string>

//using namespace std;
namespace mu2e {

  class RecoNPOTFilter : public art::EDFilter {
  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<art::InputTag>     recoNPOTTag      {Name("recoNPOTTag"         ) ,Comment("for recoNPOT")};
      fhicl::Atom<double>            minWeight        {Name("minWeight"           ) ,Comment("the largest weight we want to have a 100% sampling probability for")};
    };
    explicit RecoNPOTFilter(const art::EDFilter::Table<Config>& config);
    virtual bool filter(art::Event& event) override;


    //-----------------------------------------------------------------------------
    // helper functions
    //-----------------------------------------------------------------------------
    bool findData                                        (const art::Event& evt);




  private:

    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandFlat randomFlat_;

    //-----------------------------------------------------------------------------
    // event object labels
    //-----------------------------------------------------------------------------

    art::InputTag   _recoNPOTTag ;
    double          _minWeight;
    float           _largestInverseWeight;


    const art::Event*                  _event;
    //-----------------------------------------------------------------------------
    // collections
    //-----------------------------------------------------------------------------
    const ProtonBunchIntensity*      _recoNPOT;

    //-----------------------------------------------------------------------------
    // Data Members
    //-----------------------------------------------------------------------------

    //Hardcoded from a Histogram of nPOT with 100 bins run on mnbs2 dataset
    std::array<float, 100> _weights={1e-7, 1e-7, 0.004985, 0.019526, 0.068550, 0.154549, 0.230577, 0.341089, 0.472372, 0.530120, 0.665559, 0.755713, 0.787287, 0.863731, 0.892813, 0.929788, 0.979227, 1.000000, 0.944329, 0.985875, 0.942252, 0.959701, 0.985044, 0.970503, 0.949314, 0.891151, 0.883673, 0.847113, 0.833818, 0.840050, 0.798089, 0.753220, 0.728292, 0.704612, 0.655588, 0.688409 , 0.640216, 0.615289, 0.574574, 0.570004, 0.543831, 0.525135, 0.498961, 0.468218, 0.488990, 0.465310, 0.446199, 0.412962, 0.387619, 0.389281, 0.363939, 0.347736, 0.336934, 0.342750, 0.325717, 0.299958, 0.302036, 0.272954, 0.293311, 0.287910, 0.247611, 0.242210, 0.230162, 0.206897, 0.211467, 0.216452, 0.207312, 0.199834, 0.174491, 0.171998, 0.177815, 0.171998, 0.162027, 0.136269, 0.141670, 0.154549, 0.119651, 0.120482, 0.130868, 0.119236, 0.116743, 0.110926, 0.114250, 0.105526, 0.106356, 0.105110, 0.086415, 0.094308, 0.081429, 0.081429, 0.079767, 0.068135, 0.081014, 0.068135, 0.075613, 0.068966, 0.067304, 0.081014, 0.068135, 0.059826}; //weights of 0 are manually written as  1e-7


    //min value on x-axis from nPOT histogram
    const double _minNPOT = 0;
    //nPOT value of the maxBin / number of bins
    const double _binWidth= 1000000.0; // 100e6/100  = 1e6

  };









  RecoNPOTFilter::RecoNPOTFilter(const art::EDFilter::Table<Config>& config) :
    EDFilter{config}                                                                       ,
    eng_                   {createEngine(art::ServiceHandle<SeedService>()->getSeed())}    ,
    randomFlat_            {eng_}                                                          ,
    _recoNPOTTag           {config().recoNPOTTag()}                                        ,
    _minWeight             {config().minWeight()}


  {
    _largestInverseWeight = 1./_minWeight;
  }


  //-----------------------------------------------------------------------------
  // find data
  //-----------------------------------------------------------------------------
  bool RecoNPOTFilter::findData(const art::Event& evt) {


    //recoNPOT
    auto recoNPOTH = evt.getValidHandle<ProtonBunchIntensity>(_recoNPOTTag);
    if (recoNPOTH.product() != 0){_recoNPOT = recoNPOTH.product(); }
    else {
      _recoNPOT = 0;
      std::cout << ">>> ERROR in RecoNPOTFilter::findata: ProtonBunchIntensity, recoNPOT,  not found." << std::endl;
      return false;
    }

    return true;
  }

  //-----------------------------------------------------------------------------
  // MAIN FILTER FUNCTION
  //-----------------------------------------------------------------------------
  bool RecoNPOTFilter::filter(art::Event& event) {
    bool dataexist = findData(event);
    if (dataexist) {
      if (_recoNPOT !=nullptr){
        unsigned long long npotreco = _recoNPOT->intensity();
        //if recoNPOT exceeds our weights array, default its bin index to the last index
        unsigned long long  binIndex =std::min( (_weights.size() - 1), static_cast<size_t>((npotreco - _minNPOT) / (_binWidth)) );
        //If a bin's weight < minimum binWeight, assume it has the minWeight (and therefore 100% sampling probability)
        float weight =std::max( static_cast<float>(_minWeight), static_cast<float>(_weights[binIndex]) );

        float inverseWeight = 1.0/weight; //higher weight means we want a lower sampling probability

        float inverseWeightScaled = inverseWeight/_largestInverseWeight; //ensure between 0 and 1
        float randomVal = randomFlat_.fire(); //random integer between 0 and 1

        //Filtering Step
        if(randomVal < inverseWeightScaled) {
          return true;
        }
        else { return false;}
      }
      else {
        std::cout << "[RecoNPOTFilter] recoNPOT not found, returning false" << std::endl;
        return false;
      }
    }
    else{
      std::cout << "[RecoNPOTFilter::filter] Data does not exist, returning false" << std::endl;
      return false;
    }
    return false;
  }
}

using mu2e::RecoNPOTFilter;
DEFINE_ART_MODULE(RecoNPOTFilter)
