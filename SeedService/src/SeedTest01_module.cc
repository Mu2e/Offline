//
//  Test the SeedService.
//
//
//  Contact person Rob Kutschke
//

//
#include "SeedService/inc/SeedService.hh"

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

// C++ includes.
#include <iostream>

using namespace std;

namespace mu2e {

  class SeedTest01 : public art::EDAnalyzer {

  public:

    explicit SeedTest01(fhicl::ParameterSet const& pset);

    void analyze(const art::Event& event);

    void beginRun   ( const art::Run&    run);
    void beginSubRun( const art::SubRun& subRun);

  private:
    int testMode_;
    std::string label_;
    std::vector<std::string> instanceNames_;

  };

  SeedTest01::SeedTest01(fhicl::ParameterSet const& pSet):
    art::EDAnalyzer(pSet),
    testMode_(pSet.get<int>("testMode")),
    label_(),
    instanceNames_(pSet.get<vector<string> >("instanceNames")){

    cout << "Construct SeedTest01. Mode:  "
         << testMode_ <<  " instances:";
    for ( size_t i=0; i< instanceNames_.size(); ++i){
      cout << " " << instanceNames_[i];
    }
    cout << endl;

    art::ServiceHandle<SeedService> seeds;

    if ( testMode_ == 0 ){
      auto seed = seeds->getSeed();
      cout << "Seed is :  " << seed << endl;

      auto seed2 = seeds->getSeed();
      cout << "Seed2 is : " << seed2 << endl;
    }

    if ( testMode_ == 1 ){
      cout << "Here we go ... " << endl;
      for ( size_t i=0; i< instanceNames_.size(); ++i){
        auto seed = seeds->getSeed( instanceNames_[i] );
        cout << "Seed for : "
             << instanceNames_[i] << " is: "
             << seed
             << endl;
      }
      for ( size_t i=0; i< instanceNames_.size(); ++i){
        auto seed = seeds->getSeed( instanceNames_[i] );
        cout << "Seed2 for : "
             << instanceNames_[i] << " is: "
             << seed
             << endl;
      }
    }
  }

  void SeedTest01::analyze(const art::Event& event){
    cerr << "SeedTest01::analyze "
         << event.id()
         << endl;
    if ( testMode_ == 4 ){
      art::ServiceHandle<SeedService> seeds;
      auto seed = seeds->getSeed();
      cout << "Seed from analyze is :  " << seed << endl;
    }

  }

  void SeedTest01::beginRun( const art::Run& run){
    cerr << "SeedTest01::beginRun "
         << run.id()
         << endl;
    if ( testMode_ == 2 ){
      art::ServiceHandle<SeedService> seeds;
      auto seed = seeds->getSeed();
      cout << "Seed from beginRun is :  " << seed << endl;
    }
  }

  void SeedTest01::beginSubRun( const art::SubRun& subRun){
    cerr << "SeedTest01::beginSubRun "
         << subRun.id()
         << endl;
    if ( testMode_ == 3 ){
      art::ServiceHandle<SeedService> seeds;
      auto seed = seeds->getSeed();
      cout << "Seed from beginSubRun is :  " << seed << endl;
    }
  }

} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::SeedTest01);
