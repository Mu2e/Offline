//
// A stripped down version of ReadStrawDigi that only looks at reco info
//
// Original author Rob Kutschke.
//

// Mu2e includes.
#include "RecoDataProducts/inc/StrawDigi.hh"

#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Provenance.h"

#include "fhiclcpp/ParameterSet.h"

#include "TH1F.h"

#include <iostream>

using namespace std;

namespace mu2e {

  class ReadStrawDigiReco : public art::EDAnalyzer {
  public:
    explicit ReadStrawDigiReco(fhicl::ParameterSet const& pset);

    virtual void beginJob() override;
    virtual void analyze( art::Event const& e) override;

  private:

    int _diagLevel;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Input tag describing which hits to look at
    art::InputTag _digisTag;

    // Some diagnostic histograms.
    TH1F* _hStrawId         = nullptr;
    TH1F* _hDigiTime0       = nullptr;
    TH1F* _hDigiTime1       = nullptr;
    TH1F* _hDigiDeltaTime   = nullptr;
    TH1F* _hDigiDeltaTimeDetail  = nullptr;
    TH1F* _hNADC            = nullptr;
    TH1F* _hADC             = nullptr;
    TH1F* _hMeanADC         = nullptr;
    TH1F* _hMaxADC          = nullptr;
    TH1F* _hNDigis          = nullptr;
    TH1F* _hNDigis1         = nullptr;
    TH1F* _hNDigisPerWire   = nullptr;

    int _ncalls = 0;
  };

} // end of namespace mu2e

mu2e::ReadStrawDigiReco::ReadStrawDigiReco(fhicl::ParameterSet const& pset):
  art::EDAnalyzer(pset),
  _diagLevel(pset.get<int>("diagLevel",0)),
  _maxFullPrint(pset.get<int>("maxFullPrint",0)),
  _digisTag(pset.get<std::string>("digisTag")){
}

void mu2e::ReadStrawDigiReco::beginJob(){

  if ( _diagLevel > 0 ) {
    cout << "ReadStrawDigiReco Diaglevel: "
         << _diagLevel << " "
         << _maxFullPrint
         << endl;
  }

  art::ServiceHandle<art::TFileService> tfs;

  _hStrawId         = tfs->make<TH1F>( "StrawId",    "StrawId",                          270,   0.,   54000. );
  _hDigiTime0       = tfs->make<TH1F>( "DigiTime0",      "Digi TDC end 0",               264,   0.,  264000. );
  _hDigiTime1       = tfs->make<TH1F>( "DigiTime1",      "Digi TDC end 1",               264,   0.,  264000. );
  _hDigiDeltaTime   = tfs->make<TH1F>( "DigiDeltaTime",   "Digi Delta Time (ns)",         264, -264000., 264000. );
  _hDigiDeltaTimeDetail = tfs->make<TH1F>( "DigiDeltaTimeDetail",  "Digi Delta Time (ns)",   220, -100., 10000. );
  _hNADC            = tfs->make<TH1F>( "NADC",  "Number of Counts in Waveform",       100,  0., 100. );
  _hADC             = tfs->make<TH1F>( "ADC",   "ADC spectrum",              250,  0., 2500. );
  _hMeanADC         = tfs->make<TH1F>( "MeanADC",  "Mean ADC ",              250,  0., 2500. );
  _hMaxADC          = tfs->make<TH1F>( "MaxADC",   "Max ADC ",               250,  0., 2500. );
  _hNDigis          = tfs->make<TH1F>( "hNDigis",          "Number of straw hits",           500,   0.,   500. );
  _hNDigis1         = tfs->make<TH1F>( "hNDigis1",         "Number of straw hits",           200,   0., 12000. );
  _hNDigisPerWire   = tfs->make<TH1F>( "hNDigisPerWire",   "Number of hits per straw",        10,   0.,    10. );

  }

void mu2e::ReadStrawDigiReco::analyze(art::Event const& evt) {

  ++_ncalls;

  auto const& digis = *evt.getValidHandle<StrawDigiCollection>( _digisTag );
  _hNDigis->Fill(digis.size());
  _hNDigis1->Fill(digis.size());
  auto const& digiadcs = *evt.getValidHandle<StrawDigiADCWaveformCollection>( _digisTag );

  if ( _diagLevel > 1 ) {
    cout << "ReadStrawDigiReco: Total number of straw digis = " << digis.size() << endl;
  }

  // Counter for number of digis on each wire.
  std::map<StrawId,int> nhperwire;

  for (size_t idigi=0;idigi<digis.size();idigi++){
    StrawDigi const& digi = digis.at(idigi);
    StrawDigiADCWaveform const& digiadc = digiadcs.at(idigi);
    StrawId id = digi.strawId();

    _hStrawId->Fill(id.asUint16());

    // Calculate number of digis per wire
    ++nhperwire[id];

    auto t0 = digi.TDC(StrawEnd::cal);
    auto t1 = digi.TDC(StrawEnd::hv);
    auto const& adcs = digiadc.samples();

    _hDigiTime0 ->Fill(t0);
    _hDigiTime1 ->Fill(t1);
    _hDigiDeltaTime->Fill(t1-t0);
    _hDigiDeltaTimeDetail->Fill(t1-t0);
    _hNADC->Fill( adcs.size() );
    int sum{0};
    unsigned short maxadc{0};
    for ( auto adc : adcs ){
      sum += adc;
      maxadc = std::max( maxadc, adc);
      _hADC->Fill(adc);
    }
    int mean = ( adcs.size() != 0 ) ? sum/adcs.size() : -1.;
    _hMeanADC->Fill(mean);
    _hMaxADC->Fill(maxadc);

    if ( (int(evt.id().event()) < _maxFullPrint) && (_diagLevel > 3) ) {
      cout << "ReadStrawDigiReco: "
           << evt.id().event()        << " "
           << digi.strawId()          << " "
           << digi.TDC(StrawEnd::cal) << " "
           << digi.TDC(StrawEnd::hv)  << " "
           << digi.TOT(StrawEnd::cal) << " "
           << digi.TOT(StrawEnd::hv)  << " "
           << sum << " "
           << maxadc
           << endl;
    }

  } // end loop over digis.

  //for( std::map<StrawId,int>::iterator it=nhperwire.begin(); it!= nhperwire.end(); ++it ) {
  for ( auto const& i : nhperwire){
    _hNDigisPerWire->Fill(i.second);
  }

} // end of ::analyze.

DEFINE_ART_MODULE(mu2e::ReadStrawDigiReco);
