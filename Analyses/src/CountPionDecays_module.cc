//
// Count decays of charged pions to check that the pi e nu code is working properly.
//
// $Id: CountPionDecays_module.cc,v 1.2 2013/10/21 20:44:04 genser Exp $
// $Author: genser $
// $Date: 2013/10/21 20:44:04 $
//
// Original author Rob Kutschke
//

#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/ProcessCode.hh"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

using namespace std;

namespace mu2e {

  class CountPionDecays : public art::EDAnalyzer {
  public:

    explicit CountPionDecays(fhicl::ParameterSet const& pset);
    virtual ~CountPionDecays() { }

    void analyze(const art::Event& e);

    void endJob();

  private:

    // Module label of the g4 module that made the hits.
    std::string _g4ModuleLabel;

    // The following count pi- and pi+ separately.

    // Number of pions seen in all events.
    vector<int>_nPi;

    // Number of pions that decayed.
    vector<int>_nPiDecays;

    // Number of decays to e nu.
    vector<int>_nENu;

    // Number of decays to mu nu
    vector<int>_nMuNu;

    // Number of decays to other modes - should be zero.
    vector<int>_nOther;

  };

  CountPionDecays::CountPionDecays(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    _g4ModuleLabel(pset.get<string>("g4ModuleLabel")),
    _nPi(2,0),
    _nPiDecays(2,0),
    _nENu(2,0),
    _nMuNu(2,0),
    _nOther(2,0){
  }

  // For each event, look at tracker hits and calorimeter hits.
  void CountPionDecays::analyze(const art::Event& event) {

    art::Handle<SimParticleCollection> simsHandle;
    event.getByLabel(_g4ModuleLabel,simsHandle);
    SimParticleCollection const& sims(*simsHandle);

    for ( SimParticleCollection::const_iterator i=sims.begin();
          i != sims.end(); ++i ){
      SimParticle const& sim(i->second);

      // I want to count pi+ and pi- separately.
      int idx(0);
      if ( sim.pdgId() == PDGCode::pi_plus ){
        idx = 0;
      } else if ( sim.pdgId() == PDGCode::pi_minus ){
        idx = 1;
      } else{
        continue;
      }

      ++_nPi[idx];

      if ( sim.stoppingCode() != ProcessCode::Decay ) continue;
      ++_nPiDecays[idx];


      // Count daughters produced by decay; find electrons and muons among these.
      // There may be other daughters, such as delta rays, produced prior to decay.
      std::vector<art::Ptr<SimParticle> > const& daughters = sim.daughters();
      typedef std::vector<art::Ptr<SimParticle> >::const_iterator Iter;

      bool hasE(false);
      bool hasMu(false);
      int  nDecayDaughters(0);
      for ( Iter j=daughters.begin(); j!=daughters.end(); ++j){
        SimParticle const& dau(**j);
        if ( dau.creationCode() == ProcessCode::Decay ){
          ++nDecayDaughters;
          if ( std::abs(dau.pdgId()) == PDGCode::e_minus ){
            hasE  = true;
          } else if ( std::abs(dau.pdgId()) == PDGCode::mu_minus ){
            hasMu = true;
          }
        }
      }

      // Classify the decay as pi -> e nu, pi ->mu nu or other.
      if ( nDecayDaughters == 2 ) {
        if ( hasE && !hasMu ){
          ++_nENu[idx];
        } else if ( hasMu && !hasE ){
          ++_nMuNu[idx];
        } else{
          ++_nOther[idx];
        }
      } else{
        ++_nOther[idx];
      }

    } // end loop over SimParticles

  } // end analyze

  void CountPionDecays::endJob() {

    vector<string> names;
    names.push_back("pi+");
    names.push_back("pi-");

    cout << "\n\nSummary of Pion decays: \n" << endl;
    cout << "     Species"
         << "  Pions"
         << "  Decays"
         << "    ENu"
         << "   MuNu"
         << "  Other"
         << "    Ratio"
         << endl;
    for ( size_t i=0; i<names.size(); ++i ){

      // Ratio
      double r = ( _nPiDecays[i] > 0) ? double(_nENu[i])/double(_nPiDecays[i]) : -1.;
      cout << "    "
           << setw(8) << names[i]      << " "
           << setw(6) << _nPi[i]       << " "
           << setw(7) << _nPiDecays[i] << " "
           << setw(6) << _nENu[i]      << " "
           << setw(6) << _nMuNu[i]     << " "
           << setw(6) << _nOther[i]    << " "
           << setw(8) << r
           << endl;
    }
    cout << endl;

  } // end endJob


}  // end namespace mu2e

using mu2e::CountPionDecays;
DEFINE_ART_MODULE(CountPionDecays);
