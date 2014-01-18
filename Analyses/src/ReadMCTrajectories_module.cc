//
//
// $Id: ReadMCTrajectories_module.cc,v 1.1 2014/01/18 03:25:15 kutschke Exp $
// $Author: kutschke $
// $Date: 2014/01/18 03:25:15 $
//
// Contact person Rob Kutschke
//

// C++ includes.
#include <iostream>
#include <string>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"

// Mu2e includes.
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/MCTrajectoryCollection.hh"

#include "TNtuple.h"

using namespace std;

namespace mu2e {

  class ReadMCTrajectories : public art::EDAnalyzer {
  public:

    explicit ReadMCTrajectories(fhicl::ParameterSet const& pset);

    virtual void beginJob();
    void analyze(const art::Event& e);

  private:

    art::InputTag trajectoriesTag_;

    int maxPrint_;
    int nPrint_;

    TNtuple* ntup_;

  };

  ReadMCTrajectories::ReadMCTrajectories(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    trajectoriesTag_(pset.get<std::string>("trajectoriesTag")),
    maxPrint_(pset.get<int>("maxPrint",20)),
    nPrint_(0),
    ntup_(nullptr){
  }

  void ReadMCTrajectories::beginJob(){

    art::ServiceHandle<art::TFileService> tfs;
    ntup_ = tfs->make<TNtuple>( "ntup", "Hit ntuple", "evt:x:y:z:t");

  }


  void ReadMCTrajectories::analyze(const art::Event& event) {

    bool doPrint = ( ++nPrint_ <= maxPrint_ );

    auto trajectories = event.getValidHandle<MCTrajectoryCollection>(trajectoriesTag_);

    float nt[5];

    for ( auto const& i : *trajectories ){
      MCTrajectory const& traj = i.second;
      if ( doPrint ) {
        cout << "Trajectory: "
             << event.id().event() << " "
             << i.first   << " "
             << traj.sim() << " "
             << traj.size()
             << endl;
      }
      nt[0] = event.id().event();
      for ( auto const& pos : traj.points() ){
        nt[1] = pos.x();
        nt[2] = pos.y();
        nt[3] = pos.z();
        nt[4] = pos.t();
        ntup_->Fill(nt);
      }
    }


  }

}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::ReadMCTrajectories);
