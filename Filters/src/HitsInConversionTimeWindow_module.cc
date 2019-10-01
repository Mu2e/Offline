//
// This is not yet general purpose code; it still contains some elements that are
// specific to the job of selecting subsets of events for use in the event displays
// shown at the December 2011 PAC meeting.
//
// Find the conversion electron and all StrawHits to which it contributes.
// Select all StrawHits in a time window around the conversion electron; add these
// and their MC truth as new data products in the event.  Create new SimParticle and
// PointTrajectory collections that they contain only those tracks that either create
// one of the selected hits or are an ancestor of such a track.  Put these new
// collections into the event.
//
// For mixed events, there is only one new StrawHitCollection but there are many new
// SimParticleCollections and PointTrajectoryCollections.
//
// We anticipate a future development in which this module will look for both StrawHits
// and CaloHits that are in time with the conversion electron, and will keep all
// particles in the ancestry of any such hit.
//
// $Id: HitsInConversionTimeWindow_module.cc,v 1.3 2013/03/15 15:52:04 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/15 15:52:04 $
//
// Contact person Rob Kutschke.
//
// To do:
// 1) Check all G4Status object, not just the g4run one.
// 2) Get instance names from the pset.  Check that instance names match the branch names.
//    Not 100% sure what match means.
// 3) Remove unneeded headers.

// Mu2e includes.
#include "GeneralUtilities/inc/SequenceStatistics.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/StatusG4.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/PointTrajectoryCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "Mu2eUtilities/inc/compressSimParticleCollection.hh"
#include "Mu2eUtilities/inc/compressPointTrajectoryCollection.hh"
#include "Mu2eUtilities/inc/checkSimParticleCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"

// art includes.
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"

// Root includes
#include "TH1F.h"
#include "TNtuple.h"

// Other includes
#include "CLHEP/Units/SystemOfUnits.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ includes
#include <iostream>
#include <set>

using namespace std;

namespace mu2e {

  // A helper struct to hold information about one SimParticle.
  struct Info {
    Info():
      simId(),
      hitCount(0),
      hasTrajectory(false),
      onPathToHit(false),
      isPrimary(false){}

    Info( cet::map_vector_key id, bool primary ):
      simId(id),
      hitCount(0),
      hasTrajectory(false),
      onPathToHit(false),
      isPrimary(primary){}

    cet::map_vector_key simId;   // The id of this SimParticle.
    int hitCount;                // Does this particle have hits?
    bool hasTrajectory;          // Is there an associated PointTrajectory object?
    bool onPathToHit;            // Is this particle in the ancestry of a hit?
    bool isPrimary;              // Is this particle a primary particle?

  };

  // A collection of Info objects that is indexed in parallel to a SimParticleCollection object.
  typedef cet::map_vector<Info> SimInfoCollection;

  class HitsInConversionTimeWindow : public art::EDFilter {
  public:
    explicit HitsInConversionTimeWindow(fhicl::ParameterSet const& pset);
    virtual ~HitsInConversionTimeWindow() { }

    virtual bool filter( art::Event& event);

    virtual bool beginRun(art::Run &run);
    virtual void endJob();

  private:

    // A map type to convert a product ID into an index into the Handle collections.
    typedef map<art::ProductID,int> IdMap;

    // Cut to define the time window around the selected conversion electron.
    double timeWindow_;

    // Module label of the g4 module that made the generated particles
    std::string generatorModuleLabel_;

    // Module label of the module that made the StepPointMCCollections.
    std::string g4ModuleLabel_;

    // Module labels of the modules that made the collection of reco hits.
    std::string strawHitMakerLabel_;
    std::string crystalHitMakerLabel_;

    std::vector<std::string> instanceNames_;

    // Instance names of the StepPointMCCollections.
    std::string trackerStepPoints_;
    std::string caloStepPoints_;
    std::string caloROStepPoints_;
    std::string foilStepPoints_;
    std::string crvStepPoints_;
    std::string vDetStepPoints_;

    // Histogram pointers.
    TH1F* hNhits_;
    TH1F* hNhitsAll_;
    TH1F* hNhitsIn_;
    TH1F* hTime_;
    TH1F* hTimeAll_;
    TH1F* hEnergy_;
    TH1F* hEnergyAll_;
    TH1F* hEnergyIn_;
    TH1F* hDeltaT_;
    TH1F* htRMS_;
    TH1F* hTSpread_;
    TH1F* hTSpreadAll_;

    // Number of events that pass the filter.
    int nPassed_;

    // Some utility functions that break the complete job into smaller pieces.

    GenParticle const& findConversionGenParticle( GenParticleCollection const& gens);
    SimParticle const& findConversionSimParticle( SimParticleCollection const& sims, GenParticle const& conversion);

    void findStrawHits ( SimParticle const& sim,
                         StrawHitCollection const& strawHits,
                         PtrStepPointMCVectorCollection const& mcptrs,
                         std::vector<StrawHit const*>& conversionStrawHits,
                         SequenceStatistics& stats );

    void fillIdMap(  std::vector<art::Handle<SimParticleCollection> >const&  simHandles,
                     IdMap& idmap );

    void checkNames(  std::vector<art::Handle<SimParticleCollection> >const&  simHandles,
                      std::vector<art::Handle<PointTrajectoryCollection> >const& trajHandles );

    void doThreeJobs ( StrawHitCollection const&             strawHits,
                       StrawHitMCTruthCollection const&      strawHitTruths,
                       PtrStepPointMCVectorCollection const& mcptrs,
                       double                                t0,
                       art::Event&                           event,
                       set<art::Ptr<SimParticle> >&          contributingSims );

    void fillSimInfo( std::vector<art::Handle<SimParticleCollection> >     const& simHandles,
                      std::vector<art::Handle<PointTrajectoryCollection> > const& trajHandles,
                      set<art::Ptr<SimParticle> >                          const& inTimeSims,
                      vector<SimInfoCollection >                                & counts
                      );


  };

  HitsInConversionTimeWindow::HitsInConversionTimeWindow(fhicl::ParameterSet const& pset):
    art::EDFilter{pset},
    timeWindow_(pset.get<double>("timeWindow",50)),
    generatorModuleLabel_(pset.get<string>("generatorModuleLabel")),
    g4ModuleLabel_(pset.get<string>("g4ModuleLabel")),
    strawHitMakerLabel_(pset.get<string>("strawHitMakerLabel")),
    instanceNames_(),
    hNhits_(0),
    hNhitsAll_(0),
    hNhitsIn_(0),
    hTime_(0),
    hTimeAll_(0),
    hEnergy_(0),
    hEnergyAll_(0),
    hEnergyIn_(0),
    hDeltaT_(0),
    htRMS_(0),
    hTSpread_(0),
    hTSpreadAll_(0),
    nPassed_(0){

    instanceNames_.push_back("dioMixer");
    instanceNames_.push_back("neutronMixer");
    instanceNames_.push_back("photonMixer");
    instanceNames_.push_back("protonMixer");
    instanceNames_.push_back("conversion");

    produces<StrawHitCollection>();
    produces<StrawHitMCTruthCollection>();
    produces<PtrStepPointMCVectorCollection>("StrawHitMCPtr");

    for ( size_t i=0; i<instanceNames_.size(); ++i ){
      produces<SimParticleCollection>(instanceNames_.at(i));
      produces<PointTrajectoryCollection>(instanceNames_.at(i));
    }

  }

  // An adapter class that allows compressSimParticleCollection to access the std::map of Info objects.
  class CompressSimAdapter{
  public:
    CompressSimAdapter ( SimInfoCollection const& m ):
      m_(m){}

    bool operator[]( cet::map_vector_key key ) const {
      Info const& info(m_[key]);
      bool keepsim = info.hitCount > 0 || info.onPathToHit;
      return keepsim;
    }

  private:
    SimInfoCollection const& m_;
  };

  // An adapter class that allows compressPointTrajectoryCollection to access the std::map of Info objects.
  class CompressTrajectoryAdapter{
  public:
    CompressTrajectoryAdapter ( SimInfoCollection const& m ):
      m_(m){}

    bool operator[]( cet::map_vector_key key ) const {
      Info const& info(m_[key]);
      bool keep = info.hasTrajectory && ( info.hitCount > 0 || info.onPathToHit );
      return keep;
    }

  private:
    SimInfoCollection const& m_;
  };

  bool HitsInConversionTimeWindow::beginRun(art::Run& ){

    art::ServiceHandle<art::TFileService> tfs;

    hNhits_      = tfs->make<TH1F>( "hNhits",     "Number of Straw Hits, Conversion",   200,    0.,   200.  );
    hNhitsAll_   = tfs->make<TH1F>( "hNhitsAll",  "Number of Straw Hits, All",          200,    0.,  5000.  );
    hNhitsIn_    = tfs->make<TH1F>( "hNhitsIn",   "Number of Straw Hits, In window",    200,    0.,   500.  );
    hTime_       = tfs->make<TH1F>( "hTime",      "Hit time, Conversion;(ns)",          200,    0.,  2000.  );
    hTimeAll_    = tfs->make<TH1F>( "hTimeAll",   "Hit time, All;(ns)",                 200,    0.,  2000.  );
    hEnergy_     = tfs->make<TH1F>( "hEnergy",    "Hit EnergyDeposit Conversion;[keV]", 100,    0.,    10.  );
    hEnergyAll_  = tfs->make<TH1F>( "hEnergyAll", "Hit EnergyDeposit All;[keV]",        100,    0.,    10.  );
    hEnergyIn_   = tfs->make<TH1F>( "hEnergyIn",  "Hit EnergyDeposit All;[keV]",        100,    0.,    10.  );
    hDeltaT_     = tfs->make<TH1F>( "hDeltaT",    "Hit Delta time, Conversion;(ns)",    200,    0.,   200.  );
    htRMS_       = tfs->make<TH1F>( "htRMS",      "rms time spread, Conversion;(ns)",   200,    0.,   200.  );
    hTSpread_    = tfs->make<TH1F>( "hTSpread",   "Time-tmean, Conversion;(ns)",        200, -100.,   100.  );
    hTSpreadAll_ = tfs->make<TH1F>( "hTSpreadAll","Time-tmean, Conversion;(ns)",        200, -100.,   100.  );

    return true;
  }

  bool HitsInConversionTimeWindow::filter(art::Event& event) {

    art::Handle<StatusG4> g4StatusHandle;
    event.getByLabel( g4ModuleLabel_, g4StatusHandle);
    StatusG4 const& g4Status = *g4StatusHandle;

    // Accept only events with good status from G4.
    // Expand this to check all of the mixing status objects or the overall status object.
    if ( g4Status.status() > 1 ) {
      return false;
    }

    art::Handle<GenParticleCollection> gensHandle;
    event.getByLabel( generatorModuleLabel_, gensHandle);
    GenParticleCollection const& gens(*gensHandle);

    art::Handle<SimParticleCollection> simsHandle;
    event.getByLabel(g4ModuleLabel_,simsHandle);
    SimParticleCollection const& sims(*simsHandle);

    art::Handle<StrawHitCollection> strawHitsHandle;
    event.getByLabel(strawHitMakerLabel_, strawHitsHandle);
    StrawHitCollection const& strawHits(*strawHitsHandle);

    art::Handle<StrawHitMCTruthCollection> strawHitTruthsHandle;
    event.getByLabel(strawHitMakerLabel_, strawHitTruthsHandle);
    StrawHitMCTruthCollection const& strawHitTruths(*strawHitTruthsHandle);

    art::Handle<PtrStepPointMCVectorCollection> mcptrsHandle;
    event.getByLabel(strawHitMakerLabel_, "StrawHitMCPtr", mcptrsHandle);
    PtrStepPointMCVectorCollection const& mcptrs(*mcptrsHandle);

    // Find the GenParticle and the SimParticle corresponding to the conversion electron.
    // Expect exactly one per event; if this is no longer true, we must update the find functions.
    GenParticle const& gen(findConversionGenParticle(gens));
    SimParticle const& sim(findConversionSimParticle(sims,gen));

    // Identify the subset of StrawHits that have some energy contribution from the conversion electron.
    // Collect some statistics about them.
    std::vector<StrawHit const*> conversionStrawHits;
    conversionStrawHits.reserve(200);
    SequenceStatistics stats;
    findStrawHits ( sim, strawHits, mcptrs, conversionStrawHits, stats);

    // The set of all SimParticles that contribute to any in-time hit.
    // In-time means within +/- tWindow_ of stat.mean().
    set<art::Ptr<SimParticle> > inTimeSims;

    // Do three jobs at once to avoid repeating a double loop.
    //  - Create the three output data products, newStrawHits, newTruthHits, newMCptrs,
    //    fill them and add them to the event.
    //  - Fill the inTimeSims variable
    //  - Fill some diagnostics histograms.
    doThreeJobs ( strawHits,
                  strawHitTruths,
                  mcptrs,
                  stats.moments().mean(),
                  event,
                  inTimeSims );
    bool keep = !inTimeSims.empty();

    // Get the SimParticle and PointTrajectory Collections from the event.
    art::ProductInstanceNameSelector selector("");
    std::vector<art::Handle<SimParticleCollection> >      simHandles;
    std::vector<art::Handle<PointTrajectoryCollection> > trajHandles;
    event.getMany( selector, simHandles);
    event.getMany( selector, trajHandles);

    // Check that the SimParticle and PointTrajectory collections are properly matched up.
    checkNames( simHandles, trajHandles);

    // Fill information about each SimParticle.
    vector<SimInfoCollection > simsInfos;
    fillSimInfo( simHandles, trajHandles, inTimeSims, simsInfos );

    // Form the compressed SimParticle and PointTrajectory collections; write them out.
    for ( size_t i=0; i<simHandles.size(); ++i ){

      // These are needed to reseat the art::Ptr's inside the SimParticleCollections.
      art::ProductID simsProductId(event.getProductID<SimParticleCollection>(instanceNames_.at(i)));
      art::EDProductGetter const * productGetter = event.productGetter(simsProductId);

      // Adapters to translate from the Info map to the format needed by the compress functions.
      CompressSimAdapter        simAdapter(simsInfos.at(i));
      CompressTrajectoryAdapter trajAdapter(simsInfos.at(i));

      // Make a new data products: a compressed list of SimParticles and one of PointTrajectories.
      unique_ptr<SimParticleCollection> newsimTest( new SimParticleCollection );
      unique_ptr<PointTrajectoryCollection> newtraj( new PointTrajectoryCollection );

      compressSimParticleCollection     ( simsProductId, productGetter, *simHandles.at(i),  simAdapter,  *newsimTest);
      compressPointTrajectoryCollection ( simsProductId, productGetter, *trajHandles.at(i), trajAdapter, *newtraj   );

      // Check self-consistency.
      checkSimParticleCollection(*newsimTest,true);

      event.put(std::move(newsimTest),  instanceNames_.at(i) );
      event.put(std::move(newtraj), instanceNames_.at(i) );

    }

    if ( keep ) ++nPassed_;

    return keep;

  } // end of ::analyze.

  void HitsInConversionTimeWindow::endJob() {
    mf::LogInfo("Summary")
      << "HitsInConversionTimeWindow_module: Number of events passing the filter: "
      << nPassed_
      << "\n";
  }

  // Find the genParticle corresponding to the conversion electron.
  // This is a fast and dirty algorithm that works for the files we have now; it assumes
  // that the generated event is a single conversion electron. Need to make it more general.
  GenParticle const& HitsInConversionTimeWindow::findConversionGenParticle( GenParticleCollection const& gens){

    if ( gens.size() != 1 ){
      throw cet::exception("RANGE")
        << "HitsInConversionTimeWindow_module::findConversionGenParticle: unexpected size of GenParticle Collection: "
        << gens.size()
        << "\n";
    }
    GenParticle const& gen(gens.at(0));

    if ( !gen.generatorId().isConversion() ){
      throw cet::exception("RANGE")
        << "HitsInConversionTimeWindow_module::findConversionGenParticle: unexpected generator type: "
        << gen.generatorId()
        << "\n";
    }

    return gen;
  }  // end findConversionGenParticle

  // Find the SimParticle corresponding to the conversion electron.
  // This is a fast and dirty algorithm that works for the files we have now; it assumes
  // that the generated event is a single conversion electron. Need to make it more general.
  SimParticle const& HitsInConversionTimeWindow::findConversionSimParticle( SimParticleCollection const& sims, GenParticle const& gen){

    if ( sims.empty() ){
      throw cet::exception("RANGE")
        << "HitsInConversionTimeWindow_module::findConversionSimParticle: SimParticleCollection is empty: "
        << sims.size()
        << "\n";
    }
    SimParticle const& sim(sims[cet::map_vector_key(1)]);

    if ( sim.genParticle().get() != &gen ){
      throw cet::exception("RANGE")
        << "HitsInConversionTimeWindow_module::findConversionSimParticle: SimParticle Ptr to gen appears to be wrong: "
        << sim.genParticle().id() << " "
        << sim.genParticle().key()
        << "\n";
    }

    return sim;
  } // end findConversionSimParticle

  // Find all straw hits that contain one or more contributions from the SimParticle passed in as the first argument.
  // Also fill some histogrms and accumulate some statistics about the hits: tmin, tmax, mean, rms.
  void HitsInConversionTimeWindow::findStrawHits ( SimParticle const&                    sim,
                                                   StrawHitCollection const&             strawHitsIn,
                                                   PtrStepPointMCVectorCollection const& mcptrs,
                                                   std::vector<StrawHit const*>&         strawHitsOut,
                                                   SequenceStatistics&                   stats ){

    for ( size_t i=0; i< mcptrs.size(); ++i ){

      PtrStepPointMCVector const& mcptr(mcptrs.at(i));
      StrawHit const&               hit(strawHitsIn.at(i));

      for ( size_t j=0; j<mcptr.size(); ++j){
        StepPointMC const& step(*mcptr.at(j));
        if ( step.simParticle().get() == &sim ){
          double t=hit.time();
          double e=hit.energyDep();
          stats.accumulate(t);
          hTime_->Fill(t);
          hEnergy_->Fill(e/CLHEP::keV);
          strawHitsOut.push_back(&hit);
          break;
        }
      }
    }

    // Some diagnostic histograms.
    hNhits_    ->Fill( stats.n() );
    hNhitsAll_ ->Fill( strawHitsIn.size() );
    htRMS_     ->Fill( stats.moments().rms());
    hDeltaT_   ->Fill( stats.limits().delta() );
    for ( size_t i=0; i<strawHitsOut.size(); ++i){
      hTSpread_->Fill(strawHitsOut.at(i)->time()-stats.moments().mean());
    }


  } // end findStrawHits

  // Fill the idmap variable.
  void HitsInConversionTimeWindow::fillIdMap(  std::vector<art::Handle<SimParticleCollection> >const&  simHandles,
                               IdMap& idmap ){

    for ( size_t i=0; i<simHandles.size(); ++i ){
      art::Provenance const& simProv( *simHandles.at(i).provenance());
      idmap[simProv.productID()] = i;
    }
  }

  // Check that the SimParticleCollections and the PointTrajectoryCollections are properly paired up.
  void HitsInConversionTimeWindow::checkNames(  std::vector<art::Handle<SimParticleCollection> >const&  simHandles,
                                std::vector<art::Handle<PointTrajectoryCollection> >const& trajHandles ){

    if ( simHandles.size() !=  trajHandles.size() ){
        throw cet::exception("DATA")
          << "HitsInConversionTimeWindow_module::checkNames: SimPartile and PointTrajectory Data products have different numbers: \n"
          << "Sims:         " << simHandles.size()  << "\n"
          << "Trajectories: " << trajHandles.size() << "\n";
    }

    for ( size_t i=0; i<simHandles.size(); ++i ){
      art::Provenance const& simProv( *simHandles.at(i).provenance());
      art::Provenance const& trajProv( *trajHandles.at(i).provenance());

      // Strip of the leading element of the branch names ( the friendly name of the type ),
      // leaving 3 elements that should match.
      size_t l1(simProv.branchName().find("_"));
      size_t l2(trajProv.branchName().find("_"));

      if ( simProv.branchName().substr(l1) != trajProv.branchName().substr(l2) ){
        throw cet::exception("DATA")
          << "HitsInConversionTimeWindow_module::fillIdmapAndCheck: SimPartile and PointTrajectory Data products are not synchronized: \n"
          << "Sims:         " << simProv  << "\n"
          << "Trajectories: " << trajProv << "\n";
      }

    }

  } // end checkNames


  // Do three jobs at once to avoid repeating a double loop.
  //  - Create the three output data products, newStrawHits, newTruthHits, newMCptrs,
  //    fill them and add them to the event.
  //  - Fill the inTimeSims variable
  //  - Fill some diagnostics histograms.
  void HitsInConversionTimeWindow::doThreeJobs ( StrawHitCollection const&             strawHits,
                                 StrawHitMCTruthCollection const&      strawHitTruths,
                                 PtrStepPointMCVectorCollection const& mcptrs,
                                 double                                t0,
                                 art::Event&                           event,
                                 set<art::Ptr<SimParticle> >&          inTimeSims ){

    // A guess at the reserved size of the output StrawHit info.
    static const int initialSize = 500;

    // Create output collections for the StrawHits and their MC truth.
    unique_ptr<StrawHitCollection>             newStrawHits(new StrawHitCollection);
    unique_ptr<StrawHitMCTruthCollection>      newTruthHits(new StrawHitMCTruthCollection);
    unique_ptr<PtrStepPointMCVectorCollection> newMCptrs(new PtrStepPointMCVectorCollection);
    newStrawHits->reserve(initialSize);
    newTruthHits->reserve(initialSize);
    newMCptrs->reserve(initialSize);

    for ( size_t i=0; i<strawHits.size(); ++i ){

      StrawHit const&             hit(strawHits.at(i));
      PtrStepPointMCVector const& mcptr(mcptrs.at(i));

      // Fill histograms
      hTimeAll_->Fill(hit.time());
      hEnergyAll_->Fill(hit.energyDep()/CLHEP::keV);
      hTSpreadAll_->Fill(hit.time()-t0);

      // Select hits that are in time with the conversion electron.
      if ( std::abs(hit.time()-t0) < timeWindow_ ){

        // Fill histograms
        hEnergyIn_->Fill(hit.energyDep()/CLHEP::keV);

        // Fill output data products.
        newStrawHits->push_back(hit);
        newTruthHits->push_back( strawHitTruths.at(i));
        newMCptrs->push_back( mcptr);

        // Append to the set of SimParticles that contribute to any in-time hit.
        for ( size_t j=0; j<mcptr.size(); ++j){
          StepPointMC const& step(*mcptr.at(j));
          inTimeSims.insert(step.simParticle());
        }
      }
    }

    // A final histogram.
    hNhitsIn_->Fill(newStrawHits->size());

    // Add data products to the event.
    event.put(std::move(newStrawHits));
    event.put(std::move(newTruthHits));
    event.put(std::move(newMCptrs),"StrawHitMCPtr");


  } // end doThreeJobs

  // Fill-in the simsInfo object, that holds information about each SimParticle in the event.
  void HitsInConversionTimeWindow::fillSimInfo( std::vector<art::Handle<SimParticleCollection> >     const& simHandles,
                                                std::vector<art::Handle<PointTrajectoryCollection> > const& trajHandles,
                                                set<art::Ptr<SimParticle> >                          const& inTimeSims,
                                                vector<SimInfoCollection >                                & simsInfos
                                                ){

    // A map from ProductID into to the position within the sumHandles and trajHandles.
    IdMap idmap;
    fillIdMap( simHandles, idmap);

    for ( size_t i=0; i<simHandles.size(); ++i ){

      // One pair of (SimParticleCollection,PointTrajectoryColllection).
      SimParticleCollection     const&  sims( *simHandles.at(i));
      PointTrajectoryCollection const& trajs(*trajHandles.at(i));

      // Each SimInfoCollection corresponds to one (SimParticleCollection,PointTrajectoryColllection) pair.
      simsInfos.push_back( SimInfoCollection() );
      SimInfoCollection& infomap = simsInfos.back();
      infomap.reserve(simHandles.at(i)->size());

      // Initialize the SimInfoCollection with basic information about all SimParticles.
      for ( SimParticleCollection::const_iterator j=sims.begin(), je=sims.end();
            j !=je; ++j ){
        infomap[j->first] = Info(j->first,j->second.isPrimary());
      }

      // Record which SimParticles have PointTrajectories.
      for ( PointTrajectoryCollection::const_iterator j=trajs.begin(),je=trajs.end();
            j != je; ++j){
        cet::map_vector<Info >::iterator k = infomap.find(j->first);
        if ( k == infomap.end() ){
          throw cet::exception("DATA")
            << "HitsInConversionTimeWindow_module::fillSimInfo: PointTrajectory does not have a matching infomap entry: \n"
            << "Collection number: " << i << " "
            << trajHandles.at(i).provenance()->branchName() << "\n"
            << "Particle number: " << j->first << " " << j->second.simId() << "\n";
        }
        k->second.hasTrajectory = true;
      }
    } // end loop over simHandles


    // Update simsInfo with information about which particles made hits.
    for ( set<art::Ptr<SimParticle> >::const_iterator i=inTimeSims.begin(), e=inTimeSims.end();
          i != e; ++i){

      // Which SimParticleCollection is this SimParticle from?
      cet::map_vector_key key(i->key());
      int idx = idmap[i->id()];

      // Find the corresponding SimsInfoCollection.
      SimInfoCollection& infomap = simsInfos.at(idx);

      // Record that this SimParticle made a hit.
      SimInfoCollection::iterator k = infomap.find(key);
      if ( k == infomap.end() ){
        throw cet::exception("DATA")
          << "HitsInConversionTimeWindow_module::fillSimInfo: inTimeSims does not have a matching infomap entry: \n"
          << "Idx: " << idx << "  Particle id: " << key << "\n";
      }
      Info& info(k->second);
      ++info.hitCount;
      info.onPathToHit = true;

      // Mark all ancestors of this SimParticle.
      SimParticle const* particle = i->get();
      while ( particle->isSecondary() ){
        cet::map_vector_key parentKey = cet::map_vector_key(particle->parent().key());
        SimInfoCollection::iterator k = infomap.find(parentKey);
        if ( k == infomap.end() ){
          throw cet::exception("DATA")
            << "HitsInConversionTimeWindow_module::fillSimInfo: cannot find parent: \n"
            << "Idx: " << idx << "  Particle id: " << particle->id()
            << " Parent id: " << parentKey
            << "\n";
        }
        Info& parentInfo(k->second);
        parentInfo.onPathToHit = true;

        particle = particle->parent().get();
      }

    }
  } // end HitsInConversionTimeWindow::fillSimInfo

} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::HitsInConversionTimeWindow);
