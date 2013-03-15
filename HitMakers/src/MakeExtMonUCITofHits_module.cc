//
// An EDProducer Module that reads ExtMonUCITofHit Stepping MC objects and turns them into
// ExtMonUCITofHit objects, collection
//
// $Id: MakeExtMonUCITofHits_module.cc,v 1.9 2013/03/15 15:52:04 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/15 15:52:04 $
//  
//  

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>
#include <map>
#include <vector>
#include <utility>

// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "ExtinctionMonitorUCIGeom/inc/ExtMonUCI.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "RecoDataProducts/inc/ExtMonUCITofHitCollection.hh"
#include "MCDataProducts/inc/ExtMonUCITofHitMCTruthCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"

// Other includes.
#include "CLHEP/Vector/ThreeVector.h"

using namespace std;

namespace mu2e {

  // Utility class (structure) to hold Tof Hit info for G4 hits

  class TOFHit {
  public:

    int    _hit_id;
    double _edep;
    double _edep_corr;
    int    _charged;
    double _time;

    TOFHit(int hit_id, double edep, double edep1, int charged, double time):
      _hit_id(hit_id), _edep(edep), _edep_corr(edep1),
      _charged(charged), _time(time) { }

    // This operator is overloaded in order to time-sort the hits
    bool operator <(const TOFHit& b) const { return (_time < b._time); }

  };



  //--------------------------------------------------------------------
  //
  //
  class MakeExtMonUCITofHits : public art::EDProducer {
  public:
    explicit MakeExtMonUCITofHits(fhicl::ParameterSet const& pset) :

      // Parameters
      _diagLevel(pset.get<int>("diagLevel",0)),
      _maxFullPrint(pset.get<int>("maxFullPrint",5)),
      _stepPoints(pset.get<string>("extmonUCITofStepPoints","ExtMonUCITof")),
      _g4ModuleLabel(pset.get<string>("g4ModuleLabel")),
      _extmonUCITofModuleLabel(pset.get<std::string>("extmonUCITofModuleLabel", "ExtMonUCITofHitsMaker")),
      _messageCategory("ExtMonUCITofHitsMaker"){

      // Tell the framework what we make.
      produces<ExtMonUCITofHitCollection>();
      produces<ExtMonUCITofHitMCTruthCollection>();
      produces<PtrStepPointMCVectorCollection>("ExtMonUCITofHitsMCPtr");

    }
    virtual ~MakeExtMonUCITofHits() { }

    virtual void beginJob();

    void produce( art::Event& e);

  private:

    CLHEP::Hep3Vector _mu2eOrigin;
    SimpleConfig const * _config;

    // Diagnostics level.
    int _diagLevel;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Name of the StepPoint collection
    std::string _stepPoints;

    // Parameters
    string _g4ModuleLabel;  // Name of the module that made these hits.

    string _extmonUCITofModuleLabel; // Name of the module that made the calo hits.

    // A category for the error logger.
    const std::string _messageCategory;

    void makeExtMonUCITofHits(const art::Handle<StepPointMCCollection>&,
                              ExtMonUCITofHitCollection &,
                              ExtMonUCITofHitMCTruthCollection&,
                              PtrStepPointMCVectorCollection& );

  };

  void MakeExtMonUCITofHits::beginJob(){
  }

  void
  MakeExtMonUCITofHits::produce(art::Event& event) {

    if ( _diagLevel > 2 ) cout << "MakeExtMonUCITofHits: produce() begin event " << event.id() << endl;

    GeomHandle<WorldG4> worldGeom;
    _mu2eOrigin = worldGeom->mu2eOriginInWorld();

    if (!_config)
    { 
      // Get access to the master geometry system and its run time config.
      art::ServiceHandle<GeometryService> geom;
      _config = &(geom->config());
    }

    if ( _diagLevel > 2 )
    {
      std::cout << "Mu2e Origin : " << _mu2eOrigin << std::endl;

    }

    static int ncalls(0);
    ++ncalls;

    // Check that extmon_uci geometry description exists
    art::ServiceHandle<GeometryService> geom;
    if( ! geom->hasElement<ExtMonUCI::ExtMon>()) return;

    // A container to hold the output hits.
    unique_ptr<ExtMonUCITofHitCollection>        tofHits   (new ExtMonUCITofHitCollection);
    unique_ptr<ExtMonUCITofHitMCTruthCollection> tofMCHits (new ExtMonUCITofHitMCTruthCollection);
    unique_ptr<PtrStepPointMCVectorCollection>   tofMCptrHits (new PtrStepPointMCVectorCollection);

    // Ask the event to give us a handle to the requested hits.
    art::Handle<StepPointMCCollection> points;
    event.getByLabel(_g4ModuleLabel,_stepPoints,points);
    int nHits = points->size();

    if ( _diagLevel > 2 )
    {
      cout << "nHits " << nHits << endl;
    }

    if( nHits>0 ) {
      makeExtMonUCITofHits(points, *tofHits, *tofMCHits, *tofMCptrHits);
    }

    if ( ncalls < _maxFullPrint && _diagLevel > 2 ) {
      cout << "MakeExtMonUCITofHits: Total number of tof hits = "
           << tofHits->size()
           << endl;
      cout << "MakeExtMonUCITofHits: Total number of tof MC hits = "
           << tofMCHits->size()
           << endl;
    }

    // Add the output hit collection to the event
    event.put(std::move(tofHits));
    event.put(std::move(tofMCHits));
    event.put(std::move(tofMCptrHits),"ExtMonUCITofHitsMCPtr");

    if ( _diagLevel > 2 ) cout << "MakeExtMonUCITofHits: produce() end" << endl;

  } // end of ::analyze.

  void MakeExtMonUCITofHits::makeExtMonUCITofHits (const art::Handle<StepPointMCCollection>& steps,
                                      ExtMonUCITofHitCollection& tofHits,
                                      ExtMonUCITofHitMCTruthCollection& tofHitsMCTruth,
                                      PtrStepPointMCVectorCollection& tofHitsMCptr ) {
    // Get extmon_uci geometry description
    ExtMonUCI::ExtMon const & extmon_uci = *(GeomHandle<ExtMonUCI::ExtMon>());
    int nTofStations = extmon_uci.nTofStations();
    int nTofSegments = extmon_uci.nTofSegments();
    if ( _diagLevel >= 3)
    {
      std::cout << "MakeExtMonUCITofHits::makeExtMonUCITofHits : " << " nTofStations " << nTofStations
                << " nTofSegments " << nTofSegments << std::endl;
    }

    bool keepNoEdep = _config->getBool("extmon_uci.keepNoEdep", false);
    if ( _diagLevel >= 2)
    {
      std::cout << "MakeExtMonUCITofHits::makeExtMonUCITofHits : " << " keepNoEdep " << keepNoEdep << std::endl;
    }

    double timeGap = _config->getDouble("extmon_uci.tofTimeGap" );
    if ( _diagLevel >= 2)
    {
      std::cout << "MakeExtMonUCITofHits::makeExtMonUCITofHits : " << " timeGap " << timeGap << " ns " << std::endl;
    }

    // Organize steps by tof segments

    int nSteps   = steps->size();

    if ( _diagLevel >= 3) 
    {
      std::cout << "MakeExtMonUCITofHits::makeExtMonUCITofHits : " << " nSteps " << nSteps << std::endl;
    }

    // First vector is list of tof steps
    typedef std::map<int, std::vector<int> > HitMap;
    HitMap hitmap;

    for ( int i=0; i<nSteps; ++i){
      StepPointMC const& h = (*steps)[i];
      vector<int> &steps_id = hitmap[h.volumeId()];
      steps_id.push_back(i);
      if ( _diagLevel >= 3) 
      {
        std::cout << "MakeExtMonUCITofHits::makeExtMonUCITofHits : Tof Hit volume Id " << h.volumeId()
                  << " t " << h.time() << " e " << h.eDep() << std::endl;
      }
    }

    // Loop over all readout elements to form extmon_uci tof hits
    vector<TOFHit> tof_hits;

    for(HitMap::const_iterator aHit = hitmap.begin(); aHit != hitmap.end(); ++aHit ) {

      // Tof Segment Id
      int tofId = aHit->first;

      // Prepare info for hit creation before the next iteration
      tof_hits.clear();

      // Loop over all hits found for this Tof Segment
      vector<int> const& isteps = aHit->second;
      if ( _diagLevel >= 3)
      {
        std::cout << "isteps size " << isteps.size() << std::endl;
      } 
 
      // Loop over steps inside the segment
      for( size_t i=0; i<isteps.size(); i++ ) {

        int hitRef = isteps[i];
        StepPointMC const& h = (*steps)[hitRef];

        double edep    = h.eDep(); // each tof segment has its  energy deposit assigned to it
        if ( _diagLevel >= 3)
        {
          std::cout << "time " << h.time() << "  edep " << edep << std::endl;
        }
        if( edep<0.0 || (edep==0.0 && !keepNoEdep) ) continue; // Do not create hit if there is no energy deposition

        /*
          // Hit position in Mu2e frame
          //h_pos = h.position();

          // Hit position in local crystal frame
          // CLHEP::Hep3Vector posLocal = cal.toCrystalFrame(roid,pos);
    
          // Calculate correction for edep
          // double edep_corr = edep * (1.0+(posLocal.z()/length)*nonUniform/2.0);
        */

        tof_hits.push_back(TOFHit(hitRef, edep, 0.0, 0, h.time()));

      }  // loop of i<isteps.size()

      if (tof_hits.size() == 0) continue;

      // Sort hits by time
      sort(tof_hits.begin(), tof_hits.end());

      int stationId = tofId / nTofSegments; 
      int segmentId = tofId % nTofSegments;    

      double h_time    = tof_hits[0]._time;
      double h_edep    = tof_hits[0]._edep;
      double h_edepc   = tof_hits[0]._edep_corr;
      //int    h_charged = tof_hits[0]._charged;
      PtrStepPointMCVector mcptr_tof;
 
      mcptr_tof.push_back( art::Ptr<StepPointMC>( steps, tof_hits[0]._hit_id));

      for( size_t i=1; i<tof_hits.size(); ++i ) {
        if( (tof_hits[i]._time-tof_hits[i-1]._time) > timeGap ) {
          // Save current hit
          tofHits.push_back(       ExtMonUCITofHit( stationId, segmentId, h_time, h_edep));
          //tofHitsMCTruth.push_back(ExtMonUCITofHitMCTruth(stationId, segmentId, h_time, h_edep, h_charged));
          tofHitsMCptr.push_back(mcptr_tof);

          if ( _diagLevel >= 3)
          { 
            std::cout << "MakeExtMonUCITofHits::makeExtMonUCITofHits : Add hit "
                      << stationId << " " << segmentId << " time " << h_time << " e " << h_edep
                      << " Create a new hit due to timeGap " 
                      << std::endl;
          }


          // ...and create new hit
          mcptr_tof.clear();
          mcptr_tof.push_back( art::Ptr<StepPointMC>( steps, tof_hits[i]._hit_id));

          h_time    = tof_hits[i]._time;
          h_edep    = tof_hits[i]._edep;
          h_edepc   = tof_hits[i]._edep_corr;
          //h_charged = tof_hits[i]._charged;

        } else {

          // Append data to hit
          h_edep  += tof_hits[i]._edep;
          h_edepc += tof_hits[i]._edep_corr;

          //if( tof_hits[i]._charged>0 ) h_charged = 1; // this does not count the charge...
          mcptr_tof.push_back( art::Ptr<StepPointMC>(steps,tof_hits[i]._hit_id));
        }
      }

      //if ( h_energy > 0.0 ) tofHits.push_back( ExtMonUCITofHit(stationId, segmentId, h_time,h_energy));
      tofHits.push_back(       ExtMonUCITofHit( stationId, segmentId, h_time, h_edep));
      tofHitsMCptr.push_back(mcptr_tof);

      if ( _diagLevel >= 3)
      {
        std::cout << "MakeExtMonUCITofHits::makeExtMonUCITofHits : Add hit "
                  << stationId << " " << segmentId << " time " << h_time << " e " << h_edep
                  << std::endl;      
      }

      typedef cet::map_vector_key key_type;
      typedef std::map<key_type, std::vector<int> > TrackHitMap;

      TrackHitMap aTrackHitMap;
      aTrackHitMap.clear();

      // Loop over steps inside the segment
      for( size_t i=0; i<isteps.size(); i++ ) {

        StepPointMC const& h = (*steps)[isteps[i]];
        vector<int> &steps_id = aTrackHitMap[h.trackId()];
        steps_id.push_back(isteps[i]);

        if( _diagLevel >= 3)
        {
          std::cout << "Track Id " << h.trackId()
                    << " Step " << isteps[i]
                    << std::endl;
        }
      }
      if ( _diagLevel >= 3)
      {
        std::cout << "MakeExtMonUCITofHits::makeExtMonUCITofHits : Track size " << aTrackHitMap.size()
                  << std::endl;
      }

      for(TrackHitMap::const_iterator aTrack = aTrackHitMap.begin(); aTrack != aTrackHitMap.end(); ++aTrack ) {

        // Tof Segment Id
        key_type trackId = aTrack->first;
        int firstStepRef = -9;

        double h_time = 9e27;
        double h_energy = 0.0;
        CLHEP::Hep3Vector h_pos(0.0, 0.0, 0.0);

        vector<int> const& isteps = aTrack->second;

        // Loop over steps associated with this track
        for( size_t i=0; i<isteps.size(); i++ ) {

          int hitRef = isteps[i];
          StepPointMC const &h = (*steps)[hitRef];

          double edep    = h.eDep(); // each track hit has its  energy deposit assigned to it

          if( _diagLevel >= 3)
          {
            std::cout << "Track hit ref " << hitRef << " trackId " << h.trackId() << " edep " << edep << " time " << h.time() << std::endl;
          }
          if( edep<0.0 || (edep==0.0 && !keepNoEdep) ) continue; // Do not create hit if there is no energy deposition
          h_energy += edep;

          double time = h.time();
          if ( time < h_time ) {
            h_time = time;
            firstStepRef = hitRef;
          }
        }

        if ( h_energy == 0.0 && !keepNoEdep ) continue;

        int trackId_int = -1;
        int pdgId = -1;
        CLHEP::Hep3Vector position;
        CLHEP::Hep3Vector momentum;
        CLHEP::Hep3Vector vertex;
        CLHEP::Hep3Vector vertexMomentum;
        float vertexTime = 0.0;
        int isPrimary = -1.0;
        int orgTrackId = -1;
        int orgPdgId   = -1;
        CLHEP::Hep3Vector orgVertex;
        CLHEP::Hep3Vector orgVertexMomentum;
        float orgTime = -1.0;

        if( _diagLevel >= 3)
        {
          std::cout << " Track : Id " << trackId << " edep " << h_energy << " time " << h_time
                    << " firstStepRef " << firstStepRef
                    << std::endl;
        }

        StepPointMC const &first_h = (*steps)[firstStepRef];
        if( _diagLevel >= 3)
        {
          std::cout << " Step  : Id " << first_h.trackId() << " VolumeId " << first_h.volumeId() 
                    << " position " << first_h.position() << " momentum " << first_h.momentum()
                    << std::endl;
        }

        trackId_int  = first_h.trackId().asInt();
        position = first_h.position();
        momentum = first_h.momentum();

        art::Ptr<SimParticle> const& simp = first_h.simParticle();
        if (simp) {
          if( _diagLevel >= 3)
          {
            std::cout << " SimParticle : Id " << simp->id() << " pdgId " << simp->pdgId()
                      << " isPrimary " << simp->isPrimary()
                      << " StartGlobalTime " << simp->startGlobalTime()
                      << " StartPosition " << simp->startPosition()
                      << " StartMomentum " << simp->startMomentum()
                      << std::endl;
          }

          pdgId  = simp->pdgId();
          vertex = simp->startPosition();
          vertexMomentum = simp->startMomentum();
          vertexTime = simp->startGlobalTime();
          isPrimary  = simp->isPrimary();

          art::Ptr<SimParticle> simp_parent;
          if (simp->isPrimary())
          {
            simp_parent = simp;
          }
          else 
          {
            simp_parent = simp->parent();
            while ( simp_parent != 0 && simp_parent->parent() != 0 ) 
            {
              simp_parent = simp_parent->parent();
            }
          }

          if (simp_parent) {
            if( _diagLevel >= 3)
            {
              std::cout << " Sim Parent : Id " << simp_parent->id() << " pdgId " << simp_parent->pdgId()
                        << " isPrimary " << simp_parent->isPrimary() 
                        << " StartGlobalTime " << simp_parent->startGlobalTime()
                        << " StartPosition " << simp_parent->startPosition()
                        << " StartMomentum " << simp_parent->startMomentum()
                        << std::endl;
            }

            orgTrackId  = simp_parent->id().asInt(); 
            orgPdgId    = simp_parent->pdgId();
            orgVertex   = simp_parent->startPosition();
            orgVertexMomentum = simp_parent->startMomentum();
            orgTime     = simp_parent->startGlobalTime(); 

          }
                            
          art::Ptr<GenParticle> const& genp = simp->genParticle();
          if (genp) {
            if( _diagLevel >= 3)
            {
              std::cout << " GenParticle : pdgId " << genp->pdgId() << " generateorId " << genp->generatorId() 
                      << " time " << genp->time() << " position " << genp->position()
                      << " momentum " << genp->momentum()
                      << std::endl;
            }
          }
        }

        ExtMonUCITofHitMCTruth aHitMCTruth( stationId, segmentId, h_time, h_energy, trackId_int, pdgId, 
                                            position, momentum, vertex, vertexMomentum, vertexTime,
                                            isPrimary, orgTrackId, orgPdgId, orgVertex, orgVertexMomentum, orgTime );
        if( _diagLevel >= 3)
        {
          std::cout << "MakeExtMonUCITofHits::makeExtMonUCITofHits : Add aHitMCTruth " << std::endl;
          std::cout << aHitMCTruth << std::endl;
        }

        tofHitsMCTruth.push_back( aHitMCTruth );

      }  // loop of iterator aTrack

    }

  }
  
} // namespace mu2e

using mu2e::MakeExtMonUCITofHits;
DEFINE_ART_MODULE(MakeExtMonUCITofHits);
