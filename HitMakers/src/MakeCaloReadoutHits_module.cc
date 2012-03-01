//
// An EDProducer Module that reads StepPointMC objects and turns them into
// CaloHit, CaloHitMCTruth, CaloHitMCPtr, CrystalOnlyHit, CrystalHitMCPtr
// objects.
//
// Original author Ivan Logashenko; major changes by Krzysztof Genser and Rob Kutschke.
//
// Notes
// 1) The CrystalOnlyHitCollection is a form of MC truth on a per crystal basis.
//    It represents an idea per crystal response if no readouts were hit directly.
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
#include "art/Framework/Principal/Selector.h"
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
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/sort_functors.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "MCDataProducts/inc/CaloHitMCTruthCollection.hh"
#include "MCDataProducts/inc/CaloCrystalOnlyHitCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"

// Other includes.
#include "CLHEP/Vector/ThreeVector.h"

using namespace std;

namespace mu2e {

  // Anonymous namespace to hold some helper classes.
  namespace {

    // A helper class to hold temporary information.
    class ROHit {
    public:

      // Is this StepPointMC from a hit in the crystal or a hit in the readout device?
      enum step_type {crystal, readout};

      art::Ptr<StepPointMC> _step;    // Ptr back to the StepPointMC.
      double    _edep;                // copy of energy from the StepPointMC.
      double    _edep_corr;           // Energy corrected for absorption within the crystal.
      step_type _type;                // Is this a hit in the crystal or the readout?
      double    _time;                // copy of the time from the StepPointMC.

      ROHit(art::Ptr<StepPointMC> const& step, double edep, double edep1, step_type charged, double time):
        _step(step), _edep(edep), _edep_corr(edep1),
        _type(charged), _time(time) { }

      // This operator is overloaded in order to time-sort the hits
      bool operator <(const ROHit& b) const { return (_time < b._time); }

    }; // end class ROHit

    // A helper class to add Ptr's to the appropriate PtrStepPointMCVector collection.
    class PtrAdder{
    public:
      PtrAdder( PtrStepPointMCVector& crystals, 
                PtrStepPointMCVector& readouts ):
        _crystals(crystals),
        _readouts(readouts){}

      // Choose the appropriate collection and add the Ptr for this hit.
      void operator()( ROHit const& hit ){
        if ( hit._type == ROHit::crystal){
          _crystals.push_back( hit._step );
        } else {
          _readouts.push_back( hit._step );
        }
      }

    private:

      // The two possible output collections.
      PtrStepPointMCVector& _crystals;
      PtrStepPointMCVector& _readouts;
                   
    }; // end class PtrAdder

  } // end anonymous namespace

  //--------------------------------------------------------------------
  //
  //
  class MakeCaloReadoutHits : public art::EDProducer {
  public:

    // First vector is list of crystal steps, associated with parcular readout element.
    // Second vector is list of readout steps, associated with parcular readout element.
    typedef art::Ptr<StepPointMC> StepPtr;
    typedef std::vector<StepPtr>  StepPtrs;
    typedef std::map<int,std::pair<StepPtrs,StepPtrs> > HitMap;

    explicit MakeCaloReadoutHits(fhicl::ParameterSet const& pset) :

      // Parameters
      _diagLevel(pset.get<int>("diagLevel",0)),
      _maxFullPrint(pset.get<int>("maxFullPrint",5)),
      _stepPoints(pset.get<string>("calorimeterStepPoints","calorimeter")),
      _rostepPoints(pset.get<string>("calorimeterROStepPoints","calorimeterRO")),
      _g4ModuleLabel(pset.get<string>("g4ModuleLabel")),

      _messageCategory("CaloReadoutHitsMaker"){

      // Tell the framework what we make.
      produces<CaloHitCollection>();
      produces<CaloHitMCTruthCollection>();
      produces<CaloCrystalOnlyHitCollection>();
      produces<PtrStepPointMCVectorCollection>("CaloHitMCCrystalPtr");
      produces<PtrStepPointMCVectorCollection>("CaloHitMCReadoutPtr");

    }
    virtual ~MakeCaloReadoutHits() { }

    virtual void beginJob();

    void produce( art::Event& e);

  private:

    typedef std::vector< art::Handle<StepPointMCCollection> > HandleVector;

    // Diagnostics level.
    int _diagLevel;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Name of the StepPoint collection
    std::string _stepPoints;
    std::string _rostepPoints;

    // Parameters
    string _g4ModuleLabel;  // Name of the module that made these hits.

    // A category for the error logger.
    const std::string _messageCategory;

    void makeCalorimeterHits (HandleVector const& crystalStepsHandles,
                              HandleVector const& readoutStepsHandles,
                              CaloHitCollection &,
                              CaloHitMCTruthCollection&,
                              CaloCrystalOnlyHitCollection&,
                              PtrStepPointMCVectorCollection&,
                              PtrStepPointMCVectorCollection& );

    // Fill the hitmap, using the readout Id as the key.
    void fillMapByReadoutId( int nro,
                             HitMap& map,
                             HandleVector const& crystalStepsHandles,
                             HandleVector const& readoutStepsHandles);

    // Fill the hitmap, using the crystal Id as the key; used for making the per crystal MC truth.
    void fillMapByCrystalId( HitMap& map,
                             HandleVector const& crystalStepsHandles);

    // Print information about the data products found by the selector functions. 
    void printDataProductInfo( HandleVector const& crystalStepsHandles,
                               HandleVector const& readoutStepsHandles );

  };

  void MakeCaloReadoutHits::beginJob(){
  }

  void
  MakeCaloReadoutHits::produce(art::Event& event) {

    if ( _diagLevel > 0 ) cout << "MakeCaloReadoutHits: produce() begin" << endl;

    static int ncalls(0);
    ++ncalls;

    // Check that calorimeter geometry description exists
    art::ServiceHandle<GeometryService> geom;
    if( ! geom->hasElement<Calorimeter>() ) return;

    // A container to hold the output hits.
    auto_ptr<CaloHitCollection>               caloHits         (new CaloHitCollection);
    auto_ptr<CaloHitMCTruthCollection>        caloMCHits       (new CaloHitMCTruthCollection);
    auto_ptr<CaloCrystalOnlyHitCollection>    caloCrystalMCHits(new CaloCrystalOnlyHitCollection);
    auto_ptr<PtrStepPointMCVectorCollection>  caloMCptrHits    (new PtrStepPointMCVectorCollection);
    auto_ptr<PtrStepPointMCVectorCollection>  caloMCroptrHits  (new PtrStepPointMCVectorCollection);

    // These selectors will select data products with the given instance name, and ignore
    // all other fields of the product ID.
    art::ProductInstanceNameSelector getCrystalSteps(_stepPoints);
    art::ProductInstanceNameSelector getReadoutSteps(_rostepPoints);

    // Get the StepPointMCs from the event.
    HandleVector crystalStepsHandles, readoutStepsHandles;
    event.getMany( getCrystalSteps, crystalStepsHandles);
    event.getMany( getReadoutSteps, readoutStepsHandles);

    // Informational message on the first event.
    static bool firstEvent(true);
    if ( firstEvent ) {
      printDataProductInfo( crystalStepsHandles, readoutStepsHandles);
      firstEvent = false;
    }

    // Do the real work.
    makeCalorimeterHits(crystalStepsHandles, readoutStepsHandles,
                        *caloHits, *caloMCHits, *caloCrystalMCHits,
                        *caloMCptrHits, *caloMCroptrHits);

    if ( ncalls < _maxFullPrint && _diagLevel > 2 ) {
      cout << "MakeCaloReadoutHits: Total number of calorimeter hits = "
           << caloHits->size()
           << endl;
      cout << "MakeCaloReadoutHits: Total number of crystal MC hits = "
           << caloCrystalMCHits->size()
           << endl;
    }

    // Add the output hit collection to the event
    event.put(caloHits);
    event.put(caloMCHits);
    event.put(caloCrystalMCHits);
    event.put(caloMCptrHits,"CaloHitMCCrystalPtr");
    event.put(caloMCroptrHits,"CaloHitMCReadoutPtr");

    if ( _diagLevel > 0 ) cout << "MakeCaloReadoutHits: produce() end" << endl;

  } // end of ::analyze.

  void MakeCaloReadoutHits::makeCalorimeterHits (HandleVector const& crystalStepsHandles,
                                                 HandleVector const& readoutStepsHandles,
                                                 CaloHitCollection &caloHits,
                                                 CaloHitMCTruthCollection& caloHitsMCTruth,
                                                 CaloCrystalOnlyHitCollection& caloCrystalHitsMCTruth,
                                                 PtrStepPointMCVectorCollection& caloHitsMCCrystalPtr,
                                                 PtrStepPointMCVectorCollection& caloHitsMCReadoutPtr ){
    
    // Get calorimeter geometry description
    Calorimeter const & cal = *(GeomHandle<Calorimeter>());

    double length     = cal.crystalHalfLength();
    double nonUniform = cal.getNonuniformity();
    double timeGap    = cal.getTimeGap();
    double addEdep    = cal.getElectronEdep();
    int    nro        = cal.nROPerCrystal();

    // Organize steps by readout elements.
    HitMap hitmap;
    fillMapByReadoutId( nro, hitmap, crystalStepsHandles, readoutStepsHandles);

    // Loop over all readout elements to form ro hits

    vector<ROHit> ro_hits;

    for(HitMap::const_iterator ro = hitmap.begin(); ro != hitmap.end(); ++ro ) {

      // readout ID
      int roid = ro->first;

      // Prepare info for hit creation before the next iteration
      ro_hits.clear();

      // Loop over all hits found for this readout element
      StepPtrs const& isteps   = ro->second.first;
      StepPtrs const& irosteps = ro->second.second;

      // Loop over steps inside the crystal
      for ( StepPtrs::const_iterator i=isteps.begin(), e=isteps.end();
            i != e; ++i ){

        StepPointMC const& h = **i;

        double edep    = h.eDep()/nro; // each ro has its crystal energy deposit assigned to it
        if( edep<=0.0 ) continue; // Do not create hit if there is no energy deposition

        // Hit position in Mu2e frame
        CLHEP::Hep3Vector const& pos = h.position();
        // Hit position in local crystal frame
        CLHEP::Hep3Vector posLocal = cal.toCrystalFrame(roid,pos);
        // Calculate correction for edep
        double edep_corr = edep * (1.0+(posLocal.z()/length)*nonUniform/2.0);

        ro_hits.push_back(ROHit(*i,edep,edep_corr,ROHit::crystal,h.time()));

      }

      // Loop over steps inside the readout (direct energy deposition in APD)
      for( StepPtrs::const_iterator i=irosteps.begin(), e=irosteps.end(); 
           i!=e; ++i ){

        StepPointMC const& h = **i;

        // There is no cut on energy deposition here - may be, we need to add one?
        ro_hits.push_back(ROHit(*i,0.,0.,ROHit::readout,h.time()));

      }

      // Sort hits by time
      sort(ro_hits.begin(), ro_hits.end());

      // A buffer to collect output.
      PtrStepPointMCVector mcptr_crystal;
      PtrStepPointMCVector mcptr_readout;

      // A tool to add the art::Ptr to the correct output collection.
      PtrAdder addPtr( mcptr_crystal, mcptr_readout );

      // Loop over sorted hits and form complete ro/calorimeter hits
      double h_time    = ro_hits[0]._time;
      double h_edep    = ro_hits[0]._edep;
      double h_edepc   = ro_hits[0]._edep_corr;
      ROHit::step_type h_type = ro_hits[0]._type;
      addPtr( ro_hits[0] );

      for( size_t i=1; i<ro_hits.size(); ++i ) {
        if( (ro_hits[i]._time-ro_hits[i-1]._time) > timeGap ) {
          // Save current hit
          caloHits.push_back(       CaloHit(       roid,h_time,h_edepc+h_type*addEdep));
          caloHitsMCTruth.push_back(CaloHitMCTruth(roid,h_time,h_edep,h_type));
          caloHitsMCCrystalPtr.push_back(mcptr_crystal);
          caloHitsMCReadoutPtr.push_back(mcptr_readout);
          // ...and create new hit
          mcptr_crystal.clear();
          mcptr_readout.clear();
          addPtr( ro_hits[i] );

          h_time    = ro_hits[i]._time;
          h_edep    = ro_hits[i]._edep;
          h_edepc   = ro_hits[i]._edep_corr;
          h_type = ro_hits[i]._type;
        } else {
          // Append data to hit
          h_edep  += ro_hits[i]._edep;
          h_edepc += ro_hits[i]._edep_corr;
          if( ro_hits[i]._type != ROHit::crystal ) h_type = ROHit::readout; // this does not count the charge...
          addPtr( ro_hits[i] );
        }
      }

      caloHits.push_back(       CaloHit(roid,h_time,h_edepc+h_type*addEdep));
      caloHitsMCTruth.push_back(CaloHitMCTruth(roid,h_time,h_edep,h_type));
      caloHitsMCCrystalPtr.push_back(mcptr_crystal);
      caloHitsMCReadoutPtr.push_back(mcptr_readout);

    }

    // now CaloCrystalOnlyHit
    // Organize hits by crystal id; reject all but first RO; reject all RO hit directly
    hitmap.clear();
    fillMapByCrystalId( hitmap, crystalStepsHandles);

    // Working space.
    CaloCrystalOnlyHitCollection cr_hits;

    for(HitMap::const_iterator cr = hitmap.begin(); cr != hitmap.end(); ++cr) {

      // crystal id
      int cid = cr->first;

      // Prepare info for hit creation before the next iteration
      cr_hits.clear();

      for( size_t i=0; i<cr->second.first.size(); i++ ) {
        StepPointMC const& h2 = *cr->second.first[i];
        cr_hits.push_back(CaloCrystalOnlyHit(cid,h2.time(),h2.eDep()));
      }

      sort(cr_hits.begin(), cr_hits.end(), lessByTime<CaloCrystalOnlyHitCollection::value_type>());

      // now form final hits if they are close enough in time

      CaloCrystalOnlyHitCollection::value_type cHitMCTruth  = cr_hits[0];

      for ( size_t i=1; i<cr_hits.size(); ++i ) {

        if ( (cr_hits[i].time()-cr_hits[i-1].time()) > timeGap ) {

          // Save current hit

          caloCrystalHitsMCTruth.push_back(cHitMCTruth);

          // ...and create new hit

          cHitMCTruth  = cr_hits[i];

        } else {

          // Add energy to the old hit (keep the "earlier" time)

          cHitMCTruth.setEnergyDep(cHitMCTruth.energyDep()+cr_hits[i].energyDep());

        }

      }

      caloCrystalHitsMCTruth.push_back(cHitMCTruth);

    }


  } // end makeCalorimeterHits

  // Fill the hitmap using the readout id as the map key.
  void MakeCaloReadoutHits::fillMapByReadoutId( int nro,
                                                HitMap& hitmap,
                                                HandleVector const& crystalStepsHandles,
                                                HandleVector const& readoutStepsHandles){

    // Loop over hits from inside the crystal; make one entry for each readout attached to the crystal.
    for ( HandleVector::const_iterator i=crystalStepsHandles.begin(), e=crystalStepsHandles.end();
          i != e; ++i ){
      art::Handle<StepPointMCCollection> const& handle(*i);
      StepPointMCCollection const& steps(*handle);
      StepPointMCCollection::const_iterator j0=steps.begin();
      for ( StepPointMCCollection::const_iterator j=j0, je=steps.end(); j != je; ++j ){
        StepPointMC const& step(*j);
        for( int ro=0; ro<nro; ++ro ) {
          hitmap[step.volumeId()*nro+ro].first.push_back( StepPtr(handle,j-j0) );
        }
      }
    }

    // Loop over hits created in the readouts themselves.
    for ( HandleVector::const_iterator i=readoutStepsHandles.begin(), e=readoutStepsHandles.end();
          i != e; ++i ){
      art::Handle<StepPointMCCollection> const& handle(*i);
      StepPointMCCollection const& steps(*handle);
      StepPointMCCollection::const_iterator j0=steps.begin();
      for ( StepPointMCCollection::const_iterator j=j0, je=steps.end(); j != je; ++j ){
        StepPointMC const& step(*j);
        hitmap[step.volumeId()].second.push_back( StepPtr(handle,j-j0));
      }
    }

  }  // end fillMapByReadoutId


  // Fill the hitmap using the crystal id as the map key.
  void MakeCaloReadoutHits::fillMapByCrystalId( HitMap& hitmap,

                                                HandleVector const& crystalStepsHandles){
    for ( HandleVector::const_iterator i=crystalStepsHandles.begin(), e=crystalStepsHandles.end();
          i != e; ++i ){
      art::Handle<StepPointMCCollection> const& handle(*i);
      StepPointMCCollection const& steps(*handle);

      StepPointMCCollection::const_iterator j0=steps.begin();
      for ( StepPointMCCollection::const_iterator j=j0, je=steps.end(); j != je; ++j ){
        StepPointMC const& step(*j);
        hitmap[step.volumeId()].first.push_back( StepPtr(handle,j-j0) );
      }
    }

  } // end fillMapByCrystalId

  void MakeCaloReadoutHits::printDataProductInfo( HandleVector const& crystalStepsHandles,
                                                  HandleVector const& readoutStepsHandles ){
    mf::LogInfo log(_messageCategory);
    log << "MakeCaloReadoutHit::produce will use StepPointMCs from: \n";
    for ( HandleVector::const_iterator i=crystalStepsHandles.begin(), e=crystalStepsHandles.end();
          i != e; ++i ){
      art::Provenance const& prov(*(i->provenance()));
      log  << "   " << prov.branchName() << "\n";
    }
    log << "\nMakeCaloReadoutHit::produce will use StepPointMCs from: \n";
    for ( HandleVector::const_iterator i=readoutStepsHandles.begin(), e=readoutStepsHandles.end();
          i != e; ++i ){
      art::Provenance const& prov(*(i->provenance()));
      log  << "   " << prov.branchName() << "\n";
    }
  } // end printDataProductInfo


} // end namespace mu2e

using mu2e::MakeCaloReadoutHits;
DEFINE_ART_MODULE(MakeCaloReadoutHits);
