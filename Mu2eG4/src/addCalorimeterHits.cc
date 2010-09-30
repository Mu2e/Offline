//
// Populate output collection for calorimeter
//
// Original author Ivan Logashenko
//
#include <vector>
#include <map>

// Mu2e includes
#include "Mu2eG4/inc/StepPointG4.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "Mu2eG4/inc/addCalorimeterHits.hh"

// G4 includes
#include "G4Event.hh"
#include "G4SDManager.hh"

using namespace std;

namespace mu2e{

  // Utility class (structure) to hold calculated drift time for G4 hits

  class ROHit {
  public:

    int    _hit_id;
    double _edep;
    double _edep_corr;
    int    _charged;
    double _time;

    ROHit(int hit_id, double edep, double edep1, int charged, double time):
      _hit_id(hit_id), _edep(edep), _edep_corr(edep1), 
      _charged(charged), _time(time) { }

    // This operator is overloaded in order to time-sort the hits 
    bool ROHit::operator <(const ROHit& b) const { return (_time < b._time); }

  };

  void addCalorimeterHits ( const G4Event *g4event, 
			    CaloHitCollection &outputHits, 
			    CaloHitMCTruthCollection& mchits ) {

    // Get calorimeter geometry description
    edm::Service<GeometryService> geom;
    if( ! geom->hasElement<Calorimeter>() ) return;
    GeomHandle<Calorimeter> cg;

    // G4 Hit collections for this event.

    G4HCofThisEvent* hce = g4event->GetHCofThisEvent();
    
    // Get the collection ID for the VD hits.

    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    G4int colId = SDman->GetCollectionID("CaloCollection");
    
    // If there are not hits - do nothing

    if ( colId < 0 || hce == 0 ) return;
      
    StepPointG4Collection* hits = static_cast<StepPointG4Collection*>(hce->GetHC(colId));
    G4int nHits = hits->entries();

    if( nHits<=0 ) return;

    double length     = cg->crystalHalfLength();
    double nonUniform = cg->getNonuniformity();
    double timeGap    = cg->getTimeGap();
    double addEdep    = cg->getElectronEdep();
    int    nro        = cg->nROPerCrystal();

    // Organize hits by readout elements

    typedef std::map<int,std::vector<int> > HitMap;
    HitMap hitmap;
    for ( int i=0; i<hits->entries(); ++i){
      StepPointG4* h = (*hits)[i];
      vector<int> &hits_id = hitmap[h->volumeId()];
      hits_id.push_back(i);
    }

    // Loop over all readout elements

    vector<ROHit> ro_hits;

    for(HitMap::const_iterator ro = hitmap.begin(); ro != hitmap.end(); ++ro ) {

      // readout ID
      int roid = ro->first;

      // Prepare info for hit creation
      ro_hits.clear();

      // Loop over all hits found for this readout element

      vector<int> const& ihits = ro->second;

      for( size_t i=0; i<ihits.size(); i++ ) {

        int hitRef = ihits[i];
	StepPointG4* h = (*hits)[hitRef];
        CLHEP::Hep3Vector const& pos = h->position();
        double edep    = h->eDep()/nro;
        double hitTime = h->time();
	int charged = 0;

	// Calculate correction for edep; edep<0 means energy deposition in APD, not in crystal

	double edep_corr=0;
	if( edep>=0 ) {
	  edep_corr = edep * (1+(pos.z()/length)*nonUniform/2);
	} else { 
	  charged=1; 
	  edep=0; 
	}
	
	ro_hits.push_back(ROHit(hitRef,edep,edep_corr,charged,hitTime));

      }

      // Sort this hits in time
      sort(ro_hits.begin(), ro_hits.end());

      // Loop over sorted hits and form complete calorimeter hits

      double h_time    = ro_hits[0]._time;
      double h_edep    = ro_hits[0]._edep;
      double h_edepc   = ro_hits[0]._edep_corr;
      int    h_charged = ro_hits[0]._charged;

      for( size_t i=1; i<ro_hits.size(); ++i ) {
	if( (ro_hits[i]._time-ro_hits[i-1]._time) > timeGap ) {
	  // Save current hit
	  outputHits.push_back(CaloHit(roid,h_time,h_edepc+h_charged*addEdep));
	  mchits.push_back(CaloHitMCTruth(roid,h_time,h_edep,h_charged));
	  // ...and create new hit	  
	  h_time    = ro_hits[i]._time;
	  h_edep    = ro_hits[i]._edep;
	  h_edepc   = ro_hits[i]._edep_corr;
	  h_charged = ro_hits[i]._charged;
	} else {
	  // Append data to hit
	  h_edep  += ro_hits[i]._edep;
	  h_edepc += ro_hits[i]._edep_corr;
	  if( ro_hits[i]._charged>0 ) h_charged = 1;
	}
      }

      /*
      cout << "CaloHit: id=" << roid 
	   << " time=" << h_time
	   << " trueEdep=" << h_edep
	   << " corrEdep=" << h_edepc
	   << " chargedEdep=" << (h_edepc+h_charged*addEdep)
	   << endl;
      */

      outputHits.push_back(CaloHit(roid,h_time,h_edepc+h_charged*addEdep));
      mchits.push_back(CaloHitMCTruth(roid,h_time,h_edep,h_charged));
      
    }

  }

} // end namespace mu2e
