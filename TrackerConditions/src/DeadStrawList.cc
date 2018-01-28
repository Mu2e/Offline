//
// First version of a hit as described by Mu2e-doc-900.
//
// Original author Rob Kutschke
//

// C++ includes
#include <iostream>

// Mu2e includes
#include "TrackerConditions/inc/DeadStrawList.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "DataProducts/inc/StrawId.hh"
#include <sstream>

using namespace std;

namespace mu2e {

  namespace {

    // A helper class used to mark slots in the alive list as dead.
    class MarkAsDead{
    public:
      MarkAsDead( std::set<DeadStrawRange>& dead):
        _dead(dead){
      }

      void operator()(Straw const& s){
        _dead.insert(DeadStrawRange(s));
      }

    private:
      std::set<DeadStrawRange>& _dead;
    };

    // Out-of-class functions to deal with the parameter set work.
    // Kept out-of-class to hide implementation from the header.
    void addDeadPlanes( TTracker const& tracker,
                         fhicl::ParameterSet const& pset,
                         set<DeadStrawRange>& dead,
                         bool verbosity ){

      vector<int> devs = pset.get<vector<int> >( "deadPlanes", vector<int>() );

      MarkAsDead marker(dead);

      for ( vector<int>::const_iterator i=devs.begin(), e=devs.end();
            i != e; ++i ){
        if ( verbosity > 0 ) cout << "Deadening straws in Plane: " << *i << endl;
        tracker.getPlane(*i).forAllStraws( marker );
      }

    }

    void addDeadPanels( TTracker const& tracker,
                         fhicl::ParameterSet const& pset,
                         set<DeadStrawRange>& dead,
                         bool verbosity  ){

      vector<string> secs = pset.get<vector<string> >( "deadPanels", vector<string>() );
      vector<PanelId> secIds;

      for ( vector<string>::const_iterator i=secs.begin(), e=secs.end();
            i != e; ++i ){
        secIds.push_back( PanelId(*i) );
      }

      MarkAsDead marker(dead);

      for ( vector<PanelId>::const_iterator i=secIds.begin(), e=secIds.end();
            i != e; ++i ){
        if ( verbosity > 0 ) cout << "Deadening straws in Panel: " << *i << endl;
        tracker.getPanel(*i).forAllStraws( marker );
      }

    }

    // fixme: rewrite using panels if needed
    // void addDeadLayers( TTracker const& tracker,
    //                     fhicl::ParameterSet const& pset,
    //                     set<DeadStrawRange>& dead,
    //                      bool verbosity  ){

    //   vector<string> lays = pset.get<vector<string> >( "deadLayers", vector<string>() );
    //   vector<LayerId> layIds;

    //   for ( vector<string>::const_iterator i=lays.begin(), e=lays.end();
    //         i != e; ++i ){
    //     layIds.push_back( LayerId(*i) );
    //   }

    //   MarkAsDead marker(dead);

    //   for ( vector<LayerId>::const_iterator i=layIds.begin(), e=layIds.end();
    //         i != e; ++i ){
    //     if ( verbosity > 0 ) cout << "Deadening straws in Layer: " << *i << endl;
    //     tracker.getLayer(*i).forAllStraws( marker );
    //   }

    // }

    void addDeadStraws( TTracker const& tracker,
                        fhicl::ParameterSet const& pset,
                        set<DeadStrawRange>& dead,
                        bool verbosity  ){

      vector<string> straws = pset.get<vector<string> >( "deadStraws", vector<string>() );
      vector<StrawId> strawIds;

      for ( vector<string>::const_iterator i=straws.begin(), e=straws.end();
            i != e; ++i ){
        strawIds.push_back( StrawId(*i) );
      }

      for ( vector<StrawId>::const_iterator i=strawIds.begin(), e=strawIds.end();
            i != e; ++i ){
        if ( verbosity > 0 ) cout << "Deadening straw: " << *i << endl;
        Straw const& straw = tracker.getStraw(*i);
        dead.insert(DeadStrawRange(straw));
      }

    }

    void addPartlyDeadStraws( TTracker const& tracker,
                        fhicl::ParameterSet const& pset,
                        set<DeadStrawRange>& dead,
                        bool verbosity  ){

      vector<string> dstraws = pset.get<vector<string> >( "PartlyDeadStraws", vector<string>() );

      for ( auto dstring : dstraws){
// split the string into a StrawId (underscore delimited) and the FP range
	std::istringstream dstrings(dstring);
	double range(0.0);
	std::string sidname;
	dstrings >> sidname >> range;
	// check
	if(range == 0.0)
        throw cet::exception("CONFIG")
          << "DeadStrawList: expected StrawId and range but got "
	  << dstring << endl;
	StrawId sid(sidname);
        Straw const& straw = tracker.getStraw(sid);
        dead.insert(DeadStrawRange(straw,range));
      }
    }

  } // end anonymous namespace

  DeadStrawList::DeadStrawList( fhicl::ParameterSet const& pset ):
    _verbosity(pset.get<int>("verbosity",0)){
    // initialization is done in beginRun of user module
  }

  void DeadStrawList::reset( fhicl::ParameterSet const& pset ){

    TTracker const& tracker(*GeomHandle<TTracker>());

    // first, mark completely dead straws
    addDeadPlanes( tracker, pset, _deadstraws, _verbosity );
    addDeadPanels( tracker, pset, _deadstraws, _verbosity );
    //    addDeadLayers( tracker, pset, _deadstraws, _verbosity );
    addDeadStraws( tracker, pset, _deadstraws, _verbosity );

    // Then, add dead straw ranges.  set semantics insures a
    // partially-dead straw will never superscede a fully dead one

    addPartlyDeadStraws(tracker, pset, _deadstraws, _verbosity);


    if ( _verbosity > 0 ) {
      print(cout);
    }

    if(_verbosity > 1) {
      const auto & allstraws = tracker.getAllStraws();
      // for(const auto & straw: tracker.getAllStraws())
      for (size_t i = 0; i<tracker.nStraws(); ++i) {
        const Straw & straw = allstraws[i];
	std::cout << "Straw Index " << straw.index() << " Id " << straw.id() << endl;
      }
    }

  }

  bool DeadStrawList::isDead ( StrawIndex ind,double hitpos) const {
    bool retval(false);
    auto ifnd = _deadstraws.find(DeadStrawRange(ind));
    if(ifnd != _deadstraws.end()){
      if(ifnd->_range > 0)
	retval = fabs(hitpos) < ifnd->_range;
      else
	retval = fabs(hitpos) > -ifnd->_range;
    }
    return retval;
  }

  void DeadStrawList::print( ostream& out) const{
    TTracker const& tracker(*GeomHandle<TTracker>());

    for( auto idead: _deadstraws) {
      out << "Straw " << tracker.getStraw(idead._strawind).id() 
      << " is dead for distances < " << idead._range << endl;
    }

  }

} // namespace mu2e
