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
#include "GeneralUtilities/inc/splitLine.hh"

using namespace std;

namespace mu2e {

  namespace {

    // A helper class used to mark slots in the alive list as dead.
    class MarkAsDead{
    public:
      MarkAsDead( std::vector<bool>& alive):
        _alive(alive){
      }

      void operator()(Straw const& s){
        _alive.at(s.index().asInt()) = false;
      }

    private:
      std::vector<bool>& _alive;
    };

    // Fixme: the *FromString methods need to go into the TrackerGeom class.
    PanelId panelIdFromString ( std::string const& s ){
      vector<string> v;
      splitLine( s, "_", v);
      if ( v.size() != 2 ){
        throw cet::exception("CONFIG")
          << "panelIdFromString: expected two parts but found: "
          << v.size()
          << "\n";
      }

      istringstream sdev(v[0]);
      istringstream ssec(v[1]);
      int dev, sec;
      sdev >> dev;
      ssec >> sec;
      return PanelId(dev,sec);
    }

    LayerId layerIdFromString ( std::string const& s ){
      vector<string> v;
      splitLine( s, "_", v);
      if ( v.size() != 3 ){
        throw cet::exception("CONFIG")
          << "layerIdFromString: expected three parts but found: "
          << v.size()
          << "\n";
      }

      istringstream sdev(v[0]);
      istringstream ssec(v[1]);
      istringstream slay(v[2]);
      int dev, sec, lay;
      sdev >> dev;
      ssec >> sec;
      slay >> lay;
      return LayerId(dev,sec,lay);
    }

    StrawId strawIdFromString ( std::string const& s ){
      vector<string> v;
      splitLine( s, "_", v);
      if ( v.size() != 4 ){
        throw cet::exception("CONFIG")
          << "strawIdFromString: expected four parts but found: "
          << v.size()
          << "\n";
      }

      istringstream sdev(v[0]);
      istringstream ssec(v[1]);
      istringstream slay(v[2]);
      istringstream sstr(v[3]);
      int dev, sec, lay, str;
      sdev >> dev;
      ssec >> sec;
      slay >> lay;
      sstr >> str;
      return StrawId(dev,sec,lay,str);
    }

    // Out-of-class functions to deal with the parameter set work.
    // Kept out-of-class to hide implementation from the header.
    void addDeadPlanes( TTracker const& tracker,
                         fhicl::ParameterSet const& pset,
                         vector<bool>& alive,
                         bool verbosity ){

      vector<int> devs = pset.get<vector<int> >( "deadPlanes", vector<int>() );

      MarkAsDead marker(alive);

      for ( vector<int>::const_iterator i=devs.begin(), e=devs.end();
            i != e; ++i ){
        if ( verbosity > 0 ) cout << "Deadening straws in Plane: " << *i << endl;
        tracker.getPlane(*i).forAllStraws( marker );
      }

    }

    void addDeadPanels( TTracker const& tracker,
                         fhicl::ParameterSet const& pset,
                         vector<bool>& alive,
                         bool verbosity  ){

      vector<string> secs = pset.get<vector<string> >( "deadPanels", vector<string>() );
      vector<PanelId> secIds;

      for ( vector<string>::const_iterator i=secs.begin(), e=secs.end();
            i != e; ++i ){
        secIds.push_back( panelIdFromString(*i) );
      }

      MarkAsDead marker(alive);

      for ( vector<PanelId>::const_iterator i=secIds.begin(), e=secIds.end();
            i != e; ++i ){
        if ( verbosity > 0 ) cout << "Deadening straws in Panel: " << *i << endl;
        tracker.getPanel(*i).forAllStraws( marker );
      }

    }


    void addDeadLayers( TTracker const& tracker,
                        fhicl::ParameterSet const& pset,
                        vector<bool>& alive,
                         bool verbosity  ){

      vector<string> lays = pset.get<vector<string> >( "deadLayers", vector<string>() );
      vector<LayerId> layIds;

      for ( vector<string>::const_iterator i=lays.begin(), e=lays.end();
            i != e; ++i ){
        layIds.push_back( layerIdFromString(*i) );
      }

      MarkAsDead marker(alive);

      for ( vector<LayerId>::const_iterator i=layIds.begin(), e=layIds.end();
            i != e; ++i ){
        if ( verbosity > 0 ) cout << "Deadening straws in Layer: " << *i << endl;
        tracker.getLayer(*i).forAllStraws( marker );
      }

    }

    void addDeadStraws( TTracker const& tracker,
                        fhicl::ParameterSet const& pset,
                        vector<bool>& alive,
                        bool verbosity  ){

      vector<string> straws = pset.get<vector<string> >( "deadStraws", vector<string>() );
      vector<StrawId> strawIds;

      for ( vector<string>::const_iterator i=straws.begin(), e=straws.end();
            i != e; ++i ){
        strawIds.push_back( strawIdFromString(*i) );
      }

      for ( vector<StrawId>::const_iterator i=strawIds.begin(), e=strawIds.end();
            i != e; ++i ){
        if ( verbosity > 0 ) cout << "Deadening straw: " << *i << endl;
        Straw const& straw = tracker.getStraw(*i);
        alive.at(straw.index().asInt()) = false;
      }

    }

  } // end anonymous namespace

  DeadStrawList::DeadStrawList( fhicl::ParameterSet const& pset ):
    _verbosity(pset.get<int>("verbosity",0)),
    _alive(){
  }

  void DeadStrawList::reset( fhicl::ParameterSet const& pset ){

    TTracker const& tracker(*GeomHandle<TTracker>());

    // Default is that all straws are alive; assign implies clear.
    _alive.assign( tracker.getAllStraws().size(), true);

    // Parse the input to mark straws as dead.
    addDeadPlanes( tracker, pset, _alive, _verbosity );
    addDeadPanels( tracker, pset, _alive, _verbosity );
    addDeadLayers ( tracker, pset, _alive, _verbosity );
    addDeadStraws ( tracker, pset, _alive, _verbosity );

    if ( _verbosity > 1 ) {
      print(cout);
    }

  }

  void DeadStrawList::print( ostream& out) const{
    TTracker const& tracker(*GeomHandle<TTracker>());

    for ( deque<Straw>::const_iterator s=tracker.getAllStraws().begin(), se=tracker.getAllStraws().end();
          s != se; ++s ){

      if ( !_alive.at( s->index().asUint() ) ){
        out << "Dead straw: " << s->id() << endl;
      }

    }

  }

} // namespace mu2e
