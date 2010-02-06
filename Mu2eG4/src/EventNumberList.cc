//
// Utility class to select listed events within G4.
//
// $Id: EventNumberList.cc,v 1.2 2010/02/06 19:40:31 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/02/06 19:40:31 $
//
// Original author Rob Kutschke
//

#include "Mu2eG4/inc/EventNumberList.hh"

#include "G4RunManager.hh"

namespace mu2e {

  EventNumberList::EventNumberList( int eventNumber ):
    _eventNumbers(1,eventNumber)
  {
  }

  EventNumberList::EventNumberList( int n, int const* eventNumbers):
    _eventNumbers()
  {
    add(n, eventNumbers);
  }

  EventNumberList::EventNumberList( std::vector<int> const& eventNumbers ):
    _eventNumbers(eventNumbers)
  {
  }

  void EventNumberList::add( int n, int const* eventNumbers){
    if ( n < 0 ) {
      // print warning.
      return;
    }

    for ( int i=0; i<n; ++i){
      _eventNumbers.push_back(eventNumbers[i]);
    }
  }

  void EventNumberList::add( std::vector<int> const& eventNumbers )
  {
    add( eventNumbers.size(), &eventNumbers[0] );
  }

  void EventNumberList::add( int eventNumber){
    _eventNumbers.push_back(eventNumber);
  }

  bool EventNumberList::inList() const{
    G4Event const* event = G4RunManager::GetRunManager()->GetCurrentEvent();
    inList(event->GetEventID());
  }

  bool EventNumberList::inList( int eventNumber) const{
    for ( std::size_t i=0; i<_eventNumbers.size(); ++i ){
      if ( _eventNumbers[i] == eventNumber ) return true;
    }
    return false;
  }
  
  
}
