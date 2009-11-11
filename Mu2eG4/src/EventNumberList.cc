//
// Utility class to select listed events within G4.
//
// $Id: EventNumberList.cc,v 1.1 2009/11/11 14:40:00 kutschke Exp $
// $Author: kutschke $
// $Date: 2009/11/11 14:40:00 $
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
    Add(n, eventNumbers);
  }

  EventNumberList::EventNumberList( std::vector<int> const& eventNumbers ):
    _eventNumbers(eventNumbers)
  {
  }

  void EventNumberList::Add( int n, int const* eventNumbers){
    if ( n < 0 ) {
      // print warning.
      return;
    }

    for ( int i=0; i<n; ++i){
      _eventNumbers.push_back(eventNumbers[i]);
    }
  }

  void EventNumberList::Add( std::vector<int> const& eventNumbers )
  {
    Add( eventNumbers.size(), &eventNumbers[0] );
  }

  void EventNumberList::Add( int eventNumber){
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
