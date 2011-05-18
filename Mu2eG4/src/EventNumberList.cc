//
// Utility class to select listed events within G4.
//
// $Id: EventNumberList.cc,v 1.5 2011/05/18 02:27:17 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:17 $
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
    return inList(event->GetEventID());
  }

  bool EventNumberList::inList( int eventNumber) const{

    // Special case: if the list contains only event -1, then
    // always accept the event.
    if ( _eventNumbers.size() == 1 && _eventNumbers[0] == -1 ){
      return true;
    }

    // Normal behaviour.  Check to see if event number is in the list.
    for ( std::size_t i=0; i<_eventNumbers.size(); ++i ){
      if ( _eventNumbers[i] == eventNumber ) return true;
    }

    return false;
  }


}
