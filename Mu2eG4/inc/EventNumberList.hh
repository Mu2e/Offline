#ifndef Mu2eG4_EventNumberList_hh
#define Mu2eG4_EventNumberList_hh
//
// Utility class to select listed events within G4.
//
// $Id: EventNumberList.hh,v 1.4 2011/05/18 02:27:17 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:17 $
//
// Original author Rob Kutschke
//

#include <vector>

namespace mu2e {

  class EventNumberList{

  public:

    // Default c'tor makes an empty list.
    EventNumberList():
      _eventNumbers(){
    }

    // Constructors take an event number of a list of event numbers.
    explicit EventNumberList( int eventNumber );

    EventNumberList( int n, int const* eventNumbers );
    EventNumberList( std::vector<int> const& eventNumbers );

    // Accept compiler generated:
    // copy c'tor, d'tor and assignment operator.

    // Extend the list of event numbers.
    void add( int n, int const* eventNumbers);
    void add( std::vector<int> const& eventNumbers );
    void add( int eventNumber );

    // Is the number of the current G4 event in the list?
    bool inList() const;

    // Is the specified event number in the list?
    bool inList(int eventNumber ) const;


    // Access to the size of the list.
    std::size_t size() const {
      return _eventNumbers.size();
    }

    // Access to the list.
    std::vector<int> const& getEventNumbers() const{
      return _eventNumbers;
    }


  public:

    // The list of event numbers.
    std::vector<int> _eventNumbers;

  };
}
#endif /* Mu2eG4_EventNumberList_hh */
