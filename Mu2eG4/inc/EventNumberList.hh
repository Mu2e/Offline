#ifndef EventNumberList_HH
#define EventNumberList_HH
//
// Utility class to select listed events within G4.
//
// $Id: EventNumberList.hh,v 1.1 2009/11/11 14:40:00 kutschke Exp $
// $Author: kutschke $
// $Date: 2009/11/11 14:40:00 $
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
    void Add( int n, int const* eventNumbers);
    void Add( std::vector<int> const& eventNumbers );
    void Add( int eventNumber );

    // Is the current G4 event in the list?
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
#endif
