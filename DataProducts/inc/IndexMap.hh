#ifndef IndexMap_hh
#define IndexMap_hh

//
// Maps an index from a full collection onto the index in a condensed collection
//
// Original author Andrew Edmonds

#include <ostream>
#include <map>

namespace mu2e {

  typedef uint16_t FullIndex;
  typedef uint16_t CondensedIndex;

  class IndexMap{

  public:
    IndexMap() { };

    void addElement(FullIndex full, CondensedIndex condensed) {
      std::pair<FullIndex, CondensedIndex> newPair(full, condensed);
      _theMap.insert(newPair);
    }

    mu2e::CondensedIndex getCondensedIndex(const FullIndex& full) const {
      return _theMap.at(full);
    }

    bool checkInMap(const FullIndex& full) {
      if (_theMap.find(full) != _theMap.end()) {
        return true;
      }
      else {
        return false;
      }
    }

    // Print the information found in this hit.
    void print( std::ostream& ost, bool doEndl ) const {

      ost << "IndexMap: " << std::endl;
      for (const auto& element : _theMap) {
        ost << element.first << " --> " << element.second << std::endl;
      }
      if ( doEndl ){
        ost << std::endl;
      }
    }
    auto const& map() const { return _theMap; }

  private:
    std::map<FullIndex, CondensedIndex> _theMap; // FullIndex is first, CondensedIndex is second
  };



  inline std::ostream& operator<<( std::ostream& ost,
                                     IndexMap const& map){
    map.print(ost,false);
    return ost;
  }

}
#endif /* IndexMap_hh */
