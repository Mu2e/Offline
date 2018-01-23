// This class implements cache logic that works for overlapping maps.
//
// Andrei Gaponenko, 2012

#ifndef BFCacheManager_hh
#define BFCacheManager_hh

#include <vector>
#include <map>
#include <cassert>

#include "CLHEP/Vector/ThreeVector.h"

#include "BFieldGeom/inc/BFMap.hh"

namespace mu2e {

  // There are two classes of magnetic field maps in Mu2e: "Inner" and "Outer" maps.
  // No overlaps are allowed among any of the maps in the "Inner" set.
  // The maps in the "Outer" set may overlap with the "Inner" maps, and among themselves.
  // If a space point belongs to an "Inner" map, that map will be used to
  // compute the field value.   If a point is outside of the "Inner" map set,
  // then the "Outer" map list will be consulted in order, and the first map
  // that contains the point will be used.
  //

  class BFCacheManager {

    typedef std::vector<const BFMap*> SequenceType;
    struct MapList : public SequenceType {
      // return the first matching map or 0
      const BFMap *findMap(const CLHEP::Hep3Vector& x) const {
        for(SequenceType::const_iterator i = begin(); i!=end(); ++i) {
          if( (*i)->isValid(x)) {
            return *i;
          }
        }
        return 0;
      }
    };

    // An instance per (any) map, allows to optimize the lookup order of "inner" maps
    struct CacheElement {
      const BFMap *myMap; // the map this instance is attached to
      // A list of "inner" maps optimized for the "my" map
      // If "my" map is an inner map, it is not in the list.
      MapList inner;

      CacheElement(const BFMap *my, const MapList& in) : myMap(my), inner(in) {}
    };

    // Inner map lists optimized for the last used map
    mutable const CacheElement *innerForLastInner;
    mutable const CacheElement *innerForLastOuter; // never null

    // Outer maps in the user-specified order
    MapList outer;

    typedef std::map<const BFMap*, CacheElement> CacheType;
    CacheType innerCache; // keys are all inner maps
    CacheType outerCache; // keys are outer maps and 0
      
      int counter;
      
      
  public:

    BFCacheManager();

    void setMaps(const std::vector<BFMap>& innerMaps, const std::vector<BFMap>& outerMaps);
      
      
    // Returns pointers to an appropriate field map, or 0.
    const BFMap* findMap(const CLHEP::Hep3Vector& x) const {
        
        //std::cout << "CALLING findMap() " << std::endl;
        
      // First try to find if the point belong to any of the inner maps
      if(innerForLastInner) { // we were in an inner map last time

          //std::cout << "innerForLastInner is TRUE" << std::endl;
          
        if(innerForLastInner->myMap->isValid(x)) {
            
            //std::cout << "innerForLastInner->myMap->isValid(x)" << std::endl;
          // Cache update not needed, we are still in the same inner map            
          return innerForLastInner->myMap;
        }

        // The lookup order here is optimized
        const BFMap *newinner = innerForLastInner->inner.findMap(x);
          
        if(newinner) { // Update cache
            
            //std::cout << "found a newinner" << std::endl;
            //std::cout << "newinner key is " << newinner->getKey() << std::endl;

          CacheType::const_iterator p =  innerCache.find(newinner);
          assert(p != innerCache.end());
          innerForLastInner = &p->second;
          return newinner;
        }
          
      }
      else { // We were not in an inner map last time
          
          //std::cout << "innerForLastInner is NOT TRUE" << std::endl;

        // innerForLastOuter is never null
        const BFMap *newinner = innerForLastOuter->inner.findMap(x);
          
        if(newinner) { // Update cache
            
            //std::cout << "found a newinner" << std::endl;
            //std::cout << "newinner key is " << newinner->getKey() << std::endl;
            
          CacheType::const_iterator p =  innerCache.find(newinner);
          assert(p != innerCache.end());
          innerForLastInner = &p->second;
          return newinner;
        }
      }
        

        //std::cout << "got to pt A in findMap()" << std::endl;

      // The current point is not in any of the inner maps
      innerForLastInner = 0;
        
      // The lookup order of the outer maps is always the same
      const BFMap *newouter = outer.findMap(x);

      // Keep the inner map lookup optimized
      if(innerForLastOuter->myMap != newouter) {
        CacheType::const_iterator p =  outerCache.find(newouter);
        assert(p != outerCache.end());
        innerForLastOuter = &p->second;
      }

        //std::cout << "newouter key is " << newouter->getKey() << std::endl;
      return newouter;

    }

  };
}

#endif/*BFCacheManager_hh*/
