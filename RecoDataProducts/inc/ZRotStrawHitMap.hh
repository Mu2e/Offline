#ifndef RecoDataProducts_ZRotStrawHitMap_hh
#define RecoDataProducts_ZRotStrawHitMap_hh
//
// index to access to the StrawHit by z ID position and by AbsSttn ID, it is the input for the final pattern recognition
// stage in the fast pattern recognition
//
// $Id: ZRotStrawHitMap.hh,v 1.1 2011/10/28 00:17:18 tassiell Exp $
// $Author: tassiell $
// $Date: 2011/10/28 00:17:18 $
//
// Original author G. Tassielli
//

// C++ includes
#include <map>
//#include <functional>
//#include <utility>
//#include <algorithm>

// Mu2e includes
#include "RecoDataProducts/inc/StrawHit.hh"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

typedef art::Ptr<mu2e::StrawHit> StrawHitPtr;
typedef std::multimap<unsigned int , mu2e::StrawHitPtr, std::less<unsigned int> > AbsSectStrawHitMap;
typedef std::map<unsigned int , mu2e::AbsSectStrawHitMap, std::less<unsigned int> > ZSectStrawHitMap;
typedef std::multimap<unsigned int , mu2e::StrawHitPtr, std::less<unsigned int> > ZStrawHitMap;
typedef std::map<unsigned int , mu2e::ZStrawHitMap, std::less<unsigned int> > SectZStrawHitMap;

  class ZRotStrawHitMap{

  public:

    ZSectStrawHitMap _zsctTrackerHits;
    SectZStrawHitMap _sctZTrackerHits;
    float  _min_Time;
    float  _max_Time;
    float  _min_HStep;
    float  _max_HStep;
    size_t _NHits;

    ZRotStrawHitMap():
            _min_Time(0.00000),
            _max_Time(0.00000),
            _min_HStep(0.00000),
            _max_HStep(0.00000),
            _NHits(0) {
            _zsctTrackerHits.clear();
    }

    ZRotStrawHitMap( float min_Time, float max_Time, float min_HStep, float max_HStep):
            _min_Time(min_Time),
            _max_Time(max_Time),
            _min_HStep(min_HStep),
            _max_HStep(max_HStep),
            _NHits(0) {
            _zsctTrackerHits.clear();
    }

    ~ZRotStrawHitMap(){};

    void AddHit( unsigned int absZpos, unsigned int absSect, const StrawHitPtr &ihit ){
            ZSectStrawHitMap::iterator tmpZSectMap_it = _zsctTrackerHits.find(absZpos);
            if ( tmpZSectMap_it!=_zsctTrackerHits.end() ) {
                    tmpZSectMap_it->second.insert( AbsSectStrawHitMap::value_type(absSect,ihit) );
            }
            else {
                    AbsSectStrawHitMap tmpAbsSectMap;
                    tmpAbsSectMap.insert( AbsSectStrawHitMap::value_type(absSect,ihit) );
                    _zsctTrackerHits.insert( ZSectStrawHitMap::value_type(absZpos,tmpAbsSectMap) );
                    //tmpZSectMap_it = _zsctTrackerHits.insert( ZSectStrawHitMap::value_type( absZpos, AbsSectStrawHitMap() ) );
                    //tmpZSectMap_it->second.insert( AbsSectStrawHitMap::value_type(absSect,ihit) );
            }

            SectZStrawHitMap::iterator tmpSectZMap_it = _sctZTrackerHits.find(absSect);
            if ( tmpSectZMap_it!=_sctZTrackerHits.end() ) {
                    tmpSectZMap_it->second.insert( ZStrawHitMap::value_type(absZpos,ihit) );
            }
            else {
                    ZStrawHitMap tmpZMap;
                    tmpZMap.insert( ZStrawHitMap::value_type(absZpos,ihit) );
                    _sctZTrackerHits.insert( SectZStrawHitMap::value_type(absSect,tmpZMap) );
                    //tmpSectZMap_it = _sctZTrackerHits.insert( SectZStrawHitMap::value_type( absSect, ZStrawHitMap() ) );
                    //tmpSectZMap_it->second.insert( ZStrawHitMap::value_type(absZpos,ihit) );
            }

            ++_NHits;
    }

    size_t NHit() {
            _NHits=0;
            for ( ZSectStrawHitMap::iterator zsctTrackerHits_it=_zsctTrackerHits.begin(); zsctTrackerHits_it!=_zsctTrackerHits.end(); ++zsctTrackerHits_it ){
                    _NHits+=zsctTrackerHits_it->second.size();
            }
            return _NHits;
    }

    // Print contents of the object.
    //void print( std::ostream& ost = std::cout, bool doEndl = true ) const;
    friend std::ostream& operator<< ( std::ostream& ost,
                                      ZRotStrawHitMap const& hit){
      ost<<"Z-Rot Straw Hit map: "<<std::endl;
      ost<<"\t Time in range: "<<hit._min_Time<<" - "<<hit._max_Time<<" , expected Helix Step in range: "<<hit._min_HStep<<" - "<<hit._max_HStep<<std::endl;
      ost<<"\t number of hits: "<<hit._NHits<<std::endl;
      return ost;
    }
  };

//  inline std::ostream& operator<<( std::ostream& ost,
//                                   ZRotStrawHitMap const& hit){
//    ost<<"Z-Rot Straw Hit map: "<<std::endl;
//    ost<<"\t Time in range: "<<hit._min_Time<<" - "<<hit._max_Time<<" , expected Helix Step in range: "<<hit._min_HStep<<" - "<<hit._max_HStep<<std::endl;
//    ost<<"\t number of hits: "<<hit._NHits<<std::endl;
//    return ost;
//  }

} // namespace mu2e

#endif /* RecoDataProducts_ZRotStrawHitMap_hh */
