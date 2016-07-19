#ifndef RecoDataProducts_StrawHitFlag_hh
#define RecoDataProducts_StrawHitFlag_hh
//
// Class to describe flag bits used for straw hits
// 
// $Id: StrawHitFlag.hh,v 1.4 2013/04/04 01:08:20 brownd Exp $
// $Author: brownd $
// $Date: 2013/04/04 01:08:20 $
//
// Original author David Brown
//
// Mu2e includes
#include "GeneralUtilities/inc/BitMap.hh"
#include <string>
#include <map>
namespace mu2e {

  struct StrawHitFlagDetail {
// I need 32 bits for this class
    typedef unsigned mask_type;
// the lower 16 bits are used to flag hit properties.  bits 0-10 refer to track-related
// properties, 12-15 are for external association (like to a calorimeter cluster)
// The upper 16 bits are reserved to flag the track number (or cluster) to which these hits are associated.
    enum bit_type {stereo=0, energysel=1, radsel=2, timesel=3,  delta=6, isolated=7, outlier=8, other=9,
    tdiv=10,trksel=11,
    calosel=12, strawxtalk=13, elecxtalk=14,
    track0=16,track1=17,track2=18,track3=19,track4=20,track5=21,track6=22,track7=23,
    track8=24,track9=25,track10=26,track11=27,track12=28,track13=29,track14=30,track15=31};
// special function to return the enum value associated with a given track number
    static bit_type trackBit(unsigned itrk);
    static std::string trackBitName(unsigned itrk);
    // functions needed for the BitMap template
    static std::string const& typeName();
    static std::map<std::string,mask_type> const& bitNames();
// maximum track Id I can flag
    static unsigned _maxTrkId;
    static mask_type bit_to_mask( bit_type b){ return 1<<b; }
  };
  typedef BitMap<StrawHitFlagDetail> StrawHitFlag;

}
#endif

