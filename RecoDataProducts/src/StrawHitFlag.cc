//
// Class to describe flag bits used for straw hits
// 
// $Id: StrawHitFlag.cc,v 1.1 2013/03/02 20:48:18 brownd Exp $
// $Author: brownd $
// $Date: 2013/03/02 20:48:18 $
//
// Original author David Brown
//
// Mu2e includes
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include <stdexcept>
#include <iostream>
#include <stdio.h>

namespace mu2e {
  unsigned StrawHitFlagDetail::_maxTrkId(15);

  std::string const& StrawHitFlagDetail::typeName() {
    static std::string type("StrawHitFlag");
    return type;
  }

  std::map<StrawHitFlagDetail::mask_type,std::string> const& StrawHitFlagDetail::bitNames() {
    static std::map<mask_type,std::string> bitnames;
    if(bitnames.size()==0){
      bitnames[1<<stereo] = std::string("Stereo");
      bitnames[1<<energysel] = std::string("EnergySelection");
      bitnames[1<<radsel] = std::string("RadiusSelection");
      bitnames[1<<delta] = std::string("DeltaRay");
      bitnames[1<<isolated] = std::string("Isolated");
      bitnames[1<<outlier] = std::string("Outlier");
      bitnames[1<<calosel] = std::string("CalorimeterSelection");
      for(unsigned itrk=0;itrk<=_maxTrkId;++itrk){
	bitnames[1<<trackBit(itrk)] = trackBitName(itrk);
      }
    }
    return bitnames;
  }

  StrawHitFlagDetail::bit_type StrawHitFlagDetail::trackBit(unsigned itrk) {
    if(itrk > _maxTrkId ){
      std::ostringstream os;
      os << typeName() << " invalid track number : " << itrk;
      throw std::out_of_range( os.str() );
    }
    return static_cast<bit_type>(16+itrk);
  }

  std::string StrawHitFlagDetail::trackBitName(unsigned itrk) {
    if(itrk > _maxTrkId ){
      std::ostringstream os;
      os << typeName() << " invalid track number : " << itrk;
      throw std::out_of_range( os.str() );
    }
    static char bitname[20];
    snprintf(bitname,20,"Track%i",itrk);
    return std::string(bitname);
  }

}
