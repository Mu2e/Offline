//
// Class to describe flag bits used for straw hits
//
// $Id: StrawHitFlag.cc,v 1.4 2013/04/04 01:08:20 brownd Exp $
// $Author: brownd $
// $Date: 2013/04/04 01:08:20 $
//
// Original author David Brown
//
// Mu2e includes
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include <stdexcept>
#include <iostream>
#include <stdio.h>

namespace mu2e {
  unsigned StrawHitFlagDetail::_maxTrkId(7);

  std::string const& StrawHitFlagDetail::typeName() {
    static std::string type("StrawHitFlag");
    return type;
  }

  std::map<std::string,StrawHitFlagDetail::mask_type> const& StrawHitFlagDetail::bitNames() {
    static std::map<std::string,mask_type> bitnames;
    if(bitnames.size()==0){
      bitnames[std::string("Stereo")]               = bit_to_mask(stereo);
      bitnames[std::string("PanelCombo")]           = bit_to_mask(panelcombo);
      bitnames[std::string("ResolvedPhi")]          = bit_to_mask(resolvedphi);
      bitnames[std::string("TimeDivision")]         = bit_to_mask(tdiv);
      bitnames[std::string("EnergySelection")]      = bit_to_mask(energysel);
      bitnames[std::string("RadiusSelection")]      = bit_to_mask(radsel);
      bitnames[std::string("TimeSelection")]	    = bit_to_mask(timesel);
      bitnames[std::string("BackgroundCluster")]    = bit_to_mask(bkgclust);
      bitnames[std::string("Background")]           = bit_to_mask(bkg);
      bitnames[std::string("Isolated")]             = bit_to_mask(isolated);
      bitnames[std::string("Outlier")]              = bit_to_mask(outlier);
      bitnames[std::string("OtherBackground")]      = bit_to_mask(other);
      bitnames[std::string("TimeCluster")]	    = bit_to_mask(tclust);
      bitnames[std::string("CalorimeterPreSelection")] = bit_to_mask(calopresel);
      bitnames[std::string("CalorimeterSelection")] = bit_to_mask(calosel);
      bitnames[std::string("StrawXTalk")]	    = bit_to_mask(strawxtalk);
      bitnames[std::string("TrackerSelection")]     = bit_to_mask(trksel);
      bitnames[std::string("ElectronicsXTalk")]	    = bit_to_mask(elecxtalk);
      bitnames[std::string("Active")]		    = bit_to_mask(active);
      bitnames[std::string("GoodDOCA")]		    = bit_to_mask(doca);
      for(unsigned itrk=0;itrk<=_maxTrkId;++itrk){
	bitnames[trackBitName(itrk)] = bit_to_mask(trackBit(itrk));
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
