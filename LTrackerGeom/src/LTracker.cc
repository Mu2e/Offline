//
// Geometry and identifier info about an LTracker.
//
//
// $Id: LTracker.cc,v 1.2 2009/10/22 16:27:58 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/10/22 16:27:58 $
//
// Original author Rob Kutschke
//

#include "LTrackerGeom/inc/LTracker.hh"

using namespace std;

namespace mu2e { 

  double LTracker::zHalfLength() const{
    
    double zHalfMax = 0.;
    /*
    for ( std::vector<StrawDetail>::const_iterator i = _strawDetail.begin(),
	    e=_strawDetail.end(); 
	  i != e ; 
	  ++i ){
      double l = i->halfLength();
      zHalfMax = ( l > zHalfMax ) ? l : zHalfMax;
    }
    */
    for ( int i=0; i<_strawDetail.size(); ++i ) {
      double l = _strawDetail[i].halfLength();
      zHalfMax = ( l > zHalfMax ) ? l : zHalfMax;
    }
    return zHalfMax;
  }

// Given the _indices member in each Layer, fill the _straws member.
void LTracker::FillPointers1(){

  for ( vector<Device>::iterator idev = _devices.begin(),
	  edev = _devices.end(); 
        idev != edev;  ++idev ){
  
    for ( vector<Sector>::iterator isec = idev->_sectors.begin(),
            esec = idev->_sectors.end();
	  isec != esec; ++isec ){

      for ( vector<Layer>::iterator ilay = isec->_layers.begin(),
	      elay = isec->_layers.end();
	    ilay != elay; ++ilay ){
	ilay->_straws.clear();
	for ( vector<StrawIndex>::iterator istr = ilay->_indices.begin(),
		estr = ilay->_indices.end();
	      istr != estr ; ++istr ){
	  const Straw& str = _allStraws[(*istr).asInt()];
	  ilay->_straws.push_back( &str );
	}
      }
    }
  }

}

void LTracker::FillPointers2(){

  // Fill nearest neighbour indices and pointers from the NN Ids.
  for ( deque<Straw>::iterator i= _allStraws.begin(), 
	  e= _allStraws.end();
	i!=e; 
	++i){
    vector<const Straw *>& byPtr = i->_nearest;
    vector<StrawId>& byId        = i->_nearestById;
    vector<StrawIndex>& byIndex  = i->_nearestByIndex;
    
    byPtr.clear();
    byIndex.clear();
    
    for ( vector<StrawId>::iterator j=byId.begin(), je=byId.end();
	   j != je; ++j){
      const StrawId& id = *j;
      const Straw& straw = getStraw(id);
      byPtr.push_back( &straw);
      byIndex.push_back( straw.Index() );
    }
  }

}

} // namespace mu2e
