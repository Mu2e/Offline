// ITracker wire description
//
// $Id: ITrackerWiredata.cc,v 1.4 2012/09/25 10:08:28 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/09/25 10:08:28 $
//
// Original author G. Tassielli
//

#include "ITrackerGeom/inc/ITrackerWiredata.hh"

ITrackerWiredata::ITrackerWiredata(){
        PosMatrix = new TObjArray(0);
        NcelLayer = 0;
        epsilon = 0x0;
        alfa = 0x0;
        radius_z0 = 0x0;
}

ITrackerWiredata::~ITrackerWiredata(){
        PosMatrix->Delete();

        NcelLayer = 0;
        delete [] epsilon;
        delete [] alfa;
        delete [] radius_z0;
}
