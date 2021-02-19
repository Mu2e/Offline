#ifndef TEveMu2eProjectionInterface_h
#define TEveMu2eProjectionInterface_h

#include "TEveEventDisplay/src/dict_classes/Geom_Interface.h"
#include "TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2e2DProjection.h"
#include "TEveEventDisplay/src/shape_classes/TEveMu2eCalorimeter.h"
#include "TEveEventDisplay/src/shape_classes/TEveMu2eTracker.h"
#include "TEveEventDisplay/src/shape_classes/TEveMu2eCRV.h"

namespace mu2e{
	class TEveMu2eProjectionInterface : public TGMainFrame {
	public:
        #ifndef __CINT__	
	void CreateCRVProjection(TEveMu2e2DProjection *CRV2Dproj);
	void CreateCaloProjection(TEveMu2e2DProjection *calo2Dproj);
	void CreateTrackerProjection(TEveMu2e2DProjection *tracker2DProj);
	#endif
	ClassDef(TEveMu2eProjectionInterface,0);

    }; //end class def

}//end namespace mu2e

#endif /*TEveMu2eProjectionInterface.h*/
