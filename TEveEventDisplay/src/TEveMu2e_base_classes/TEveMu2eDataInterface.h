#ifndef TEveMu2eDataInterface_h
#define TEveMu2eDataInterface_h
//ROOT
#include <TApplication.h>
#include <TEvePad.h>
#include <TObject.h>
#include <TSystem.h>
//libGeom
#include <TGeoManager.h>
//TEve
#include <TEveManager.h>
#include <TEveStraightLineSet.h>
//Mu2e General:
#include "Offline/GeometryService/inc/DetectorSystem.hh"
//TEveMu2e
#include "Offline/TEveEventDisplay/src/dict_classes/Collection_Filler.h"
#include "Offline/TEveEventDisplay/src/dict_classes/Geom_Interface.h"
#include "Offline/TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2e2DProjection.h"
#include "Offline/TEveEventDisplay/src/shape_classes/TEveMu2eCalorimeter.h"
#include "Offline/TEveEventDisplay/src/shape_classes/TEveMu2eTracker.h"
#include "Offline/TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eHit.h"
#include "Offline/TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eCRVEvent.h"
#include "Offline/TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eCluster.h"
#include "Offline/TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eCustomHelix.h"
#include "Offline/TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eStraightTrack.h"
#include "Offline/TEveEventDisplay/src/dict_classes/GeomUtils.h"
#include "Offline/CalorimeterGeom/inc/CaloGeomUtil.hh"
//C++
#include <limits>
#include <vector>
#include <tuple>
#include <algorithm>

namespace mu2e{
    class TEveMu2eDataInterface {
    public:
      #ifndef __CINT__
      TEveMu2eDataInterface() : fHitsList2DXY(0),fHitsList2DXZ(0),fHitsList3D(0),fTrkHitsList2DXY(0),fTrkHitsList2DXZ(0),fTrkHitsList3D(0),fTCHitsList2DXY(0),fTCHitsList2DXZ(0),fTCHitsList3D(0),
                                fCrystalHitList(0),fTrackList2DXY(0),fTrackList2DXZ(0),fTrackList3D(0), fClusterList2D_disk0(0), fClusterList2D_disk1(0), fClusterList3D(0), fCrvList2DXY(0),
                                fCrvList2DYZ(0),  fCrvList3D(0), fExtTrackList2D(0), fExtTrackList3D(0){};
      TEveMu2eDataInterface(const TEveMu2eDataInterface &);
      TEveMu2eDataInterface& operator=(const TEveMu2eDataInterface &);
      virtual ~TEveMu2eDataInterface(){};
      #endif
      bool show2D;
      TEveElementList *fHitsList2DXY;
      TEveElementList *fHitsList2DXZ;
      TEveElementList *fHitsList3D;
      TEveElementList *fTrkHitsList2DXY;
      TEveElementList *fTrkHitsList2DXZ;
      TEveElementList *fTrkHitsList3D;
      TEveElementList *fTCHitsList2DXY;
      TEveElementList *fTCHitsList2DXZ;
      TEveElementList *fTCHitsList3D;
      TEveElementList *fCrystalHitList;
      TEveElementList *fTrackList2DXY;
      TEveElementList *fTrackList2DXZ;
      TEveElementList *fTrackList3D;
      TEveElementList *fClusterList2D_disk0;
      TEveElementList *fClusterList2D_disk1;
      TEveElementList *fClusterList3D;
      TEveElementList *fCrvList2DXY;
      TEveElementList *fCrvList2DYZ;
      TEveElementList *fCrvList3D;
      TEveElementList *fExtTrackList2D;
      TEveElementList *fExtTrackList3D;

      std::vector<double> getTimeRange(bool firstloop, const ComboHitCollection *chcol, const CrvRecoPulseCollection *crvcoincol, const CaloClusterCollection *clustercol, const CaloHitCollection *cryHitcol, bool addCRV,  bool addHits, bool addCalo);

      void AddCRVInfo(bool firstloop, const CrvRecoPulseCollection *crvcoincol, double min_time, double max_time, TEveMu2e2DProjection *CRV2Dproj, bool Redraw, bool accumulate, TEveProjectionManager *TXYMgr,TEveProjectionManager *TRZMgr, TEveScene *scene1, TEveScene *scene2);
      std::vector<double> AddComboHits(bool firstloop, const ComboHitCollection *chcol, TEveMu2e2DProjection *tracker2Dproj, bool Redraw, double min_energy, double max_energy, double min_time, double max_time, bool accumulate, TEveProjectionManager *TXYMgr, TEveProjectionManager *TRZMgr, TEveScene *scene1, TEveScene *scene2);

      void AddTrkHits(bool firstloop, const ComboHitCollection *chcol, std::tuple<std::vector<std::string>, std::vector<const KalSeedCollection*>> track_tuple, TEveMu2e2DProjection *tracker2Dproj, bool Redraw, double min_energy, double max_energy, double min_time, double max_time, bool accumulate, TEveProjectionManager *TXYMgr, TEveProjectionManager *TRZMgr, TEveScene *scene1, TEveScene *scene2);
      void AddTimeClusters(bool firstloop, const TimeClusterCollection *tccol, TEveMu2e2DProjection *tracker2Dproj, bool Redraw, bool accumulate, TEveProjectionManager *TXYMgr, TEveProjectionManager *TRZMgr, TEveScene *scene1, TEveScene *scene2);
      template<class KTRAJ> void AddKinKalTrajectory( std::unique_ptr<KTRAJ> const& trajectory, TEveMu2eCustomHelix *line, TEveMu2eCustomHelix *line_twoDXY, TEveMu2eCustomHelix *line_twoDXZ);
      void FillKinKalTrajectory(bool firstloop, std::tuple<std::vector<std::string>, std::vector<const KalSeedCollection*>> track_tuple, TEveMu2e2DProjection *tracker2Dproj, TEveProjectionManager *TXYMgr, TEveProjectionManager *TRZMgr, TEveScene *scene1, TEveScene *scene2);

      std::vector<double> AddCaloClusters(bool firstloop, const CaloClusterCollection *clustercol,TEveMu2e2DProjection *calo2Dproj,  bool Redraw, double min_energy, double max_energy, double min_time, double max_time, bool accumulate, TEveProjectionManager *CfXYMgr, TEveProjectionManager *CfRZMgr, TEveScene *scene1, TEveScene *scene2);

      void AddCrystalHits(bool firstloop, const CaloHitCollection *cryHitcol, TEveMu2e2DProjection *calo2Dproj,  double min_time, double max_time, bool Redraw, bool accumulate, TEveProjectionManager *CfXYMgr, TEveProjectionManager *CfRZMgr, TEveScene *scene1, TEveScene *scene2);

      void AddCosmicTrack(bool firstloop, const CosmicTrackSeedCollection *cosmiccol, TEveMu2e2DProjection *tracker2Dproj, double min_time, double max_time, bool Redraw, bool accumulate, TEveProjectionManager *TXYMgr, TEveProjectionManager *TRZMgr, TEveScene *scene1, TEveScene *scene2);

      void AddHelixPieceWise3D(bool firstloop, std::tuple<std::vector<std::string>, std::vector<const KalSeedCollection*>> track_tuple, TEveMu2e2DProjection *tracker2Dproj, double min_time, double max_time, bool Redraw, bool accumulate, TEveProjectionManager *TXYMgr, TEveProjectionManager *TRZMgr, TEveScene *scene1, TEveScene *scene2);

      void AddHelixPieceWise2D(bool firstloop, std::tuple<std::vector<std::string>, std::vector<const KalSeedCollection*>> track_tuple, TEveMu2e2DProjection *tracker2Dproj, double min_time, double max_time, bool Redraw, bool accumulate, TEveProjectionManager *TXYMgr, TEveProjectionManager *TRZMgr, TEveScene *scene1, TEveScene *scene2);

      ClassDef(TEveMu2eDataInterface,0);

  }; //end class def

}//end namespace mu2e

#endif /*TEveMu2eDataInterface.h*/
