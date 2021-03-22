#ifndef TEveMu2eCluster_h
#define TEveMu2eCluster_h

#include <TObject.h>
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eHit.h"
#include <string.h>
#include <string>
#include <TEvePointSet.h>
#include <TEveCalo.h>

namespace mu2e {
  class TEveMu2eCluster: public TEvePointSet{
      CaloCluster fCaloCluster;   
      
    public:
      #ifndef __CINT__
      explicit TEveMu2eCluster();
      TEveMu2eCluster(CaloCluster cluster) : fCaloCluster(cluster){};
      virtual ~TEveMu2eCluster(){};
      #endif 
      void DrawCluster(const std::string &pstr, CLHEP::Hep3Vector COG, int energylevel, TEveElementList *list, std::vector<CLHEP::Hep3Vector> hits, bool addHits);
      void DrawCrystalHits(const std::string &pstr, CLHEP::Hep3Vector COG, TEveElementList *list);
      const  CLHEP::Hep3Vector GetPositon() { return fCaloCluster.cog3Vector() ;}
      double GetEnergy() { return fCaloCluster.energyDep();}
      std::string DataTitle(const std::string &pstr, double edep);
      ClassDef(TEveMu2eCluster, 0);
  };
}
#endif
