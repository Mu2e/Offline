#ifndef TEveMu2eCluster_h
#define TEveMu2eCluster_h

#include <TObject.h>
#include <TEvePointSet.h>
#include <TEveCalo.h>
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eHit.h"
#include "Offline/TEveEventDisplay/src/dict_classes/GeomUtils.h"
#include <string.h>
#include <string>


namespace mu2e {
  class TEveMu2eCluster: public TEvePointSet{

    public:
      #ifndef __CINT__
      explicit TEveMu2eCluster();
      TEveMu2eCluster(CaloCluster cluster) : fCaloCluster_(cluster){};
      virtual ~TEveMu2eCluster(){};

      CaloCluster fCaloCluster_;

      void DrawCluster(const std::string &pstr, CLHEP::Hep3Vector COG, int energylevel, TEveElementList *list, std::vector<CLHEP::Hep3Vector> hits, bool addHits);
      void DrawCrystalHits(const std::string &pstr, CLHEP::Hep3Vector COG, TEveElementList *list);
      const  CLHEP::Hep3Vector GetPositon() { return fCaloCluster_.cog3Vector() ;}
      double GetEnergy() { return fCaloCluster_.energyDep();}
      std::string DataTitle(const std::string &pstr, double edep);
      #endif
      ClassDef(TEveMu2eCluster, 0);
  };
}
#endif
