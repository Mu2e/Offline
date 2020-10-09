#ifndef TEveMu2eCluster_h
#define TEveMu2eCluster_h

#include <TObject.h>
#include "RecoDataProducts/inc/CaloCluster.hh"
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
      void DrawCluster(const std::string &pstr, CLHEP::Hep3Vector COG, int energylevel, TEveElementList *list); 
      const  CLHEP::Hep3Vector GetPositon() { return fCaloCluster.cog3Vector() ;}
      double GetEnergy() { return fCaloCluster.energyDep();}
      inline std::string DataTitle(const std::string &pstr, double edep){
        std::string dstr= "\nLayer: ";
        std::string strlst=pstr+dstr+std::to_string(edep);
        return(strlst);
      }
      ClassDef(TEveMu2eCluster, 0);
  };
}
#endif
