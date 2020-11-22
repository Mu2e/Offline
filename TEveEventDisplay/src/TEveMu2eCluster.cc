#include "TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eCluster.h"
#include "TEveEventDisplay/src/dict_classes/GeomUtils.h"
using namespace mu2e;
namespace mu2e{

  TEveMu2eCluster::TEveMu2eCluster(){}
  
  std::string TEveMu2eCluster::DataTitle(const std::string &pstr, double edep){
        std::string dstr= "\nLayer: ";
        std::string strlst=pstr+dstr+std::to_string(edep);
        return(strlst);
  }

  void TEveMu2eCluster::DrawCluster(const std::string &pstr,  CLHEP::Hep3Vector cog, int energylevel, TEveElementList *ClusterList)
  {
    double edep = fCaloCluster.energyDep();
    this->SetTitle((DataTitle(pstr, edep)).c_str());
    hep3vectorTocm(cog);
    /*bool addHits = true;
    if(addHits){
       Calorimeter const &cal = *(GeomHandle<Calorimeter>());
       TEvePointSet *teve_hit2D = new TEvePointSet();
       for(unsigned h =0 ; h < fCaloCluster.caloCrystalHitsPtrVector().size();h++){
        art::Ptr<CaloCrystalHit>  crystalhit = fCaloCluster.caloCrystalHitsPtrVector()[h] ;
        int diskId = cal.crystal(crystalhit->id()).diskId();
        CLHEP::Hep3Vector HitPos(cal.geomUtil().mu2eToDiskFF(diskId, cal.crystal(crystalhit->id()).position()));
        CLHEP::Hep3Vector pointInMu2e = PointToCalo(HitPos,diskId);
        hep3vectorTocm(pointInMu2e);
        teve_hit2D->SetMarkerColor(kGreen);
        teve_hit2D->SetNextPoint(pointInMu2e.x(), pointInMu2e.y(), pointInMu2e.z()); 
        teve_hit2D->SetMarkerSize(2);
        teve_hit2D->SetPickable(kTRUE);
        ClusterList->AddElement(teve_hit2D);
      }
    }*/
    Int_t mSize = 3;
    int colors[] = {+10, +5, +7, +8, -3, +1, -5, 0, -2, -4, +6, -9};
    this->SetMarkerColor(kViolet + colors[energylevel]);
    this->SetNextPoint(cog.x(), cog.y(), cog.z()); 
    this->SetMarkerStyle(9);
    this->SetMarkerSize(mSize);
    this->SetPickable(kTRUE);
    ClusterList->AddElement(this);
  }
  
  void TEveMu2eCluster::DrawCrystalHits(const std::string &pstr, CLHEP::Hep3Vector cog, TEveElementList *ClusterList){
    hep3vectorTocm(cog);
    Int_t mSize = 2;
    this->SetMarkerColor(kGreen);
    this->SetNextPoint(cog.x(), cog.y(), cog.z()); 
    this->SetMarkerStyle(31);
    this->SetMarkerSize(mSize);
    this->SetPickable(kTRUE);
    ClusterList->AddElement(this);
  }
}

