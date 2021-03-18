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

  void TEveMu2eCluster::DrawCluster(const std::string &pstr,  CLHEP::Hep3Vector cog, int energylevel, TEveElementList *ClusterList,  std::vector<CLHEP::Hep3Vector> hits, bool addHits)
  {
    double edep = fCaloCluster.energyDep();
    this->SetTitle((DataTitle(pstr, edep)).c_str());
    hep3vectorTocm(cog);
    Int_t mSize = 3;
    int colors[] = {+10, +5, +7, +8, -3, +1, -5, 0, -2, -4, +6, -9};
    this->SetMarkerColor(kViolet + colors[energylevel]);
    this->SetNextPoint(cog.x(), cog.y(), cog.z()); 
    this->SetMarkerStyle(9);
    this->SetMarkerSize(mSize);
    this->SetPickable(kTRUE);

    if(addHits){
      
       TEvePointSet *teve_hit2D = new TEvePointSet();
       for(unsigned int h =0 ; h < hits.size();h++) {
        teve_hit2D->SetNextPoint(hits[h].x(), hits[h].y(), hits[h].z());
        teve_hit2D->SetMarkerSize(2);
        teve_hit2D->SetMarkerStyle(47);
        teve_hit2D->SetMarkerColor(kViolet + colors[energylevel]);
        ClusterList->AddElement(teve_hit2D);
      }
    }
    
   
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

