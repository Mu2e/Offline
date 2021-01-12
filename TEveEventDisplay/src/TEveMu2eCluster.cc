#include "TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eCluster.h"
#include "TEveEventDisplay/src/dict_classes/GeomUtils.h"
using namespace mu2e;
namespace mu2e{

  TEveMu2eCluster::TEveMu2eCluster(){}

  void TEveMu2eCluster::DrawCluster(const std::string &pstr,  CLHEP::Hep3Vector cog, int energylevel, TEveElementList *ClusterList)
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
    ClusterList->AddElement(this);
  }
}

