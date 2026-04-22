 //
// An EDAnalyzer module that reads back the hits created by the calorimeter and produces an ntuple
//
//
// Original author Bertrand Echenard
//

#include "CLHEP/Units/SystemOfUnits.h"

#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/VirtualDetector.hh"


#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/PtrStepPointMCVector.hh"
#include "Offline/MCDataProducts/inc/GenId.hh"
#include "Offline/DataProducts/inc/VirtualDetectorId.hh"

#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Provenance.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/InputTag.h"

#include "TDirectory.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1F.h"

#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <vector>


namespace mu2e {


  class CaloClusterCheck : public art::EDAnalyzer {

     public:
       explicit CaloClusterCheck(fhicl::ParameterSet const& pset);

       virtual void beginJob();
       virtual void analyze(const art::Event& e);


     private:
       std::string caloCrystalModuleLabel_;
       std::string caloClusterModuleLabel_;
       int diagLevel_;
       int nProcess_;

       TH2F *crIn_,*crOut_;
       TH1F *dtime_;
  };



  CaloClusterCheck::CaloClusterCheck(fhicl::ParameterSet const& pset) :
     art::EDAnalyzer(pset),
     caloCrystalModuleLabel_(pset.get<std::string>("caloCrystalModuleLabel")),
     caloClusterModuleLabel_(pset.get<std::string>("caloClusterModuleLabel")),
     diagLevel_(pset.get<int>("diagLevel",0)),
     nProcess_(0),
     crIn_(0),crOut_(0),dtime_(0)
  {}



  void CaloClusterCheck::beginJob()
  {
     art::ServiceHandle<art::TFileService> tfs;
     crIn_  = tfs->make<TH2F>("crIn","dist diff vs time diff for crystal in",100,0.,200.,100,0,100);
     crOut_ = tfs->make<TH2F>("crOut","dist diff vs time diff for crystal out",100,0.,200.,100,0,100);
     dtime_ = tfs->make<TH1F>("dtime","delta time in crystal",100,0.,100);
  }



  void CaloClusterCheck::analyze(const art::Event& event)
  {
      ++nProcess_;
      if (nProcess_%10==0 && diagLevel_ > 0) std::cout<<"Processing event from CaloClusterCheck =  "<<nProcess_ << std::endl;


      //Get handle to the calorimeter
      art::ServiceHandle<GeometryService> geom;
      if( ! geom->hasElement<Calorimeter>() ) return;
      Calorimeter const & cal = *(GeomHandle<Calorimeter>());

      //Get calo crystal hits (average from readouts)
      art::Handle<CaloHitCollection> CaloHitsHandle;
      event.getByLabel(caloCrystalModuleLabel_, CaloHitsHandle);
      CaloHitCollection const& CaloHits(*CaloHitsHandle);

      //Get calo cluster
      art::Handle<CaloClusterCollection> caloClustersHandle;
      event.getByLabel(caloClusterModuleLabel_, caloClustersHandle);
      CaloClusterCollection const& caloClusters(*caloClustersHandle);


      for (const auto& hit : CaloHits)
      {
          if (hit.time() < 500 || hit.energyDep() < 0.9) continue;

          int nClu(0),nCr(0);
          std::vector<double> dmin,dtime;
          for (auto const& clusterIt : caloClusters)
          {
               for (auto const& crystalIncluster : clusterIt.caloHitsPtrVector())
               {
                    if (cal.crystal(hit.crystalID()).diskID() != cal.crystal(crystalIncluster->crystalID()).diskID()) break;

                    if (&(*crystalIncluster) == &hit) ++nCr;
                    if (nCr > 1) std::cout<<"Warning, crystal associated to more than one cluster "<<hit.crystalID()<<" for cluster "<<nClu<<std::endl;

                    CLHEP::Hep3Vector crystalPos2 = cal.crystal(crystalIncluster->crystalID()).position();
                    double deltaDist = (cal.crystal(hit.crystalID()).position()-crystalPos2).mag();
                    double deltaTime = abs(crystalIncluster->time() - hit.time());

                    dmin.push_back(deltaDist);
                    dtime.push_back(deltaTime);
               }
               ++nClu;
           }

           if (nCr!=0) continue;
           for (unsigned int i=0;i<dmin.size();++i) crOut_->Fill(dmin[i],dtime[i]);
       }


       for (auto const& clusterIt : caloClusters)
       {
           double dtime(0);
           for (int i=0;i<clusterIt.size()-1;++i)
           {
              for (int j=i+1;j<clusterIt.size();++j)
              {
                  double deltaTime = abs(clusterIt.caloHitsPtrVector().at(i)->time()- clusterIt.caloHitsPtrVector().at(j)->time());
                  if (deltaTime > dtime) dtime = deltaTime;
              }
           }
           dtime_->Fill(dtime);
       }

  }

}

DEFINE_ART_MODULE(mu2e::CaloClusterCheck)


