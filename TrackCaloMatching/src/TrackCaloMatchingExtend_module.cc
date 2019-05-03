//
// Original author B. Echenard
//
// A generic track calo matching algorithm. 
//
// The matching is performed in the calorimeter front face section frame, hence only x and y coordinate matter (z is always zero in that frame)
//
// As of today, use ad-hoc errors, will work on including the event-by-event error later
// As of today, use x-y coordinates, might want to shift to u-v
//
//


#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"


#include "BTrk/BbrGeom/HepPoint.h"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "RecoDataProducts/inc/TrkCaloIntersectCollection.hh"
#include "RecoDataProducts/inc/TrkCaloMatchCollection.hh"

#include "BTrk/TrkBase/HelixParams.hh"
#include "BTrk/TrkBase/TrkRep.hh"
#include "BTrk/TrkBase/HelixTraj.hh"

#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <functional>
#include <numeric>

#include "TH1D.h"
#include "TMVA/Reader.h"





namespace mu2e {

   class TrackCaloMatchingExtend : public art::EDProducer {

       public:

	   explicit TrackCaloMatchingExtend(fhicl::ParameterSet const& pset):
             art::EDProducer{pset},
	     caloClusterModuleLabel_    (pset.get<std::string>("caloClusterModuleLabel")),
	     trkIntersectModuleLabel_   (pset.get<std::string>("trkIntersectModuleLabel")),
	     trackExtend_               (pset.get<double>("trackExtend")),  
	     dtOffset_                  (pset.get<double>("dtOffset")),  
	     sigmaT_                    (pset.get<double>("sigmaT")),  
	     sigmaXY_                   (pset.get<double>("sigmaXY")),  
	     chi2Cut_                   (pset.get<double>("chi2Cut")),  
	     diagLevel_                 (pset.get<int>   ("diagLevel",0))
	   {
               produces<TrkCaloMatchCollection>();
	   }

	   virtual ~TrackCaloMatchingExtend() {}

	   void beginJob();
	   void endJob  () {}
	   void produce (art::Event& e);



       private:

           void matchMe(TrkCaloMatchCollection& trackClusterMatch, art::Handle<TrkCaloIntersectCollection>  trjIntersectHandle, 
                        art::Handle<CaloClusterCollection> caloClustersHandle);

           double calcChi2Pos(const Calorimeter& cal, const CaloCluster& cluster, 
                              const CLHEP::Hep3Vector& trkPos, double trkTime, double cellsize);
                                            
           CLHEP::Hep3Vector recalculateCog(const Calorimeter& cal, const CaloCluster& cluster);
          

	   std::string      caloClusterModuleLabel_;
	   std::string      trkIntersectModuleLabel_;
           double           trackExtend_;
           double           dtOffset_;
  	   double           sigmaT_;
	   double           sigmaXY_;
	   double           chi2Cut_;
          
	   int              diagLevel_;
           
	   TH1F *_hXY,*_hT,*_hChi2,*_hChi2Time,*_hChi2Pos;


   };


   void TrackCaloMatchingExtend::beginJob()
   {

       art::ServiceHandle<art::TFileService> tfs;
       _hXY       = tfs->make<TH1F>("XY_diff","Clu - TRK X-Y dist",100,-100,100.);
       _hT        = tfs->make<TH1F>("T_diff","Clu - TRK Time dist",100,-10,10.);
       _hChi2     = tfs->make<TH1F>("Chi2",    "Chi2",200,0,400.);
       _hChi2Time = tfs->make<TH1F>("Chi2Pos", "Chi2",200,0,400.);
       _hChi2Pos  = tfs->make<TH1F>("Chi2Time","Chi2",200,0,400.);
              
   }



   //-----------------------------------------------------------------------------
   void TrackCaloMatchingExtend::produce(art::Event& event)
   {
	 
         art::Handle<CaloClusterCollection> caloClustersHandle;
	 event.getByLabel(caloClusterModuleLabel_, caloClustersHandle);
	 const CaloClusterCollection& caloClusters(*caloClustersHandle);

         art::Handle<TrkCaloIntersectCollection>  trjIntersectHandle;
         event.getByLabel(trkIntersectModuleLabel_, trjIntersectHandle);
         const TrkCaloIntersectCollection& trkIntersects(*trjIntersectHandle);

	 //output of TrackCaloMatchingExtend
         std::unique_ptr<TrkCaloMatchCollection> trackClusterMatch(new TrkCaloMatchCollection);


	 if (diagLevel_) std::cout<<"Event Number: "<< event.event()<<"\nStart TrkCaloMatching with nExtrapol = "
                                  <<trkIntersects.size()<<", nCluster="<<caloClusters.size()<<std::endl;

	 if (trkIntersects.size() && caloClusters.size()) matchMe(*trackClusterMatch, trjIntersectHandle, caloClustersHandle);

	 event.put(std::move(trackClusterMatch));
   }



   //-----------------------------------------------------------------------------
   void TrackCaloMatchingExtend::matchMe(TrkCaloMatchCollection& trackClusterMatch, art::Handle<TrkCaloIntersectCollection>  trjIntersectHandle, 
                                      art::Handle<CaloClusterCollection> caloClustersHandle)
   {

       const Calorimeter& cal = *(GeomHandle<Calorimeter>());

       const TrkCaloIntersectCollection& trkIntersects(*trjIntersectHandle);
       const CaloClusterCollection& caloClusters(*caloClustersHandle);

       
       double cellsize = cal.caloInfo().getDouble("crystalXYLength")+2.0*cal.caloInfo().getDouble("wrapperThickness");
       const auto* trkIntersectBase = &trkIntersects.front();
       const auto* caloClusterBase = &caloClusters.front();

       
       if (diagLevel_ > 1)
	   for (const auto& cluster: caloClusters) std::cout<<"In TrackClusterMatching, cluster energy="
                                                             <<cluster.energyDep()<<"  time="<<cluster.time()
                                                             <<"  cog="<<cluster.cog3Vector()<<std::endl;


       std::vector<double> chi2vec,chi2Timevec,chi2Posvec;
       std::vector<int>    itrkvec,icluvec;
       
       for (const auto& trkIntersect : trkIntersects)
       {
           double pathLength             = trkIntersect.pathLengthEntrance();	       
	   double trkTime                = trkIntersect.trk()->arrivalTime(pathLength);
	   CLHEP::Hep3Vector trkMomentum = trkIntersect.trk()->momentum(pathLength);
	   HepPoint          point       = trkIntersect.trk()->position(pathLength);
                      	            
           CLHEP::Hep3Vector posTrkInTracker(point.x(),point.y(),point.z());	     
	   CLHEP::Hep3Vector posTrkInSectionFF = cal.geomUtil().mu2eToDiskFF(trkIntersect.diskId(),cal.geomUtil().trackerToMu2e(posTrkInTracker));


           //needed to compute the trajectory inside the disk
           HelixTraj trkHel(trkIntersect.trk()->helix(pathLength).params(),trkIntersect.trk()->helix(pathLength).covariance());

           double phi0          = trkHel.phi0();
           double omega         = trkHel.omega();
           double radius        = 1.0/omega;
           double cosDip        = trkHel.cosDip();
           double sinDip        = sqrt(1-cosDip*cosDip);

           double centerCircleX = (trkHel.d0() + radius)*sin(phi0);
           double centerCircleY = (trkHel.d0() + radius)*cos(phi0);	      
           double length = (posTrkInTracker.z()-trkHel.z0())/sinDip + trackExtend_;
           double x =  radius*sin(phi0+omega*cosDip*length) - centerCircleX;
           double y = -radius*cos(phi0+omega*cosDip*length) + centerCircleY;
           double z = posTrkInSectionFF.z()+trackExtend_/sinDip;

           CLHEP::Hep3Vector posTrkMatch(x,y,z);
           

 	   for (const auto& cluster : caloClusters)
           {
               
               if (trkIntersect.diskId() != cluster.diskId()) continue;               

               CLHEP::Hep3Vector diff = cluster.cog3Vector()-posTrkInSectionFF;
               double deltaTime       = std::abs(cluster.time()-trkTime-dtOffset_);
               
               if (deltaTime > 50) continue;
               //if (sqrt(diff.x()*diff.x()+diff.y()*diff.y()) > 200) continue;

               double chi2Time = pow(deltaTime/sigmaT_,2);                                            
               double chi2Pos  = calcChi2Pos(cal, cluster, posTrkMatch, trkTime, cellsize);
               double chi2     = chi2Time + chi2Pos;
	       
	       if (chi2Time<100 && diagLevel_ > 2) {_hChi2->Fill(chi2);_hChi2Pos->Fill(chi2Pos);}
               if (chi2Pos<100  && diagLevel_ > 2) {_hT->Fill(cluster.time()-trkTime-dtOffset_);_hChi2Time->Fill(chi2Time);}
               
               
               if (chi2 > chi2Cut_) continue;
               
               size_t itrk = (&trkIntersect - trkIntersectBase);    	       
               size_t iclu = (&cluster - caloClusterBase);

               itrkvec.push_back(itrk);
               icluvec.push_back(iclu);
               chi2vec.push_back(chi2);
               chi2Timevec.push_back(chi2Time);
               chi2Posvec.push_back(chi2Pos);

               if (diagLevel_ > 0) std::cout <<"TrackCaloMatchingExtend inserted match for icl="<<iclu<<" trkMAtch="<<itrk
	                                     <<"  cluster pos= "<<cluster.cog3Vector()<<"   trk pos="<<posTrkInSectionFF
  	                                     <<"  with chi2="<<chi2<<" with Chi2Time= "<<chi2Time<<" and Chi2Pos = "<<chi2Pos<<std::endl;	     
           }                          	      
       }
       

       //At this point, save only the best match, but you can easily change to save them all
       
       if (chi2vec.empty()) return; 
       size_t idx = std::min_element(chi2vec.begin(),chi2vec.end())-chi2vec.begin();
       
       art::Ptr<TrkCaloIntersect> extrapolPtr = art::Ptr<TrkCaloIntersect>(trjIntersectHandle,itrkvec[idx]);
       art::Ptr<CaloCluster> clusterPtr = art::Ptr<CaloCluster>(caloClustersHandle,icluvec[idx]);	       

       TrkCaloMatch match(extrapolPtr,clusterPtr, icluvec[idx], chi2vec[idx], chi2Timevec[idx], chi2Posvec[idx]);
       trackClusterMatch.push_back(match);     
          
   }



   //---------------------------------------------------------------------------------------------------------
   double TrackCaloMatchingExtend::calcChi2Pos(const Calorimeter& cal, const CaloCluster& cluster, 
                                            const CLHEP::Hep3Vector& trkPos, double trkTime, double cellsize)
   {
        
        //need to recalculate the COG here
        CLHEP::Hep3Vector newCog = recalculateCog(cal,cluster);
                        
        double dx = newCog.x()-trkPos.x();
        double dy = newCog.y()-trkPos.y();
        
        double chi2 = (dx*dx + dy*dy)/sigmaXY_/sigmaXY_;
                                     
	if (diagLevel_ > 2) _hXY->Fill(dx);
	if (diagLevel_ > 2) _hXY->Fill(dy);
        
        return chi2;        
    }


    //---------------------------------------------------------------------------------------------------------
    CLHEP::Hep3Vector TrackCaloMatchingExtend::recalculateCog(const Calorimeter& cal, const CaloCluster& cluster)
    {
         auto const& main = cluster.caloCrystalHitsPtrVector();

         double sxi(0),syi(0),swi(0);
         for (auto it = main.begin(); it !=main.end(); ++it)
         {
              int    crId((*it)->id());
              double energy((*it)->energyDep());

              CLHEP::Hep3Vector crystalPos = cal.geomUtil().mu2eToDiskFF(cal.crystal(crId).diskId(), cal.crystal(crId).position());

              double weight = energy - 4.939;
              //double weight = -5.45 + 2.63*log(energy);
              if (weight < 0) weight = 0;

              sxi  += crystalPos.x()*weight;
              syi  += crystalPos.y()*weight;
              swi  += weight;
        }


        if (swi > 1e-3) return CLHEP::Hep3Vector(sxi/swi,syi/swi,0);
        return CLHEP::Hep3Vector(0,0,0);
    }
    
 
 







}





using mu2e::TrackCaloMatchingExtend;
DEFINE_ART_MODULE(TrackCaloMatchingExtend);
