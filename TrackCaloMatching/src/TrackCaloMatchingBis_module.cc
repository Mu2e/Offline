//
// Original author B. Echenard
//
// A generic track calo matching algorithm. 
//
// The matching is performed in the calorimeter front face section frame, hence only x and y coordinate matter (z is always zero in that frame)
//
// As of today, use ad-hoc errors, will work on including the event-by-event error later
// As of today, use x-y coordinates, might want to shift to u-v


// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Handle.h"

// From the art tool-chain
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "KalmanTests/inc/KalRepPtrCollection.hh"
#include "KalmanTests/inc/TrkFitDirection.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "RecoDataProducts/inc/TrkCaloIntersectCollection.hh"
#include "RecoDataProducts/inc/TrkCaloMatchCollection.hh"

#include "BTrk/TrkBase/HelixParams.hh"
#include "BTrk/TrkBase/TrkRep.hh"
#include "BTrk/TrkBase/HelixTraj.hh"


//CLHEP includes
#include "CLHEP/Vector/TwoVector.h"
#include "BTrk/BbrGeom/HepPoint.h"

// Other includes.
#include "cetlib/exception.h"

#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <functional>

#include "TH1D.h"



namespace mu2e {

  class TrackCaloMatchingBis : public art::EDProducer {

      public:

	  explicit TrackCaloMatchingBis(fhicl::ParameterSet const& pset):
	    _caloClusterModuleLabel(pset.get<std::string>("caloClusterModuleLabel")),
	    _trkIntersectModuleLabel(pset.get<std::string>("trkIntersectModuleLabel")),
	    _trkFitterModuleLabel(pset.get<std::string>("fitterModuleLabel")),
	    _tpart((TrkParticle::type)(pset.get<int>("fitparticle"))),
	    _fdir((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection"))),
	    _diagLevel             (pset.get<int>   ("diagLevel",0)),
	    _minClusterEnergy      (pset.get<double>("minClusterEnergy")),  
	    _maxDeltaT             (pset.get<double>("maxDeltaT")),  
	    _meanInteractionLength (pset.get<double>("meanInteractionLength")),
	    _dtOffset              (pset.get<double>("dtOffset")),  
	    _sigmaT                (pset.get<double>("sigmaT")),  
	    _sigmaX                (pset.get<double>("sigmaX")),  
	    _sigmaY                (pset.get<double>("sigmaY")),  
	    _saveBestOnly          (pset.get<bool>("saveBestOnly")),  
	    _hXY(0),_hT(0)
	  {
              _trkfitInstanceName = _fdir.name() + _tpart.name();
              produces<TrkCaloMatchCollection>();
	  }

	  virtual ~TrackCaloMatchingBis() {}

	  void beginJob();
	  void endJob  () {}
	  void produce (art::Event& e);



      private:

	  void doMatching(TrkCaloMatchCollection& trackClusterMatch, KalRepPtrCollection const& trksPtrColl, 
                          art::Handle<CaloClusterCollection> caloClustersHandle, art::Handle<TrkCaloIntersectCollection>  trjIntersectHandle);

          void matchBest(TrkCaloMatchCollection& trackClusterMatch, art::Handle<CaloClusterCollection> caloClustersHandle, 
                         art::Ptr<TrkCaloIntersect> const& extrapolPtr, CLHEP::Hep3Vector const& posTrkMatch, double trkTime);

          void matchAll(TrkCaloMatchCollection& trackClusterMatch,  art::Handle<CaloClusterCollection> caloClustersHandle,  
                        art::Ptr<TrkCaloIntersect> const& extrapolPtr, CLHEP::Hep3Vector const& posTrkMatch, double trkTime);

	  std::string      _caloClusterModuleLabel;
	  std::string      _trkIntersectModuleLabel;
	  std::string      _trkFitterModuleLabel;
	  std::string      _trkfitInstanceName;
	  TrkParticle      _tpart;
	  TrkFitDirection  _fdir;

	  int              _diagLevel;
	  double           _minClusterEnergy;  
	  double           _maxDeltaT;		
	  double           _meanInteractionLength;
	  double           _dtOffset;
	  double           _sigmaT;
	  double           _sigmaX;
	  double           _sigmaY;
	  bool             _saveBestOnly;

          
	  TH1F *_hXY,*_hT,*_hChi2,*_hChi2Time,*_hChi2Pos;


  };


  void TrackCaloMatchingBis::beginJob()
  {

      art::ServiceHandle<art::TFileService> tfs;
      _hXY = tfs->make<TH1F>("XY_diff","Clu - TRK X-Y dist",100,-400,400.);
      _hT  = tfs->make<TH1F>("T_diff","Clu - TRK Time dist",100,-10,10.);
      _hChi2 = tfs->make<TH1F>("Chi2","Chi2",100,0,1000.);
      _hChi2Time = tfs->make<TH1F>("Chi2Pos","Chi2",100,0,1000.);
      _hChi2Pos = tfs->make<TH1F>("Chi2Time","Chi2",100,0,1000.);

  
  }
  
  //-----------------------------------------------------------------------------
  void TrackCaloMatchingBis::produce(art::Event& event)
  {
	
	art::Handle<KalRepPtrCollection> trksHandle;
	event.getByLabel(_trkFitterModuleLabel,_trkfitInstanceName,trksHandle);
	KalRepPtrCollection const& trksPtrColl(*trksHandle);

	art::Handle<CaloClusterCollection> caloClustersHandle;
	event.getByLabel(_caloClusterModuleLabel, caloClustersHandle);
	CaloClusterCollection const& caloClusters(*caloClustersHandle);

        art::Handle<TrkCaloIntersectCollection>  trjIntersectHandle;
        event.getByLabel(_trkIntersectModuleLabel, trjIntersectHandle);
        TrkCaloIntersectCollection const& trkIntersects(*trjIntersectHandle);
		

	//output of TrackCaloMatchingBis
        std::unique_ptr<TrkCaloMatchCollection> trackClusterMatch(new TrkCaloMatchCollection);


	if (_diagLevel) std::cout<<"Event Number: "<< event.event()<<"\nStart TrkCaloMatching with nTrk="<<trksPtrColl.size()
	                         <<", nExtrapol = "<<trkIntersects.size()<<", nCluster="<<caloClusters.size()<<std::endl;
	
	if (trksPtrColl.size() && caloClusters.size()) doMatching(*trackClusterMatch, trksPtrColl, caloClustersHandle, trjIntersectHandle);

	event.put(std::move(trackClusterMatch));
  }



  //-----------------------------------------------------------------------------
  void TrackCaloMatchingBis::doMatching(TrkCaloMatchCollection& trackClusterMatch, KalRepPtrCollection const& trksPtrColl, 
                                     art::Handle<CaloClusterCollection> caloClustersHandle, art::Handle<TrkCaloIntersectCollection>  trjIntersectHandle)
  {
  
	 Calorimeter const& cal = *(GeomHandle<Calorimeter>());

 	 TrkCaloIntersectCollection const& trkIntersects(*trjIntersectHandle);
	 CaloClusterCollection const& caloClusters(*caloClustersHandle);
	 
	 if (_diagLevel > 1)
	     for (auto const& cluster: caloClusters) std::cout<<"In TrackClusterMatching, cluster energy="<<cluster.energyDep()<<"  time="<<cluster.time()<<"  cog="<<cluster.cog3Vector()<<std::endl;

	 
	   for (unsigned int itrk=0;itrk<trkIntersects.size();++itrk)
	   {
      
               TrkCaloIntersect const& extrapol = trkIntersects.at(itrk);

	       double pathLength             = extrapol.pathLengthEntrance();	       
	       double trkTime                = extrapol.trk()->arrivalTime(pathLength);
	       CLHEP::Hep3Vector trkMomentum = extrapol.trk()->momentum(pathLength);
	       HepPoint          point       = extrapol.trk()->position(pathLength);

	       CLHEP::Hep3Vector posTrkInTracker(point.x(),point.y(),point.z());	     
	       CLHEP::Hep3Vector posTrkInSectionFF = cal.toSectionFrameFF(extrapol.sectionId(),cal.fromTrackerFrame(posTrkInTracker));


	       HelixTraj trkHel(extrapol.trk()->helix(pathLength).params(),extrapol.trk()->helix(pathLength).covariance());

	       double phi0          = trkHel.phi0();
	       double omega         = trkHel.omega();
	       double radius        = 1.0/omega;
	       double cosDip        = trkHel.cosDip();
	       double sinDip        = sqrt(1-cosDip*cosDip);
	       double centerCircleX = (trkHel.d0() + radius)*sin(phi0);
	       double centerCircleY = (trkHel.d0() + radius)*cos(phi0);	      

	       double length = (posTrkInTracker.z()-trkHel.z0())/sinDip;
                      length += _meanInteractionLength;	      

	       double x =  radius*sin(phi0+omega*cosDip*length) - centerCircleX;
	       double y = -radius*cos(phi0+omega*cosDip*length) + centerCircleY;
               double z = posTrkInTracker.z();
               CLHEP::Hep3Vector posTrkMatch(x,y,z);

 	       art::Ptr<TrkCaloIntersect> extrapolPtr = art::Ptr<TrkCaloIntersect>(trjIntersectHandle,itrk);	       
              
	       if (_saveBestOnly) matchBest(trackClusterMatch, caloClustersHandle, extrapolPtr, posTrkMatch,trkTime);
               else               matchAll( trackClusterMatch, caloClustersHandle, extrapolPtr, posTrkMatch,trkTime);
	   
	   }

  }



  void TrackCaloMatchingBis::matchBest(TrkCaloMatchCollection& trackClusterMatch, art::Handle<CaloClusterCollection> caloClustersHandle, 
                                    art::Ptr<TrkCaloIntersect> const& extrapolPtr, CLHEP::Hep3Vector const& posTrkMatch, double trkTime)
  {
        
	CaloClusterCollection const& caloClusters(*caloClustersHandle);
	
	unsigned int iBest(0);double chi2Best(1e6);
	for (unsigned int i=0;i<caloClusters.size();++i)
	{
	    CaloCluster cluster = caloClusters.at(i);

	    double cluX = cluster.cog3Vector().x();
	    double cluY = cluster.cog3Vector().y();
	    double cluT = cluster.time();
	    //double cluE = cluster.energyDep();

	    double chi2 =  pow((cluT-trkTime-_dtOffset)/_sigmaT,2);
		   //chi2 += pow((cluE-trkMomentum.mag())/(0.1*cluE),2);
		   chi2 += pow((cluX-posTrkMatch.x())/_sigmaX,2);
		   chi2 += pow((cluY-posTrkMatch.y())/_sigmaY,2);

	    if (chi2 < chi2Best) {chi2Best=chi2; iBest=i;} 
	 }


	 //noe recalculate the values for the best chi2
	 CaloCluster const* cluster = &(caloClusters.at(iBest));

	 double chi2Time =  pow((cluster->time()-trkTime-_dtOffset)/_sigmaT,2);
	 double chi2Pos  =  pow((cluster->cog3Vector().x()-posTrkMatch.x())/_sigmaX,2);
		chi2Pos  += pow((cluster->cog3Vector().y()-posTrkMatch.y())/_sigmaY,2);


	 CaloCluster const* caloClusterBase = &caloClusters.front();
	 size_t idx = (cluster - caloClusterBase);
         art::Ptr<CaloCluster> clusterPtr = art::Ptr<CaloCluster>(caloClustersHandle,idx);	       

	 TrkCaloMatch match(extrapolPtr,clusterPtr,iBest, chi2Best, chi2Time, chi2Pos);
         trackClusterMatch.push_back(match);     

	 if (_diagLevel > 0) std::cout<<"Best chi2 for track is "<<chi2Best<<" "<<iBest
		                      <<" with Chi2Time= "<<chi2Time<<" and Chi2Pos = "<<chi2Pos<<std::endl;

	 	       
   }
 
 
   void TrackCaloMatchingBis::matchAll(TrkCaloMatchCollection& trackClusterMatch,  art::Handle<CaloClusterCollection> caloClustersHandle,  
                                    art::Ptr<TrkCaloIntersect> const& extrapolPtr, CLHEP::Hep3Vector const& posTrkMatch, double trkTime)
   {

         CaloClusterCollection const& caloClusters(*caloClustersHandle);
	 CaloCluster const* caloClusterBase = &caloClusters.front();

	 for (unsigned int i=0;i<caloClusters.size();++i)
	 {
	     CaloCluster const* cluster = &(caloClusters.at(i));

	     double chi2Time =  pow((cluster->time()-trkTime-_dtOffset)/_sigmaT,2);
	     double chi2Pos  =  pow((cluster->cog3Vector().x()-posTrkMatch.x())/_sigmaX,2);
		    chi2Pos  += pow((cluster->cog3Vector().y()-posTrkMatch.y())/_sigmaY,2);
	     double chi2 = chi2Time + chi2Pos;

	     size_t idx = (cluster - caloClusterBase);
             art::Ptr<CaloCluster> clusterPtr = art::Ptr<CaloCluster>(caloClustersHandle,idx);	       

	     TrkCaloMatch match(extrapolPtr,clusterPtr, i, chi2, chi2Time, chi2Pos);
             trackClusterMatch.push_back(match);     
             
	     if (_diagLevel > 0) std::cout <<"TrackCaloMatchingBis inserted match for icl="<<i<<" trkId="<<extrapolPtr->trkId()
	                                   <<"  cluster pos= "<<cluster->cog3Vector()<<"   trk pos="<<posTrkMatch
  	                                   <<"  with chi2="<<chi2<<" with Chi2Time= "<<chi2Time<<" and Chi2Pos = "<<chi2Pos<<std::endl;
	     
	     _hXY->Fill(cluster->cog3Vector().x()-posTrkMatch.x());
	     _hXY->Fill(cluster->cog3Vector().y()-posTrkMatch.y());
	     _hT->Fill(cluster->time()-trkTime-_dtOffset);
	     _hChi2->Fill(chi2);
	     _hChi2Time->Fill(chi2Time);
	     _hChi2Pos->Fill(chi2Pos);
	  
	  }



    }
   

}

using mu2e::TrackCaloMatchingBis;
DEFINE_ART_MODULE(TrackCaloMatchingBis);
