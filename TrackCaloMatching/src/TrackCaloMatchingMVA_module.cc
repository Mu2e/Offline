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
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"


#include "BTrk/BbrGeom/HepPoint.h"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"

#include "Offline/RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/TrkCaloIntersect.hh"
#include "Offline/RecoDataProducts/inc/TrkCaloMatch.hh"


#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <functional>
#include <numeric>

#include "TH1D.h"
#include "TMVA/Reader.h"




namespace {

   struct BDTinput
   {
      float e0,e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20;
      float t0,t1,r0,z0;
   };
}



namespace mu2e {

   class TrackCaloMatchingMVA : public art::EDProducer {

       public:

           explicit TrackCaloMatchingMVA(fhicl::ParameterSet const& pset):
             art::EDProducer{pset},
             caloClusterModuleLabel_    (pset.get<std::string>("caloClusterModuleLabel")),
             trkIntersectModuleLabel_   (pset.get<std::string>("trkIntersectModuleLabel")),
             BDTWeightsX_               (pset.get<std::string>("BDTWeightsX")),
             BDTWeightsY_               (pset.get<std::string>("BDTWeightsY")),
             dtOffset_                  (pset.get<double>("dtOffset")),
             sigmaT_                    (pset.get<double>("sigmaT")),
             sigmaXY_                   (pset.get<double>("sigmaXY")),
             chi2Cut_                   (pset.get<double>("chi2Cut")),
             diagLevel_                 (pset.get<int>   ("diagLevel",0)),
             readerX_("!Color:!Silent"),
             readerY_("!Color:!Silent"),
             input_()
           {
               produces<TrkCaloMatchCollection>();
           }

           virtual ~TrackCaloMatchingMVA() {}

           void beginJob();
           void endJob  () {}
           void produce (art::Event& e);



       private:

           void setReader(TMVA::Reader& reader, std::string filename);

           void matchMe(TrkCaloMatchCollection& trackClusterMatch, art::Handle<TrkCaloIntersectCollection>  trjIntersectHandle,
                        art::Handle<CaloClusterCollection> caloClustersHandle);

           double calcChi2Pos(const Calorimeter& cal, const CaloCluster& cluster, const CLHEP::Hep3Vector& trkMomentum,
                              const CLHEP::Hep3Vector& posTrkMatch, double trkTime, double cellsize);


           std::string      caloClusterModuleLabel_;
           std::string      trkIntersectModuleLabel_;
           std::string      BDTWeightsX_;
           std::string      BDTWeightsY_;
           double           dtOffset_;
             double           sigmaT_;
           double           sigmaXY_;
           double           chi2Cut_;

           int              diagLevel_;

           TMVA::Reader     readerX_;
           TMVA::Reader     readerY_;
           BDTinput         input_;

           TH1F *_hXY,*_hT,*_hChi2,*_hChi2Time,*_hChi2Pos;


   };


   //------------------------------------
   void TrackCaloMatchingMVA::beginJob()
   {

       art::ServiceHandle<art::TFileService> tfs;
       _hXY       = tfs->make<TH1F>("XY_diff","Clu - TRK X-Y dist",100,-100,100.);
       _hT        = tfs->make<TH1F>("T_diff","Clu - TRK Time dist",100,-10,10.);
       _hChi2     = tfs->make<TH1F>("Chi2",    "Chi2",200,0,400.);
       _hChi2Time = tfs->make<TH1F>("Chi2Pos", "Chi2",200,0,400.);
       _hChi2Pos  = tfs->make<TH1F>("Chi2Time","Chi2",200,0,400.);

       setReader(readerX_,BDTWeightsX_);
       setReader(readerY_,BDTWeightsY_);
   }

   //-----------------------------------------------------------------------------
   void TrackCaloMatchingMVA::produce(art::Event& event)
   {

         art::Handle<CaloClusterCollection> caloClustersHandle;
         event.getByLabel(caloClusterModuleLabel_, caloClustersHandle);
         const CaloClusterCollection& caloClusters(*caloClustersHandle);

         art::Handle<TrkCaloIntersectCollection>  trjIntersectHandle;
         event.getByLabel(trkIntersectModuleLabel_, trjIntersectHandle);
         const TrkCaloIntersectCollection& trkIntersects(*trjIntersectHandle);

         //output of TrackCaloMatchingMVA
         std::unique_ptr<TrkCaloMatchCollection> trackClusterMatch(new TrkCaloMatchCollection);


         if (diagLevel_) std::cout<<"Event Number: "<< event.event()<<"\nStart TrkCaloMatching with nExtrapol = "
                                  <<trkIntersects.size()<<", nCluster="<<caloClusters.size()<<std::endl;

         if (trkIntersects.size() && caloClusters.size()) matchMe(*trackClusterMatch, trjIntersectHandle, caloClustersHandle);

         event.put(std::move(trackClusterMatch));
   }



   //-----------------------------------------------------------------------------
   void TrackCaloMatchingMVA::matchMe(TrkCaloMatchCollection& trackClusterMatch, art::Handle<TrkCaloIntersectCollection>  trjIntersectHandle,
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

            for (const auto& cluster : caloClusters)
           {
               if (trkIntersect.diskId() != cluster.diskID()) continue;

               CLHEP::Hep3Vector diff = cluster.cog3Vector()-posTrkInSectionFF;
               double deltaTime       = std::abs(cluster.time()-trkTime-dtOffset_);

               if (deltaTime > 50) continue;
               //if (sqrt(diff.x()*diff.x()+diff.y()*diff.y()) > 200) continue;

               double chi2Time = pow(deltaTime/sigmaT_,2);
               double chi2Pos  = calcChi2Pos(cal, cluster, trkMomentum, posTrkInSectionFF, trkTime, cellsize);
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

               if (diagLevel_ > 0) std::cout <<"TrackCaloMatchingMVA inserted match for icl="<<iclu<<" trkMAtch="<<itrk
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



   //-------------------------------------------------------------------------------------------
   double TrackCaloMatchingMVA::calcChi2Pos(const Calorimeter& cal, const CaloCluster& cluster,
                                            const CLHEP::Hep3Vector& trkMomentum,
                                            const CLHEP::Hep3Vector& posTrkMatch,
                                            double trkTime, double cellsize)
   {

        const auto& hit0 = cluster.caloHitsPtrVector().at(0);
        CLHEP::Hep3Vector center = cal.geomUtil().mu2eToDiskFF(cal.crystal(hit0->crystalID()).diskID(), cal.crystal(hit0->crystalID()).position());

        std::vector<int> neighbors1 = cal.neighbors(hit0->crystalID());
        std::vector<int> neighbors2 = cal.nextNeighbors(hit0->crystalID());

        double eCells(0);
        std::vector<double> evec(19,0.0);
        evec[0]=hit0->energyDep();

        for (size_t i=0;i<neighbors1.size();++i)
        {
            double eCell(0);
            for (const auto& hit : cluster.caloHitsPtrVector() ){
                if (hit->crystalID() == neighbors1[i]) {eCell=hit->energyDep(); break;}
            }
            evec[i+1] = eCell;
            eCells += eCell;
        }

        for (size_t i=0;i<neighbors2.size();++i)
        {
            double eCell(0);
            for (const auto& hit : cluster.caloHitsPtrVector() ){
                if (hit->crystalID() == neighbors1[i]) {eCell=hit->energyDep(); break;}
            }
            evec[i+7] = eCell;
            eCells += eCell;
        }


        input_.e0 = evec[0]/100;
        input_.e1 = evec[1]/50;
        input_.e2 = evec[2]/50;
        input_.e3 = evec[3]/50;
        input_.e4 = evec[4]/50;
        input_.e5 = evec[5]/50;
        input_.e6 = evec[6]/50;
        input_.e7 = evec[7];
        input_.e8 = evec[8];
        input_.e9 = evec[9];
        input_.e10 = evec[10];
        input_.e11 = evec[11];
        input_.e12 = evec[12];
        input_.e13 = evec[13];
        input_.e14 = evec[14];
        input_.e15 = evec[15];
        input_.e16 = evec[16];
        input_.e17 = evec[17];
        input_.e18 = evec[18];
        input_.e19 = (cluster.energyDep()-eCells)/100;
        input_.e20 = 0;

        input_.t0  = trkMomentum.phi();
        input_.t1  = trkMomentum.theta();
        input_.r0  = sqrt(center.x()*center.x()+center.y()*center.y())/1000;
        input_.z0  = posTrkMatch.z()/200;


        double val1 = (readerX_.EvaluateRegression("BDTG method"))[0];
        double val2 = (readerY_.EvaluateRegression("BDTG method"))[0];
        double valX = cellsize*val1+center.x();
        double valY = cellsize*val2+center.y();
        double chi2 = ( (valX-posTrkMatch.x())*(valX-posTrkMatch.x()) +
                        (valY-posTrkMatch.y())*(valY-posTrkMatch.y()) )/sigmaXY_/sigmaXY_;

        if (diagLevel_ >3)
        {
            std::cout<<"[TrackCaloMatchingMVA::calcChi2Pos] center"<<center<<std::endl;
            std::cout<<"[TrackCaloMatchingMVA::calcChi2Pos]"<<input_.e2<<" "<<input_.e3<<" "<<input_.e4<<" "<<input_.e5
                     <<input_.e6<<" "<<input_.e7<<" "<<input_.e8<<" "<<input_.e9<<" "<<input_.e10<<" "<<input_.e11<<" "<<input_.e12
                     <<" "<<input_.e13<<" "<<input_.e14<<" "<<input_.e15<<" "<<input_.e16<<" "<<input_.e17<<" "<<input_.e18
                     <<" "<<input_.e19<<"    "
                     <<input_.t1<<" "<<input_.t1<<" "<<input_.r0<<std::endl;
            std::cout<<"[TrackCaloMatchingMVA::calcChi2Pos]"<<val1<<" "<<val2<<std::endl;
        }

        if (diagLevel_ > 2) _hXY->Fill(valX-posTrkMatch.x());
        if (diagLevel_ > 2) _hXY->Fill(valY-posTrkMatch.y());

        return chi2;

    }



    //---------------------------------------------------------------------------
    void TrackCaloMatchingMVA::setReader(TMVA::Reader& reader, std::string filename)
    {
        reader.AddVariable(  "e0", &input_.e0);
        reader.AddVariable(  "e1", &input_.e1);
        reader.AddVariable(  "e2", &input_.e2);
        reader.AddVariable(  "e3", &input_.e3);
        reader.AddVariable(  "e4", &input_.e4);
        reader.AddVariable(  "e5", &input_.e5);
        reader.AddVariable(  "e6", &input_.e6);
        reader.AddVariable(  "e7", &input_.e7);
        reader.AddVariable(  "e8", &input_.e8);
        reader.AddVariable(  "e9", &input_.e9);
        reader.AddVariable(  "e10", &input_.e10);
        reader.AddVariable(  "e11", &input_.e11);
        reader.AddVariable(  "e12", &input_.e12);
        reader.AddVariable(  "e13", &input_.e13);
        reader.AddVariable(  "e14", &input_.e14);
        reader.AddVariable(  "e15", &input_.e15);
        reader.AddVariable(  "e16", &input_.e16);
        reader.AddVariable(  "e17", &input_.e17);
        reader.AddVariable(  "e18", &input_.e18);
        reader.AddVariable(  "e19", &input_.e19);
        reader.AddVariable(  "t0", &input_.t0);
        reader.AddVariable(  "t1", &input_.t1);
        reader.AddVariable(  "r0", &input_.r0);
        reader.AddVariable(  "z0", &input_.z0);
        reader.BookMVA("BDTG method", filename );
    }



   /*

        std::vector<double> matrixX,matrixY,matrixE;
        for (const auto& hit : cluster.caloHitsPtrVector() )
        {
           CLHEP::Hep3Vector crystalPos = cal.geomUtil().mu2eToDiskFF(cal.crystal(hit->crystalID()).diskID(), cal.crystal(hit->crystalID()).position());

           matrixX.push_back( (crystalPos.x()-center.x())/cellsize );
           matrixY.push_back( (crystalPos.y()-center.y())/cellsize );
           matrixE.push_back( hit->energyDep()/105.0 );
        }

        if (matrixX.size()<9) for (unsigned int i=0;i<9-matrixX.size();++i){matrixX.push_back(0);matrixY.push_back(0);matrixE.push_back(0);}

        input_.e2  = matrixE[0];input_.e5  = matrixE[1];input_.e8  = matrixE[2];input_.e11 = matrixE[3];
        input_.e14 = matrixE[4];input_.e17 = matrixE[5];input_.e20 = matrixE[6];input_.e23 = matrixE[7];
        input_.e3  = matrixX[1];input_.e6  = matrixX[2];input_.e9  = matrixX[3];input_.e12 = matrixX[4];
        input_.e15 = matrixX[5];input_.e18 = matrixX[6];input_.e21 = matrixX[7];
        input_.e4  = matrixY[1];input_.e7  = matrixY[2];input_.e10 = matrixY[3];input_.e13 = matrixY[4];
        input_.e16 = matrixY[5];input_.e19 = matrixY[6];input_.e2 = matrixY[7];
        input_.t0  = trkMomentum.phi();
        input_.t1  = trkMomentum.theta();
        input_.r0  = sqrt(center.x()*center.x()+center.y()*center.y())/1000;

    */


}





using mu2e::TrackCaloMatchingMVA;
DEFINE_ART_MODULE(TrackCaloMatchingMVA)
