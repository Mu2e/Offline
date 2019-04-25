//Note : could take advantage of sorting straw hits by time in the mergeCluster module, 
//                                     
// Object to perform helix fit to straw hits
//
#include "TrkReco/inc/TNTClusterer.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "BFieldGeom/inc/BFieldConfig.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include "GeneralUtilities/inc/Angles.hh"

#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "cetlib_except/exception.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <functional>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/weighted_median.hpp>

#include "TTree.h"
#include "Math/VectorUtil.h"

using namespace boost::accumulators;
using namespace ROOT::Math::VectorUtil;

// Note 1: The list allows you to take pointer of its elements, a vector will reallocate and invalidate the pointers
//         The code with the vector version is much more complicated and has a tiny performwnce improvement.
//         The initDB also makes use of the list specificities to delete empty clusters


namespace mu2e
{
// std::sort sorts in ascending order, so the following will put the best hits first
  struct BkgComp : public std::binary_function<BkgClusterHit,BkgClusterHit, bool> {
    BkgComp(ComboHitCollection const& chcol) : _chcol(chcol){}
    bool operator()(BkgClusterHit const& x, BkgClusterHit const& y) {
      return _chcol[x.index()].wireRes() < _chcol[y.index()].wireRes();
//      return _chcol[x.index()].nStrawHits() > _chcol[y.index()].nStrawHits();
    }
    ComboHitCollection const& _chcol;
  };

  int ClusterStraw::idCounter_=0;

   ClusterStraw::ClusterStraw(BkgClusterHit& hit, const ComboHit& chit) : 
      _id(idCounter_++),_pos(chit.pos()),_time(chit.time()),_hasChanged(true),_hitsPtr(std::vector<BkgClusterHit*>(1,&hit))
   {}


   void ClusterStraw::updateCache(const ComboHitCollection& chcol, float maxwt)
   {             
       if (_hitsPtr.size()==1)
       {
          int idx = _hitsPtr.at(0)->index();          
          _pos = chcol[idx].pos();
          _pos.SetZ(0);
          _time = chcol[idx].time();
          return;
       } 
              
       accumulator_set<float, stats<tag::weighted_median>, unsigned> racc, pacc, tacc;      
             
       float crho  = sqrtf(_pos.perp2());
       float cphi  = _pos.phi();
//       XYZVec rdir = PerpVector(_pos,Geom::ZDir()).unit();
//       XYZVec pdir(-rdir.y(),rdir.x(),0.0);

       for (auto& hitPtr : _hitsPtr)
       {
          int idx         = hitPtr->index();
          float hphi      = chcol[idx].phi();
          
          float dt  = chcol[idx].time() -_time;
	  float dr  = sqrtf(chcol[idx].pos().perp2()) - crho;
	  float dp = Angles::deltaPhi(hphi,cphi);

          // weight according to the # of hits

	  tacc(dt,weight=chcol[idx].nStrawHits());
	  racc(dr,weight=chcol[idx].nStrawHits());
	  pacc(dp,weight=chcol[idx].nStrawHits());
       }

       crho  += extract_result<tag::weighted_median>(racc);
       cphi  += extract_result<tag::weighted_median>(pacc);
       _time += extract_result<tag::weighted_median>(tacc);       
       _pos   = XYZVec(crho*cos(cphi),crho*sin(cphi),0.0);
   }


   //---------------------------------------------------------------------------------------
   TNTClusterer::TNTClusterer(const fhicl::ParameterSet& pset) :
     _diag(pset.get<int>(                         "diagLevel",0)),
     _testflag(pset.get<bool>(                    "TestFlag")),
     _bkgmask(pset.get<std::vector<std::string> >("BackgroundMask",std::vector<std::string>())),
     _sigmask(pset.get<std::vector<std::string> >("SignalMask",std::vector<std::string>())),
     _comboInit(pset.get<bool>(                   "ComboInit",true)),         
     _mergeInit(pset.get<bool>(                   "MergeInit",false)),        
     _dseed(pset.get<float>(                      "SeedDistance")),           
     _dhit(pset.get<unsigned>(                    "HitDistance")),            
     _dd(pset.get<float>(                         "ClusterDiameter",10.0)),   
     _dt(pset.get<float>(                         "TimeDifference",30.0)),    
     _maxdt(pset.get<float>(                      "MaxTimeDifference")),      
     _maxdsum(pset.get<float>(                    "MaxDistanceSum",100.0)),   
     _maxNiter(pset.get<unsigned>(                "MaxNIterations")),
     _maxNchanged(pset.get<unsigned>(             "MaxNChanged",2)),
     _hitIndex(200,std::vector<ClusterStraw*>()),
     _cptrs()
   {
       // cache some values
       float minerr(pset.get<float>( "MinHitError",5.0));
       float maxdist(pset.get<float>("MaxDistance",50.0));
       float trms(pset.get<float>(   "TimeRMS",2.0));            

       _trms2inv  = 1.0/trms/trms;     
       _ditime = int(_maxdt/10.0)+1;
       _dd2 = _dd*_dd;
       _maxwt = 1.0/minerr;
       _md2 = maxdist*maxdist;
   }


   //---------------------------------------------------------------------------------------
   void TNTClusterer::init()
   {
      if (_diag > 0)
      {
	art::ServiceHandle<art::TFileService> tfs;
	_idiag = tfs->make<TTree>("idiag","iteration diagnostics");
	_idiag->Branch("nhits",      &nhits_,        "nhits/I");
	_idiag->Branch("hitRad",     &hitRad_,       "hitRad[nhits]/F");
	_idiag->Branch("hitPhi",     &hitPhi_,       "hitPhi[nhits]/F");
	_idiag->Branch("hitTime",    &hitTime_,      "hitTime[nhits]/F");      
	_idiag->Branch("hitNcombo",  &hitNcombo_,    "hitNcombo[nhits]/I");      

	_idiag->Branch("ncluIter",   &ncluIter_,     "ncluIter/I");
	_idiag->Branch("cluId",      &cluId_,        "cluId[ncluIter]/I");
	_idiag->Branch("cluNpass",   &cluNpass_,     "cluNpass[ncluIter]/I");      
	_idiag->Branch("cnhi",       &cluNhit_,      "cnhi[ncluIter]/I");      
	_idiag->Branch("cluRad",     &cluRad_,       "cluRad[ncluIter]/F");
	_idiag->Branch("cluPhi",     &cluPhi_,       "cluPhi[ncluIter]/F");
	_idiag->Branch("cluTime",    &cluTime_,      "cluTime[ncluIter]/F");	
	_idiag->Branch("cluR2diff",  &cluR2diff_,    "cluR2diff[ncluIter]/F");	
	_idiag->Branch("cluRdiff",   &cluRdiff_,     "cluRdiff[ncluIter]/F");	
	_idiag->Branch("cluPdiff",   &cluPdiff_,     "cluPdiff[ncluIter]/F");	
	_idiag->Branch("cluTdiff",   &cluTdiff_,     "cluTdiff[ncluIter]/F");	
	_idiag->Branch("cluTRdiff",  &cluTRdiff_,    "cluTRdiff[ncluIter]/F");	

	_idiag->Branch("nhitClu",    &nhitClu_,      "nhitClu/I");
	_idiag->Branch("hcIdxClu",   &hcIdxClu_,     "hcIdxClu[nhitClu]/I");
	_idiag->Branch("hcIdxHit",   &hcIdxHit_,     "hcIdxHit[nhitClu]/I");
	_idiag->Branch("hcNpass",    &hcNpass_,      "hcNpass[nhitClu]/I");      

	_idiag->Branch("niter",      &niter_,        "niter/I");
        _idiag->Branch("nclu",       &nclu_,         "nclu[niter]/I");
	_idiag->Branch("nChanged",   &nChanged_,     "nChanged[niter]/I");
	_idiag->Branch("odist",      &odist_,        "odist[niter]/F");
	_idiag->Branch("tdist",      &tdist_,        "tdist[niter]/F");

      }
   }









   //-------------------------------------------------------------------------------------------------------------------
   void TNTClusterer::findClusters(BkgClusterCollection& clusterColl, const ComboHitCollection& chcol)
   {
 
      ClusterStraw::resetCounter();      
      if (_diag) {ncluIter_ = nhitClu_ = 0;fillHitTree(chcol);}

      
      //reset stuff
      _cptrs = std::vector<ClusterStraw*>(chcol.size(),nullptr);
      for (auto& vec: _hitIndex) vec.clear();
      std::list<ClusterStraw> clusters; //see Note 1


      // loop over the straw hits and create ClusterHits          
      std::vector<BkgClusterHit> chits;
      chits.reserve(chcol.size());

      
      algo1(chcol, clusters, chits);
      //algo2(chcol, clusters, chits);
      
      if (_diag) _idiag->Fill();

      //form and insert final BkgClusters here
      for (auto& cluster: clusters)
      {            
         if (cluster.hits().size()==0) continue;
         clusterColl.emplace_back(BkgCluster(cluster.pos(), cluster.time()));

         BkgCluster& sclust = clusterColl.back();
         sclust._hits.reserve(cluster.hits().size());
         for (const auto& hitPtr : cluster.hits()) sclust._hits.emplace_back(std::move(*hitPtr));
      }
   }
   
  
   //----------------------------------------------------------------------------------------------
   void TNTClusterer::algo1(const ComboHitCollection& chcol, std::list<ClusterStraw>& clusters, std::vector<BkgClusterHit>& chits)
   {
                            
       if (_mergeInit) initCluMerge(chcol,chits, clusters);
       else            initClu(chcol,chits);
                            
       float odist(2.0*_maxdsum); 
       float tdist(0.0);

       unsigned  niter(0);
       while (std::abs(odist - tdist) > _maxdsum  && niter < _maxNiter)
       {        
	   unsigned nChanged = formClusters(chcol, chits, clusters);

           odist = tdist;      
           tdist = 0.0;
           for (auto& cluster: clusters)
             for ( auto& hitPtr : cluster.hits()) tdist += hitPtr->distance();

	   if (_diag) fillCluTree(chcol,clusters,niter,odist,tdist,nChanged);
	   ++niter;        
       }     
      //if you want to do a final merge only, and comment the merging step in formclusters
      //mergeClusters(clusters, _dt, _dd2);

   }
   
   
   //----------------------------------------------------------------------------------------------
   void TNTClusterer::algo2(const ComboHitCollection& chcol, std::list<ClusterStraw>& clusters, std::vector<BkgClusterHit>& chits)
   {                            
       if (_mergeInit) initCluMerge(chcol,chits, clusters);
       else            initClu(chcol,chits);
                            
       unsigned  niter(0),nChanged(clusters.size());
       while (nChanged > _maxNchanged  && niter < _maxNiter)
       {        
	   nChanged = formClusters(chcol, chits,clusters);
	   if (_diag) fillCluTree(chcol, clusters, niter, 0, 0, nChanged);
	   ++niter;        
       }     
      //if you want to do a final merge only, and comment the merging step in formclusters
      //mergeClusters(clusters, _dt, _dd2);
   }
   


   //-------------------------------------------------------------------------------------------------------------------
   void TNTClusterer::initClu(const ComboHitCollection& chcol, std::vector<BkgClusterHit>& chits) 
   {      
     for (size_t ish=0; ish<chcol.size(); ++ish){         
       if (_testflag && (!chcol[ish].flag().hasAllProperties(_sigmask) || chcol[ish].flag().hasAnyProperty(_bkgmask))) continue;           
       chits.emplace_back(BkgClusterHit(ish,chcol[ish].flag()));         
     }
     if (_comboInit) 
     {
       BkgComp chsort(chcol);
       std::sort(chits.begin(),chits.end(),chsort);
     } 
   }  


   //-------------------------------------------------------------------------------------------------------------------
   void TNTClusterer::initCluMerge(const ComboHitCollection& chcol, std::vector<BkgClusterHit>& chits, std::list<ClusterStraw>& clusters) 
   {
       //merge init only uses combohit with two hits as starting point
       for (size_t ish=0;ish<chcol.size();++ish)
       {
          if (_testflag && (!chcol[ish].flag().hasAllProperties(_sigmask) || chcol[ish].flag().hasAnyProperty(_bkgmask))) continue;           
          chits.emplace_back(BkgClusterHit(ish,chcol[ish].flag()));   
          if (chcol.at(ish).nCombo()>1) clusters.emplace_back(ClusterStraw(chits.back(),chcol[ish]));
       }          
       
       mergeClusters(clusters, chcol, _maxdt, _md2); 
      
       for (auto& cluster : clusters)
       {
           _hitIndex[cluster.itime()].emplace_back(&cluster);

           if (cluster.hits().size()==1) 
             cluster.hits().at(0)->distance(0); 
           else 
             for (auto& hit : cluster.hits()) {hit->distance(distance(cluster,chcol[hit->index()])); _cptrs[hit->index()] = &cluster;}                     
       }
   }        
   


   //-------------------------------------------------------------------------------------------------------------------
   void TNTClusterer::mergeClusters(std::list<ClusterStraw>& clusters, const ComboHitCollection& chcol, float dt, float dd2)
   {
       unsigned niter(0);    
       while (niter < _maxNiter)
       {
          int nchanged(0);
          for (auto it1 = clusters.begin(); it1 != std::prev(clusters.end()); ++it1)
          {
              for (auto it2 = std::next(it1); it2!=clusters.end(); ++it2)
              {
                 //if (niter==0 && it2->time() - it1->time() >dt) break; //for sorted times 
                 if (std::abs(it1->time() - it2->time()) > dt) continue;
		 if ((it1->pos() - it2->pos()).perp2() > dd2) continue;

                 it1->flagChanged(true); 
                 it2->flagChanged(true); 
                 ++nchanged;
                 if (it1->hits().size() >= it2->hits().size()) 
                     {mergeTwoClu(*it1,*it2); it2=clusters.erase(it2); --it2;} //remove it2 from list
                 else 
                     {mergeTwoClu(*it2,*it1); it1=clusters.erase(it1); --it1; break;} //remove it1 from list
	      }             
           }	

           ++niter;
           if (nchanged==0) break;
           
           clusters.remove_if([](auto& cluster){return cluster.hits().empty();});
           for (auto& cluster : clusters )
           {
              if (cluster.hasChanged()) cluster.updateCache(chcol, _maxwt);
              cluster.flagChanged(false); 
           }
           if (_diag>0) std::cout<<"Merge "<<niter<<" "<<nchanged<<std::endl;
       }

       return;    
   }

   void TNTClusterer::mergeTwoClu(ClusterStraw& clu1, ClusterStraw& clu2 )
   {
       if (clu1.hits().empty() || clu2.hits().empty()) return;
       clu1.hits().insert(clu1.hits().end(),clu2.hits().begin(), clu2.hits().end());
       clu2.hits().clear();
       clu1.flagChanged(true);
   }


   //-------------------------------------------------------------------------------------------------------------------
   // loop over hits, re-affect them to their original cluster if they are still within the radius, otherwise look at 
   // candidate clusters to check if they could be added. If not, make a new cluster.
   // to speed up, do not update clusters who haven't changed and keep a list of clusters within a given time window
   unsigned  TNTClusterer::formClusters(const ComboHitCollection& chcol, std::vector<BkgClusterHit>& chits, std::list<ClusterStraw>& clusters)
   {     
       unsigned  nchanged(0);
       for (auto& cluster : clusters) cluster.hits().clear();

       for (auto& chit : chits)
       {                        
          //if (chit._nchanged > 5) continue;

          //check if the current hit is still ok or find the best cluster
          if (_cptrs[chit.index()] != nullptr && chit.distance() < _dhit) 
          { 
              _cptrs[chit.index()]->hits().emplace_back(&chit); 
              continue;
          }
          

          ClusterStraw* minc(nullptr);                  
          float mindist(FLT_MAX);         
          int hitIdx  = chit.index();
          int itime = int(chcol[hitIdx].time()/10.0);

          for (int i=std::max(0,itime-_ditime);i<itime+_ditime;++i)
          {
             for (const auto& ic : _hitIndex[i])
             {                
                 float dist = distance(*ic,chcol[hitIdx]);
                 if (dist < mindist) {mindist = dist; minc = ic;}
                 if (mindist < _dhit) break;               
             }          
             if (mindist < _dhit) break;               
          }

          
          
          if (mindist < _dhit) 
          {
              minc->hits().emplace_back(&chit); //add to best cluster
          } 
          else if (mindist > _dseed)
          {
              clusters.emplace_back(ClusterStraw(chit, chcol[hitIdx])); 
              minc = &clusters.back();              
              _hitIndex[minc->itime()].emplace_back(minc);
          } 
          else 
          {
               mindist =10000;
               minc = nullptr;
          }

          
          
          if (minc != nullptr)
          {             
               //associated to a new cluster, need to flag old/new  cluster accordingly
               if (_cptrs[hitIdx] != minc) 
               { 
                  ++nchanged; 
                  //++chit._nchanged; 
                  if (_cptrs[hitIdx] != nullptr) _cptrs[hitIdx]->flagChanged(true); 
                  minc->flagChanged(true);
               }
               _cptrs[hitIdx] = minc;
           } 
           else 
           {
               //unassociated to previous cluster, need to flag old cluster accordingly
               if (_cptrs[hitIdx] != nullptr) 
               {
                  ++nchanged; 
                  //++chit._nchanged; 
                  _cptrs[hitIdx]->flagChanged(true);
               }
               _cptrs[hitIdx] = nullptr;
           }      
       }

       
       
       //clear cluster timing vector, update clusters and hit distance and refil timing vector
       for (auto& vec: _hitIndex) vec.clear(); 
       
       for(auto& cluster : clusters)
       {
           if (cluster.hasChanged())
           {
              cluster.updateCache(chcol, _maxwt);
              cluster.flagChanged(false); 

              if (cluster.hits().size()==1) 
                cluster.hits().at(0)->distance(0); 
              else
                for (auto& hit : cluster.hits()) hit->distance(distance(cluster,chcol[hit->index()]));                      
           }

           _hitIndex[cluster.itime()].emplace_back(&cluster);
       }
       
       return nchanged;
   }








   //---------------------------------------------------------------------------------------
   // only count differences if they are above the natural hit size (drift time, straw size)      
   float TNTClusterer::distance(const ClusterStraw& cluster, const ComboHit& hit) const 
   {     
       float dt = std::abs(hit.time()-cluster.time());
       if (dt > _maxdt) return _dseed+1.0;           

       XYZVec psep = PerpVector(hit.pos()-cluster.pos(),Geom::ZDir());
       float d2 = psep.mag2();
       if (d2 > _md2) return _dseed+1.0; 

       float retval(0.0);
       if (dt > _dt) {float tdist = dt -_dt;  retval = tdist*tdist*_trms2inv;}      
       if (d2 > _dd2) 
       {	
	   XYZVec that(-hit.wdir().y(),hit.wdir().x(),0.0);
           float dw = std::max(0.0f,hit.wdir().Dot(psep)-_dd)/hit.posRes(ComboHit::wire);
	   float dp = std::max(0.0f,that.Dot(psep)-_dd)*_maxwt;  //maxwt = 1/minerr
	   retval += dw*dw + dp*dp;
       }      
       return retval;
   }





   //-------------------------------------------------------------------------------------------------------------------
   void TNTClusterer::dump(std::list<ClusterStraw> clusters)
   {
       int iclu(0);      
       clusters.sort([](ClusterStraw& a, ClusterStraw& b){return a.hits().at(0)->index() < b.hits().at(0)->index() ;});
       for (auto& cluster: clusters)
       { 
           std::cout<<"Cluster "<<iclu<<" "<<cluster.pos()<<" "<<cluster.time()<<"  "<<cluster.hits().size()<<"  - ";
           for (auto& hit : cluster.hits()) std::cout<<hit->index()<<" ";
           std::cout<<std::endl;
           ++iclu;
       }
   }

   //---------------------------------------------------------------------------------------
   void TNTClusterer::fillHitTree(const ComboHitCollection& chcol)
   {
      nhits_ = 0;
      
      for (size_t ish = 0; ish < chcol.size(); ++ish)
      {	
	 hitRad_[nhits_]    = sqrtf(chcol.at(ish).pos().perp2());
	 hitPhi_[nhits_]    = chcol.at(ish).pos().phi();
	 hitTime_[nhits_]   = chcol.at(ish).time();
	 hitNcombo_[nhits_] = chcol.at(ish).nCombo();	
	 hmap_[ish]         = nhits_;
	 ++nhits_;  
      }
      
   }


   //---------------------------------------------------------------------------------------
   void TNTClusterer::fillCluTree(const ComboHitCollection& chcol, std::list<ClusterStraw>& clusters, 
                                   int npass, float odist, float tdist, int nChanged)
   {

      int nclu(0);
      for (const auto& cluster : clusters)
      {         
	 float rdiff(-1.0),r2diff(-1.0),pdiff(-1.0),tdiff(-1.0),trdiff(-1.0);
	 for (unsigned i=0;i<cluster.hits().size() ;++i)
	 {
	    int ishIdx = cluster.hits().at(i)->index();
	    
            float dr = sqrtf(chcol.at(ishIdx).pos().perp2()) - sqrtf(cluster.pos().perp2());
	    float dp = chcol.at(ishIdx).pos().phi()-cluster.pos().phi();
	    float dt = chcol.at(ishIdx).time()-cluster.time();
	    rdiff    = std::max(rdiff,std::abs(dr)); 
  	    pdiff    = std::max(pdiff,std::abs(dp)); 
  	    tdiff    = std::max(tdiff,std::abs(dt)); 
  	    r2diff   = std::max(r2diff,sqrt((chcol.at(ishIdx).pos()-cluster.pos()).perp2())); 

            XYZVec psep = PerpVector(chcol.at(ishIdx).pos()-cluster.pos(),Geom::ZDir());
            float dw    = std::max(float(0.0),chcol.at(ishIdx).wdir().Dot(psep)-_dd)/chcol.at(ishIdx).posRes(ComboHit::wire);
	    XYZVec that(-chcol.at(ishIdx).wdir().y(),chcol.at(ishIdx).wdir().x(),0.0);
	    float dpe   = std::max(float(0.0),that.Dot(psep)-_dd)*_maxwt;  //maxwt = 1/minerr
	    float dtr   =  dw*dw + dpe*dpe;  	    
            trdiff = std::max(trdiff,dtr); 

	    hcIdxClu_[nhitClu_] = ncluIter_;
	    hcIdxHit_[nhitClu_] = hmap_[ishIdx];
	    hcNpass_[nhitClu_]  = npass;
	    ++nhitClu_;
         }
	 if (cluster.hits().size()>0) ++nclu;
          
         cluId_[ncluIter_]    = cluster.id();
         cluNpass_[ncluIter_] = npass;
	 cluNhit_[ncluIter_]  = cluster.hits().size();
         cluRad_[ncluIter_]   = sqrtf(cluster.pos().perp2());
         cluPhi_[ncluIter_]   = cluster.pos().phi();
         cluTime_[ncluIter_]  = cluster.time();
         cluR2diff_[ncluIter_] = r2diff;
         cluRdiff_[ncluIter_]  = rdiff;
         cluPdiff_[ncluIter_]  = pdiff;
         cluTdiff_[ncluIter_]  = tdiff;
         cluTRdiff_[ncluIter_] = trdiff;

         ++ncluIter_;
      }
      niter_           = npass;
      nclu_[npass]     = nclu;
      odist_[npass]    = odist;
      tdist_[npass]    = tdist;
      nChanged_[npass] = nChanged;      

   }



}


// Note 1: The list allows you to take pointer of its elements, a vector will reallocate and invalidate the pointers
//         The code with the vector version is much more complicated and has a very small performwnce improvement.


