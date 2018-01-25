//                                     
// Object to perform helix fit to straw hits
//
#include "TrkReco/inc/TNTBClusterer.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "BFieldGeom/inc/BFieldConfig.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include "GeneralUtilities/inc/Angles.hh"

#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "cetlib_except/exception.h"

#include <iostream>
#include <vector>
#include <algorithm>
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

   ClusterStraw2::ClusterStraw2(ClusterStrawHit2& hit) : 
      _pos(hit._pos),_time(hit._time),_itime(int(_time/10.0)), _hitsPtr(std::vector<ClusterStrawHit2*>(1,&hit)),
      _hasChanged(true)
   {}

   void ClusterStraw2::updateCache(double _maxwt)
   {      
       if (_hitsPtr.size()==1)
       {
          _pos = _hitsPtr.at(0)->_pos;
          _time = _hitsPtr.at(0)->_time;
          _itime = int(_time/10.0);
          return;
       } 
       
       accumulator_set<double, stats<tag::weighted_median(with_p_square_quantile) >, double > racc, pacc, tacc;      

       double crho = sqrt(_pos.perp2());
       double cphi = _pos.phi();
       XYZVec rdir = PerpVector(_pos,Geom::ZDir()).unit();
       XYZVec pdir(-rdir.y(),rdir.x(),0.0);

       for (auto& hitPtr : _hitsPtr)
       {
          double dt  = hitPtr->_time -_time;
	  double dr  = sqrt(hitPtr->_pos.perp2()) - crho;
	  double dp  = hitPtr->_phi-cphi;
          while(dp > 3.1415926)  dp -= 6.2831853;
          while(dp < -3.1415926) dp += 6.2831853;

          // weight according to the wire direction error, linearly for now
	  double twt = std::min(_maxwt,hitPtr->_posResInv);	
	  double rwt = std::min(_maxwt,std::abs(rdir.Dot(hitPtr->_wdir))*hitPtr->_posResInv);
	  double pwt = std::min(_maxwt,std::abs(pdir.Dot(hitPtr->_wdir))*hitPtr->_posResInv);
	  tacc(dt,weight=twt);
	  racc(dr,weight=rwt);
	  pacc(dp,weight=pwt);
       }

       crho  += extract_result<tag::weighted_median>(racc);
       cphi  += extract_result<tag::weighted_median>(pacc);
       _time += extract_result<tag::weighted_median>(tacc);

       _pos = XYZVec(crho*cos(cphi),crho*sin(cphi),0.0);
       _itime = int(_time/10.0);
   }



   TNTBClusterer::TNTBClusterer(const fhicl::ParameterSet& pset) :
     _diag(pset.get<int>("diagLevel",0)),
     _bkgmask(pset.get<std::vector<std::string> >("BackgroundMask",std::vector<std::string>())),
     _sigmask(pset.get<std::vector<std::string> >("SignalMask",std::vector<std::string>())),
     _stereoInit(pset.get<bool>("StereoInit",true)),        // # of sigma to define a new cluster
     _dseed(pset.get<double>("SeedDistance")),     // # of sigma to define a new cluster
     _dhit(pset.get<unsigned>("HitDistance")),       // # of sigma to add hits
     _dd(pset.get<double>("ClusterDiameter",10.0)),   // mm: the natural cluster size
     _dt(pset.get<double>("TimeDifference",30.0)),    // nsec: the natural time spread
     _maxdt(pset.get<double>("MaxTimeDifference")), // Maximum time difference
     _maxdsum(pset.get<double>("MaxDistanceSum",100.0)),   // iteration convergence
     _maxdist(pset.get<double>("MaxDistance")),       // Maximum transverse distance (mm)    
     _maxniter(pset.get<unsigned>("MaxNIterations")),
     _maxnchanged(pset.get<unsigned>("MaxNChanged",5)),
     _hitIndex(200,std::vector<ClusterStraw2*>()),
     _cptrs()
   {
       // cache a maximum of values
       double minerr(pset.get<double>("MinHitError",5.0));
       double maxdist(pset.get<double>("MaxDistance",50.0));
       double trms(pset.get<double>("TimeRMS",2.0));            

       _trms2inv  = 1.0/trms/trms;     
       _ditime = int(_maxdt/10.0)+1;
       _dd2 = _dd*_dd;
       _maxwt = 1.0/minerr;
       _md2 = maxdist*maxdist;
   }


   void TNTBClusterer::init()
   {
      if (_diag > 0){
         art::ServiceHandle<art::TFileService> tfs;
         _idiag = tfs->make<TTree>("idiag","iteration diagnostics");
         _idiag->Branch("niter", &_niter,"niter/i");
         _idiag->Branch("nclu",  &_nclu, "nclu/i");
         _idiag->Branch("nhits", &_nhits,"nhits/i");
         _idiag->Branch("nhits", &_nchits,"nchits/i");
         _idiag->Branch("odist", &_odist,"odist/D");
         _idiag->Branch("tdist", &_tdist,"tdist/D");
      }
   }



   // distance for pre-clustering
   double TNTBClusterer::distance2(const ClusterStraw2& cluster, ClusterStrawHit2& hit) const 
   {     
       double dt = std::abs(hit._time-cluster.time());
       if (dt > _maxdt) return _dseed+1.0;           

       double d2 = (hit._pos-cluster.pos()).perp2();
       if (d2 > _md2) return _dseed+1.0; 

       return 0;
   }

   // only count differences if they are above the natural hit size (drift time, straw size)      
   double TNTBClusterer::distance(const ClusterStraw2& cluster, ClusterStrawHit2& hit) const 
   {     

       double dt = std::abs(hit._time-cluster.time());
       if (dt > _maxdt) return _dseed+1.0;           

       XYZVec psep = PerpVector(hit._pos-cluster.pos(),Geom::ZDir());
       double d2 = psep.mag2();
       if (d2 > _md2) return _dseed+1.0; 

       double retval(0.0);
       if (dt > _dt) {double tdist = dt -_dt;  retval = tdist*tdist*_trms2inv;}      
       if (d2 > _dd2) 
       {	
           double dw = std::max(0.0,hit._wdir.Dot(psep)-_dd)*hit._posResInv;
	   XYZVec that(-hit._wdir.y(),hit._wdir.x(),0.0);
	   double dp = std::max(0.0,that.Dot(psep)-_dd)*_maxwt;  //maxwt = 1/minerr
	   retval += dw*dw + dp*dp;
       }      
       return retval;
   }



   void TNTBClusterer::findClusters(BkgClusterCollection& clusterColl,
                                   const ComboHitCollection& chcol)
   {
      StrawHitFlag stflag(StrawHitFlag::stereo);    
      //reset stuff
      _cptrs = std::vector<ClusterStraw2*>(chcol.size(),nullptr);
      for (auto& vec: _hitIndex) vec.clear();
      std::list<ClusterStraw2> clusters; //see Note 1


      // loop over the straw hits and create ClusterHits          
      std::vector<ClusterStrawHit2> chits;
      chits.reserve(chcol.size());

      
      for(size_t ish=0;ish<chcol.size();++ish)
         if (chcol[ish].flag().hasAllProperties(_sigmask) && !chcol[ish].flag().hasAnyProperty(_bkgmask) && chcol[ish].flag().hasAllProperties(stflag))
                 chits.emplace_back(ClusterStrawHit2(ish,chcol[ish],_srms2inv));
      
      if (chits.size()==0) return;

      _niter = 0;
      _odist=2.0*_maxdsum; 
      _tdist=0.0;
      //int nChanged(chits.size());
      
      while (std::abs(_odist - _tdist) > _maxdsum  && _niter < _maxniter)
      //while (std::abs(_odist - _tdist) > _maxdsum  && nChanged > 10  && _niter < _maxniter)
      {        
          /*nChanged =*/ formClusters(chits,clusters, true);
          _odist = _tdist;      
          _tdist = 0.0;
          for (auto& cluster: clusters)
            for ( auto& hitPtr : cluster.hits()) _tdist += hitPtr->_dist;
          //std::cout<<"Iter "<<_niter<<" "<<nChanged<<std::endl;
          ++_niter;        
      } 






      for(size_t ish=0;ish<chcol.size();++ish)
         if (chcol[ish].flag().hasAllProperties(_sigmask) && !chcol[ish].flag().hasAnyProperty(_bkgmask) && !chcol[ish].flag().hasAllProperties(stflag))
                 chits.emplace_back(ClusterStrawHit2(ish,chcol[ish],_nsrms2inv));

      _nhits = chits.size();
      _niter = 0;
      _odist=2.0*_maxdsum; 
      _tdist=0.0;
      //int nChanged(chits.size());
      while (std::abs(_odist - _tdist) > _maxdsum  && _niter < _maxniter)
      //while (std::abs(_odist - _tdist) > _maxdsum  && nChanged > 10  && _niter < _maxniter)
      {        
          /*nChanged =*/ formClusters(chits,clusters, false);

          _odist = _tdist;      
          _tdist = 0.0;
          for (auto& cluster: clusters)
            for ( auto& hitPtr : cluster.hits()) _tdist += hitPtr->_dist;

          if (_diag>0)
          {
             _nclu = clusters.size();
             _nchits = 0;
             for (auto& cluster: clusters) _nchits += cluster.hits().size();
             _idiag->Fill();
          }
          //std::cout<<"Iter "<<_niter<<" "<<nChanged<<std::endl;
          ++_niter;        
      } 


      //form and insert final BkgClusters here
      for (auto& cluster: clusters)
      {            
         if (cluster.hits().size()==0) continue;
         clusterColl.emplace_back(BkgCluster(cluster.pos(), cluster.time()));

         BkgCluster& sclust = clusterColl.back();
         for (const auto& hitPtr : cluster.hits())
           sclust._hits.emplace_back(BkgClusterHit(hitPtr->_dist, hitPtr->_index, chcol[hitPtr->_index].flag()));
      }
   }





   //-------------------------------------------------------------------------------------------------------------------
   // loop over hits, re-affect them to their original cluster if they are still within the radius, otherwise look at 
   // candidate clusters to check if they could be added. If not, make a new cluster.
   // to speed up, do not update clusters who haven't changed and keep a list of clusters within a given time window
   unsigned  TNTBClusterer::formClusters(std::vector<ClusterStrawHit2>& chits, std::list<ClusterStraw2>& clusters, bool init=false)
   {     
       unsigned  nchanged(0);
       for (auto& cluster : clusters) cluster.hits().clear();


       for (auto& chit : chits)
       {                        
          if (chit._nchanged > 5) continue;

          //check if the current hit is still ok or find the best cluster
          if (_cptrs[chit._index] != nullptr && chit._dist < _dhit) 
          { 
              _cptrs[chit._index]->hits().emplace_back(&chit); 
              continue;
          }

          double mindist(FLT_MAX);         
          ClusterStraw2* minc(nullptr);                  
          for (int i=std::max(0,chit._itime-_ditime);i<chit._itime+_ditime;++i)
          {
             for (const auto& ic : _hitIndex[i])
             {                
                 double dist = (init) ? distance2(*ic,chit) :  distance(*ic,chit);
                 if (dist < mindist) {mindist = dist; minc = ic;}
                 if (mindist < _dhit) break;               
             }          
             if (mindist < _dhit) break;               
          }

          // check if need to append hit form new cluster
          if (mindist < _dhit) 
          {
              minc->hits().emplace_back(&chit); 
          } 
          else if (mindist > _dseed)
          {
              clusters.emplace_back(ClusterStraw2(chit)); 
              minc = &clusters.back();
              _hitIndex[chit._itime].emplace_back(minc);
          } 
          else 
          {
               mindist =10000;
               minc = nullptr;
          }

          if (minc != nullptr)
          {             
               //associated to a new cluster, need to flag old/new  cluster accordingly
               if (_cptrs[chit._index] != minc) 
               { 
                  ++nchanged; 
                  ++chit._nchanged; 
                  if (_cptrs[chit._index] != nullptr) _cptrs[chit._index]->flagChanged(true); 
                  minc->flagChanged(true);
               }
               _cptrs[chit._index] = minc;
           } 
           else 
           {
               //unassociated to previous cluster, need to flag old cluster accordingly
               if (_cptrs[chit._index] != nullptr) 
               {
                  ++nchanged; 
                  ++chit._nchanged; 
                  _cptrs[chit._index]->flagChanged(true);
               }
               _cptrs[chit._index] = nullptr;
           }      
       }

       //clear cluster timing vector, update clusters and hit distance and refil timing vector
       for (auto& vec: _hitIndex) vec.clear(); 

       for(auto& cluster : clusters)
       {
           if (cluster.hasChanged())
           {
              cluster.updateCache(_maxwt);
              cluster.flagChanged(false); 

              if (cluster.hits().size()==1) 
                cluster.hits().at(0)->_dist = 0; 
              else
                for (auto& hit : cluster.hits()) hit->_dist = distance(cluster,*hit);                      
           }

           _hitIndex[cluster.itime()].emplace_back(&cluster);
       }      

       return nchanged;
   }


  //-------------------------------------------------------------------------------------------------------------------
  void TNTBClusterer::dump(std::list<ClusterStraw2>& clusters)
  {
      int iclu(0);      
      for (auto& cluster: clusters)
      { 
          //std::cout<<"Cluster "<<iclu<<" "<<cluster.pos()<<" "<<cluster.time()<<"  - ";
          std::cout<<"Cluster "<<iclu<<" "<<cluster.hits().size()<<"  - ";
          for (auto& hit : cluster.hits()) std::cout<<hit->_index<<" ";
          std::cout<<std::endl;
          ++iclu;
      }
  }

}


// Note 1: The list allows you to take pointer of its elements, a vector will reallocate and invalidate the pointers
//         The code with the vector version is much more complicated and has a very small performwnce improvement.


