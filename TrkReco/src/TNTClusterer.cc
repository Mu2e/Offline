#include "TrkReco/inc/TNTClusterer.hh"
#include <vector>
#include <algorithm>

namespace mu2e
{
   TNTClusterer::TNTClusterer(const Config& config) :
      hitDtIdx_(),
      dhit_      (config.hitDistance()),           
      dseed_     (config.seedDistance()),           
      dd_        (config.clusterDiameter()),   
      dt_        (config.clusterTime()),    
      tbinMin_   (config.deltaTimeBinMin()),    
      maxHitdt_  (config.maxHitTimeDiff()),      
      maxDistSum_(config.maxSumDistance()),   
      maxNiter_  (config.maxCluIterations()),
      useMedian_ (config.medianCentroid()),         
      comboInit_ (config.comboInit()),         
      bkgmask_   (config.bkgmsk()),
      sigmask_   (config.sigmsk()),
      testflag_  (config.testflag()),
      diag_      (config.diag())
   {
       // cache some values
       float minerr (config.minHitError());
       float maxdist(config.maxDistance());
       float trms   (config.timeRMS());            

       dd2_      = dd_*dd_;
       maxwt_    = 1.0f/minerr;
       md2_      = maxdist*maxdist;
       trms2inv_ = 1.0f/trms/trms;                 
   }

        
   //---------------------------------------------------------------------------------------
   void TNTClusterer::init(){}


   //----------------------------------------------------------------------------------------------------------
   void TNTClusterer::findClusters(BkgClusterCollection& clusters, const ComboHitCollection& chcol, float mbtime, int iev)
   {        
        std::vector<BkgHit> BkgHits;
        BkgHits.reserve(chcol.size());

        //adjust the time binning to index clusters in the clustering algo
        float tbin = std::max(mbtime/float(numBuckets)+0.001f,tbinMin_);
        if (int(mbtime/tbin) >= numBuckets) throw cet::exception("RECO")<< "Too many bucket bins requested for TNTCluster!"<< std::endl;        
        int  ditime = int(maxHitdt_/tbin);
        hitDtIdx_.clear();
        for (int i=0;i<=ditime;++i) {hitDtIdx_.push_back(i); if (i>0) hitDtIdx_.push_back(-i);}                 

        initClu(chcol, clusters, BkgHits);
        clusterAlgo(chcol, clusters, BkgHits, tbin);

        //a final merge does almost nothing but we leave it in case of need
        //mergeClusters(clusters, chcol, dt_, dd2_);

        clusters.erase(std::remove_if(clusters.begin(),clusters.end(),[](auto& cluster){return cluster.hits().empty();}),clusters.end());

        //Transform BkgHits indices into ComboHit indices
        for (auto& cluster: clusters)
            std::transform(cluster.hits().begin(),cluster.hits().end(),cluster.hits().begin(),
                           [&BkgHits] (const int i){return BkgHits[i].chidx_;});
   }
   

   //------------------------------------------------------------------------------------------------------------------------
   void TNTClusterer::initClu(const ComboHitCollection& chcol, std::vector<BkgCluster>& clusters, std::vector<BkgHit>& BkgHits) 
   {      
        for (size_t ish=0; ish<chcol.size(); ++ish)
        {         
             if (testflag_ && (!chcol[ish].flag().hasAllProperties(sigmask_) || chcol[ish].flag().hasAnyProperty(bkgmask_))) continue;           
             BkgHits.emplace_back(BkgHit(ish));                 
        }

        if (comboInit_) std::sort(BkgHits.begin(),BkgHits.end(),
                                  [&chcol](const BkgHit& x, const BkgHit& y)
                                  {return chcol[x.chidx_].wireRes() < chcol[y.chidx_].wireRes();});      
   }  
   
   
   //----------------------------------------------------------------------------------------------------------------------
   void TNTClusterer::clusterAlgo(const ComboHitCollection& chcol, std::vector<BkgCluster>& clusters, std::vector<BkgHit>& BkgHits, float tbin)
   {                            
        arrayVecBkg hitIndex;
        for (auto& vec : hitIndex) vec.reserve(16);
              
        unsigned niter(0);
        float odist(2.0f*maxDistSum_),tdist(0.0f); 
        while (std::abs(odist - tdist) > maxDistSum_ && niter < maxNiter_)
        {        
	    ++niter;
	    formClusters(chcol, clusters,tbin,  BkgHits, hitIndex);

            odist = tdist;      
            tdist = 0.0f;
            for (const auto& cluster: clusters)
                for (const auto& cidx : cluster.hits()) tdist += BkgHits[cidx].distance_;           
        }           
   }

   
   //-------------------------------------------------------------------------------------------------------------------
   // loop over hits, re-affect them to their original cluster if they are still within the radius, otherwise look at 
   // candidate clusters to check if they could be added. If not, make a new cluster.
   // to speed up, do not update clusters who haven't changed and keep a list of clusters within a given time window
   //
   unsigned TNTClusterer::formClusters(const ComboHitCollection& chcol, std::vector<BkgCluster>& clusters, float tbin,  
                                       std::vector<BkgHit>& BkgHits, arrayVecBkg& hitIndex)
   {                         
       unsigned nchanged(0);       
       for (auto& cluster : clusters) cluster.clearHits();

       for (size_t ihit=0; ihit<BkgHits.size(); ++ihit)
       {                        
           // -- if hit is ok, reassign it right away
           auto& hit = BkgHits[ihit];                    
           if (hit.distance_ < dhit_) 
           {               
               clusters[hit.clusterIdx_].addHit(ihit); 
               continue;
           }

           // -- Find cluster closest to hit
           int minc(-1);                  
           float mindist(dseed_+1.0f);         
           int itime = int(chcol[hit.chidx_].time()/tbin);

           for (auto i : hitDtIdx_)
           {
               for (const auto& ic : hitIndex[itime+i])
               {                
                   float dist = distance(clusters[ic],chcol[hit.chidx_]);
                   if (dist < mindist) {mindist = dist; minc = ic;}
                   if (mindist < dhit_) break;               
               }          
               if (mindist < dhit_) break;               
           }

           // -- Form new cluster, add hit to new cluster or do nothing
           if (mindist < dhit_) 
           {
               clusters[minc].addHit(ihit);
           }
           else if (mindist > dseed_)
           {
               minc = clusters.size();
               clusters.emplace_back(BkgCluster(chcol[hit.chidx_].pos(),chcol[hit.chidx_].time())); 
               clusters[minc].addHit(ihit);         
               int itimeClu  = int(chcol[hit.chidx_].time()/tbin);
               hitIndex[itimeClu].emplace_back(minc);
           } 
           else 
           {
                BkgHits[ihit].distance_ = 10000.0f;
                minc = -1;
           }

           // -- Update cluster flag and hit->cluster pointer if associated to new cluster or removed from previous cluster
           if (minc != -1)
           {             
               if (hit.clusterIdx_ != minc)
               { 
                   ++nchanged; 
                   if (hit.clusterIdx_ != -1) clusters[hit.clusterIdx_]._flag = BkgClusterFlag::update; 
                   clusters[minc]._flag = BkgClusterFlag::update;
               }
               hit.clusterIdx_ = minc;
            } 
            else 
            {
               if (hit.clusterIdx_ != -1)
               {
                   ++nchanged; 
                   clusters[hit.clusterIdx_]._flag = BkgClusterFlag::update;
               }
               hit.clusterIdx_ = -1;
            }      
       }


       //update cluster, hit distance and maps
       for (auto& vec: hitIndex) vec.clear(); 

       for (unsigned ic=0;ic<clusters.size();++ic)
       {
           BkgCluster& cluster = clusters[ic];
           if (cluster._flag == BkgClusterFlag::update) 
           {
              cluster._flag = BkgClusterFlag::unchanged; 
              updateCluster(cluster, chcol, BkgHits); 

              if (cluster.hits().size()==1)       
                  BkgHits[cluster.hits().at(0)].distance_ = 0.0f; 
              else 
                  for (auto& hit : cluster.hits()) BkgHits[hit].distance_ = distance(cluster,chcol[BkgHits[hit].chidx_]);                                                                

           }

           int itimeClu  = int(cluster.time()/tbin);
           hitIndex[itimeClu].emplace_back(ic);             
       }

       return nchanged;
   }


   //-----------------------------------------------------------------------------------------------
   void TNTClusterer::mergeClusters(std::vector<BkgCluster>& clusters, const ComboHitCollection& chcol, 
                                    std::vector<BkgHit>& BkgHits, float dt, float dd2)
   {
        unsigned niter(0);    
        while (niter < maxNiter_)
        {
            int nchanged(0);
            for (auto it1 = clusters.begin(); it1 != std::prev(clusters.end()); ++it1)
            {
                 if (it1->hits().empty()) continue;
                 for (auto it2 = std::next(it1); it2!=clusters.end(); ++it2)
                 {
                     if (it2->hits().empty()) continue;
                     if (std::abs(it1->time() - it2->time()) > dt) continue;
		     if ((it1->pos() - it2->pos()).perp2() > dd2)  continue;

                     ++nchanged;
                     mergeTwoClu(*it1,*it2);                     
	         }             
             }	

             ++niter;
             if (diag_>0) std::cout<<"Merge "<<niter<<" "<<nchanged<<"  "<<clusters.size()<<std::endl;

             if (nchanged==0) break;

             std::remove_if(clusters.begin(),clusters.end(),[](auto& cluster){return cluster.hits().empty();});
             for (auto& cluster : clusters ) updateCluster(cluster, chcol, BkgHits);
        }
        return;    
   }

   void TNTClusterer::mergeTwoClu(BkgCluster& clu1, BkgCluster& clu2 )
   {
       if (clu1.hits().empty() || clu2.hits().empty()) return;
       clu1.hits().insert(clu1.hits().end(),clu2.hits().begin(), clu2.hits().end());
       clu2.hits().clear();
   }


   //---------------------------------------------------------------------------------------
   // only count differences if they are above the natural hit size (drift time, straw size)      
   float TNTClusterer::distance(const BkgCluster& cluster, const ComboHit& hit) const 
   {     
       float psep_x = hit.pos().x()-cluster.pos().x();
       float psep_y = hit.pos().y()-cluster.pos().y();
       float d2     = psep_x*psep_x+psep_y*psep_y;

       if (d2 > md2_) return dseed_+1.0f; 

       float dt = std::abs(hit.time()-cluster.time());
       if (dt > maxHitdt_) return dseed_+1.0f;


       float retval(0.0f);
       if (dt > dt_) {float tdist = dt -dt_;  retval = tdist*tdist*trms2inv_;}      
       if (d2 > dd2_) 
       {	
           //This is equivalent to but faster than - yes, that's nice
           //XYZVec that(-hit.wdir().y(),hit.wdir().x(),0.0);
           //float dw = std::max(0.0f,hit.wdir().Dot(psep)-dd_)/hit.posRes(ComboHit::wire);
	   //float dp = std::max(0.0f,that.Dot(psep)-dd_)*maxwt_;  //maxwt = 1/minerr
           float hwx = hit.wdir().x();
           float hwy = hit.wdir().y();
           float dw  = std::max(0.0f,(psep_x*hwx+psep_y*hwy-dd_)/hit.posRes(ComboHit::wire));
           float dp  = std::max(0.0f,(hwx*psep_y-hwy*psep_x-dd_)*maxwt_);
	   retval += dw*dw + dp*dp;
       }      
       return retval;
   }


   //-------------------------------------------------------------------------------------------------------------------
   void TNTClusterer::updateCluster(BkgCluster& cluster, const ComboHitCollection& chcol, std::vector<BkgHit>& BkgHits)
   {                    
       
       if (cluster.hits().empty()) {cluster.time(0.0f); cluster.pos(XYZVec(0.0f,0.0f,0.0f));return;}  
      
       if (cluster.hits().size()==1)
       {
           int idx = BkgHits[cluster.hits().at(0)].chidx_;          
           cluster.time(chcol[idx].time());
           cluster.pos(XYZVec(chcol[idx].pos().x(),chcol[idx].pos().y(),0.0f));
           return;
       } 
              
       float crho  = sqrtf(cluster.pos().perp2());
       float cphi  = cluster.pos().phi();
       float ctime = cluster.time();

       if (useMedian_) 
       {
           std::vector<float> racc,pacc,tacc;
           for (auto& hit : cluster.hits())
           {
              int idx  = BkgHits[hit].chidx_;
              float dt = chcol[idx].time() - ctime;
	      float dr = sqrtf(chcol[idx].pos().perp2()) - crho;          
              float dp = chcol[idx].phi() - cphi;
              if (dp > M_PI)  dp -= 2*M_PI;
              if (dp < -M_PI) dp += 2*M_PI;

              // weight according to the # of hits
              for (int i=0;i<chcol[idx].nStrawHits();++i)
              {
                 racc.emplace_back(dr); 
                 pacc.emplace_back(dp);
                 tacc.emplace_back(dt);
              }
           }

           size_t vecSize = racc.size()/2;
           std::nth_element(racc.begin(),racc.begin()+vecSize,racc.end());
           std::nth_element(pacc.begin(),pacc.begin()+vecSize,pacc.end());
           std::nth_element(tacc.begin(),tacc.begin()+vecSize,tacc.end());

           crho  += racc[vecSize];
           cphi  += pacc[vecSize];
           ctime += tacc[vecSize];
       }
       else
       {  
           float sumWeight(0), deltaT(0), deltaP(0), deltaR(0);
           for (auto& hit : cluster.hits())
           {
               int   idx    = BkgHits[hit].chidx_;          
               float weight = chcol[idx].nStrawHits();
               float dt     = chcol[idx].time()-ctime;
               float dr     = sqrtf(chcol[idx].pos().perp2()) - crho;

               float dp     = chcol[idx].phi()-cphi;
               if (dp > M_PI)  dp -= 2*M_PI;
               if (dp < -M_PI) dp += 2*M_PI;

               deltaT    += dt*weight;
	       deltaR    += dr*weight;
               deltaP    += dp*weight;
               sumWeight += weight;
           }
           crho  += deltaR/sumWeight;
           cphi  += deltaP/sumWeight;
           ctime += deltaT/sumWeight;
       }
       
       cluster.time(ctime);       
       cluster.pos(XYZVec(crho*cos(cphi),crho*sin(cphi),0.0f));
   }



   //-------------------------------------------------------------------------------------------
   void TNTClusterer::dump(const std::vector<BkgCluster>& clusters, std::vector<BkgHit>& BkgHits)
   {
       int iclu(0);      
       for (auto& cluster: clusters)
       { 
           std::cout<<"Cluster "<<iclu<<" "<<cluster.pos()<<" "<<cluster.time()<<"  "<<cluster.hits().size()<<"  - ";
           for (auto& hit : cluster.hits()) std::cout<<BkgHits[hit].chidx_<<" ";
           std::cout<<std::endl;
           ++iclu;
       }
   }
   
}


