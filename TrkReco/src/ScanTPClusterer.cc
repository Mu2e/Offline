#include "TrkReco/inc/ScanTPClusterer.hh"
#include <vector>
#include <queue>
#include <algorithm>

namespace mu2e
{
  
   ScanTPClusterer::ScanTPClusterer(const Config& config) :
      tmin_       (config.tmin()),
      tmax_       (config.tmax()),
      tbin_       (config.tbin()),
      pbin_       (config.pbin()),     
      minPeakHit_ (config.minPeakHit()),     
      minSeedHit_ (config.minSeedHit()),     
      minRadHit_  (config.minRadHit()),     
      dhit_       (config.hitDistance()),           
      dseed_      (config.seedDistance()),           
      dd_         (config.clusterDiameter()),   
      dt_         (config.clusterTime()),    
      tbinMinCluster_(config.deltaTimeBinMin()),    
      maxHitdt_   (config.maxHitTimeDiff()),      
      filterAlgo_ (config.filterAlgo()),         
      comboInit_  (config.comboInit()),         
      bkgmask_    (config.bkgmsk()),
      sigmask_    (config.sigmsk()),
      testflag_   (config.testflag()),
      diag_       (config.diag())
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
   void ScanTPClusterer::init(){}


   //----------------------------------------------------------------------------------------------------------
   void ScanTPClusterer::findClusters(BkgClusterCollection& preFilterClusters, BkgClusterCollection& postFilterClusters,
                                      const ComboHitCollection& chcol, float mbtime, int iev)
   {             
        std::vector<unsigned> hitSel; 
        hitSel.reserve(chcol.size());
                
        if (filterAlgo_==1) fastFilter1(preFilterClusters,chcol,hitSel);
        if (filterAlgo_==2) fastFilter2(preFilterClusters,chcol,hitSel);

        if (comboInit_) std::sort(hitSel.begin(),hitSel.end(),
                                 [&chcol](const auto& i, const auto& j) 
                                 {return chcol[i].wireRes() < chcol[j].wireRes();});      
         
        if (addOnePassClu_) fastCluster(postFilterClusters, chcol, hitSel, mbtime);
   }
   
   
   //----------------------------------------------------------------------------------------------------------------------
   void ScanTPClusterer::fastFilter1(BkgClusterCollection& clusters, const ComboHitCollection& chcol, std::vector<unsigned>& hitSel)
   {                            
       //NEED TO ADJUST THE NUMBER OF TIME BINS FOR OFF-SPILL
              
       const unsigned nTimeBins = unsigned((tmax_-tmin_)/tbin_);
       const unsigned nPhiBins  = unsigned(2*M_PI/pbin_+1e-5)+1;
       const unsigned nTotBins  = nTimeBins*nPhiBins;

       std::vector<unsigned> timePhiHist(nTotBins,0), blindIdx(nTotBins,0);
       std::vector<float>    radHist(nTotBins,0);      
       for (unsigned ich=0; ich<chcol.size();++ich)
       {
           const ComboHit& hit = chcol[ich];          
           if (testflag_ && (!hit.flag().hasAllProperties(sigmask_) || hit.flag().hasAnyProperty(bkgmask_))) continue;            

           unsigned pOffset = unsigned( (hit.phi()+M_PI)/pbin_);
           unsigned tOffset = unsigned( (hit.time()-tmin_)/tbin_);
           unsigned idx     =  tOffset + pOffset*nTimeBins;
           timePhiHist[idx] += 1;  //we could make it nStrawHits here instead
           radHist[idx] += sqrtf(hit.pos().perp2());
       }

       unsigned numCluster(1);
       std::vector<float> radii{0};
       for (unsigned idx=0;idx<nTotBins-1;++idx)
       {          
           unsigned numHitBox(timePhiHist[idx]);
           if (idx%nTimeBins>0) numHitBox += timePhiHist[idx-1];
           if (idx%nTimeBins+1<nTimeBins) numHitBox += timePhiHist[idx+1];
           if (numHitBox<minSeedHit_) continue;
           
           //if (timePhiHist[idx]<minSeedHit_) continue;
          

           std::queue<unsigned> toProcess;
           toProcess.emplace(idx);

           float sumRad(0),sumNorm(0);
           while (!toProcess.empty())
           {             
              unsigned j = toProcess.front();                          
              if (timePhiHist[j] >= minPeakHit_)
              {
                  blindIdx[j] = numCluster;
                  sumRad  += radHist[j];
                  sumNorm += timePhiHist[j];    

                  if (j%nTimeBins>0)           toProcess.emplace(j-1);
                  if (j%nTimeBins+1<nTimeBins) toProcess.emplace(j+1);
                  toProcess.emplace((j+nTimeBins+nTotBins)%nTotBins);
                  toProcess.emplace((j-nTimeBins+nTotBins)%nTotBins);
              } 

              timePhiHist[j]=0;                     
              toProcess.pop();
           }
           radii.push_back(sumRad/sumNorm);                    
           ++numCluster;           
       }

       //collect all preFiltered hits in a single cluster
       clusters.emplace_back(BkgCluster(XYZVec(0,0,0), 0));
       for (unsigned ich=0; ich<chcol.size();++ich)
       {
           const ComboHit& hit = chcol[ich];          
           unsigned pOffset = unsigned( (hit.phi()+M_PI)/pbin_);
           unsigned tOffset = unsigned( (hit.time()-tmin_)/tbin_);
           unsigned idx     =  tOffset + pOffset*nTimeBins;
           
           if (blindIdx[idx]==0) {hitSel.emplace_back(ich); continue;}
           
           float dr = std::abs(sqrtf(hit.pos().perp2())-radii[blindIdx[idx]]);
           
           if (dr > minRadHit_) hitSel.emplace_back(ich); 
           else                 clusters.back().addHit(ich);      
       }        
   }







   //----------------------------------------------------------------------------------------------------------------------
   void ScanTPClusterer::fastFilter2(BkgClusterCollection& clusters, const ComboHitCollection& chcol, std::vector<unsigned>& hitSel)
   {                            
       const float rmax_ = 650;
       const float rmin_ = 350;
       const float rbin_ = 50;
     
       const unsigned nTimeBins = unsigned((tmax_-tmin_)/tbin_);
       const unsigned nPhiBins  = unsigned(2*M_PI/pbin_+1e-5)+1;
       const unsigned nRadBins  = unsigned((rmax_-rmin_)/rbin_);
       const unsigned nRTBins   = nTimeBins*nRadBins;
       const unsigned nTotBins  = nTimeBins*nRadBins*nPhiBins;
    
       std::vector<short unsigned> timePhiHist(nTotBins,0), blindIdx(nTotBins,0);

       for (unsigned ich=0; ich<chcol.size();++ich)
       {
           const ComboHit& hit = chcol[ich];          
           if (testflag_ && (!hit.flag().hasAllProperties(sigmask_) || hit.flag().hasAnyProperty(bkgmask_))) continue;            

           unsigned pOffset = unsigned( (hit.phi()+M_PI)/pbin_);
           unsigned tOffset = unsigned( (hit.time()-tmin_)/tbin_);
           unsigned rOffset = unsigned( (sqrtf(hit.pos().perp2())-rmin_)/rbin_);
           unsigned idx     = tOffset + rOffset*nTimeBins + pOffset*nRTBins;
           timePhiHist[idx] += 1;
       }
 
       int numCluster(1);
       for (unsigned idx=0;idx<nTotBins-1;++idx){

           if (timePhiHist[idx]<minSeedHit_) continue;

           std::queue<unsigned> toProcess;
           toProcess.emplace(idx);

           while (!toProcess.empty()){          

              unsigned j = toProcess.front();                          
              if (timePhiHist[j] >= minPeakHit_) {

                  blindIdx[j] = numCluster;

                  //look at neighboring time cells
                  unsigned tidx = j%nTimeBins;
                  if (tidx>0)           toProcess.emplace(j-1);
                  if (tidx+1<nTimeBins) toProcess.emplace(j+1);

                  //neighboring in radius cells
                  unsigned ridx = (j%nRTBins)/nTimeBins;
                  if (ridx>0)          toProcess.emplace(j-nTimeBins);
                  if (ridx+1<nRadBins) toProcess.emplace(j+nTimeBins);

                  //neighboring in phi cells with cyclic couner
                  toProcess.emplace((j+nRTBins+nTotBins)%nTotBins);
                  toProcess.emplace((j-nRTBins+nTotBins)%nTotBins);
             } 

              timePhiHist[j]=0;                     
              toProcess.pop();
           }
           ++numCluster;           
       }

       //collect all preFiltered hits in a single cluster
       clusters.emplace_back(BkgCluster(XYZVec(0,0,0), 0));
       for (unsigned ich=0; ich<chcol.size();++ich)
       {
           const ComboHit& hit = chcol[ich];          
           unsigned pOffset = unsigned( (hit.phi()+M_PI)/pbin_);
           unsigned tOffset = unsigned( (hit.time()-tmin_)/tbin_);
           unsigned rOffset = unsigned( (sqrtf(hit.pos().perp2())-rmin_)/rbin_);
           unsigned idx     = tOffset + rOffset*nTimeBins + pOffset*nRTBins;
           
           if (blindIdx[idx]>0) clusters.back().addHit(ich);
           else                 hitSel.emplace_back(ich); 
        }        
   }







       
   //----------------------------------------------------------------------------------------------------------------------
   void ScanTPClusterer::fastCluster(BkgClusterCollection& clusters, const ComboHitCollection& chcol, const std::vector<unsigned>& hitSel, float mbtime)
   {
       //adjust the time binning to index clusters in the clustering algo
       float tbin = std::max(mbtime/float(numBuckets)+0.001f,tbinMinCluster_);
       if (int(mbtime/tbin) >= numBuckets) throw cet::exception("RECO")<< "Too many bucket bins requested for TNTCluster!"<< std::endl;        
       
       //calculate the lookup indices for the clustering algorithm
       int ditime = int(maxHitdt_/tbin);
       std::vector<int> hitDtIdx;
       for (int i=0;i<=ditime;++i) {hitDtIdx.push_back(i); if (i>0) hitDtIdx.push_back(-i);}                 

       //structure to hold the cluster indices vs time bucket
       std::array<std::vector<int>,numBuckets> hitIndex;
       for (auto& vec : hitIndex) vec.reserve(16);
       
       // -- then do one pass clustering
       for (auto ihit : hitSel)
       {                        
           // -- Find cluster closest to hit
           int minc(-1);                  
           float mindist(dseed_+1.0f);         
           int itime = int(chcol[ihit].time()/tbin);

           for (auto i : hitDtIdx)
           {
               for (const auto& ic : hitIndex[itime+i])
               {                
                   float dist = distance(clusters[ic],chcol[ihit]);
                   if (dist < mindist) {mindist = dist; minc = ic;}
                   if (mindist < dhit_) break;               
               }          
               if (mindist < dhit_) break;               
           }

           // either form new cluster, add hit to existing cluster or do nothing (two level threshold)
           if (mindist < dhit_) clusters[minc].addHit(ihit);
           if (mindist > dseed_)
           {
               clusters.emplace_back(BkgCluster(chcol[ihit].pos(),chcol[ihit].time())); 
               clusters.back().addHit(ihit);         
               
	       int itimeClu  = int(chcol[ihit].time()/tbin);
               hitIndex[itimeClu].emplace_back(clusters.size()-1);
           } 
       }
   }
   
   
   //---------------------------------------------------------------------------------------
   // only count differences if they are above the natural hit size (drift time, straw size)      
   float ScanTPClusterer::distance(const BkgCluster& cluster, const ComboHit& hit) const 
   {     
       float psep_x = hit.pos().x()-cluster.pos().x();
       float psep_y = hit.pos().y()-cluster.pos().y();
       float d2     = psep_x*psep_x+psep_y*psep_y;

       if (d2 > md2_) return dseed_+1.0f; 

       float dt = std::abs(hit.time()-cluster.time());
       if (dt > maxHitdt_) return dseed_+1.0f;


       float retval(0.0);
       if (dt > dt_) {float tdist = dt -dt_;  retval = tdist*tdist*trms2inv_;}      
       if (d2 > dd2_) 
       {	
           //This is equivalent to but faster than these lines - yes, that's nice
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

   
}


/*
   //------------------------------------------------------------------------------------------------------------------------
   // The simplest clustering algorithm in phi/time
   void ScanTPClusterer::cluAlgo1(const ComboHitCollection& chcol) 
   {      
        //FIND PEAK IN Phi/TIME DISTRIBUTION
        const unsigned minPeakHit_(6U);
        const unsigned minSumHit_(12U);

        const unsigned nTimeBins = unsigned((tmax_-tmin_)/tbin_)+2;
        const unsigned nPhiBins  = unsigned(2*M_PI/pbin_+1e-5)+1;
        const unsigned nTotBins  = nTimeBins*nPhiBins;
        std::vector<unsigned> timePhiHist(nTotBins,0);

        for (unsigned ich=0; ich<chcol.size();++ich)
        {
            const ComboHit& hit = chcol[ich];          
            if (testflag_ && (!hit.flag().hasAllProperties(sigmask_) || hit.flag().hasAnyProperty(bkgmask_))) continue;            
            unsigned pOffset = unsigned( (hit.phi()+M_PI)/pbin_);
            unsigned tOffset = unsigned( (hit.time()-tmin_)/tbin_)+1;
            unsigned idx     = tOffset + pOffset*nTimeBins;
            timePhiHist[idx] += 1;
        }

        for (unsigned idx=0;idx<nTotBins-1;++idx)
        {
            if (timePhiHist[idx]<minPeakHit_) continue;

            unsigned idxUp   = (idx+nTimeBins+nTotBins)%nTotBins;
            unsigned idxDown = (idx-nTimeBins+nTotBins)%nTotBins;
            unsigned sum     = timePhiHist[idx]+timePhiHist[idx-1]+timePhiHist[idx+1]+timePhiHist[idxUp]+timePhiHist[idxUp-1]+
                               timePhiHist[idxUp+1]+timePhiHist[idxDown]+timePhiHist[idxDown-1]+timePhiHist[idxDown+1];

            if (sum >= minSumHit_ && timePhiHist[idx]>=minPeakHit_ ) 
            {
               blindIdx[idx]     = 1; blindIdx[idx+1]     = 1; blindIdx[idx-1]     = 1;
               blindIdx[idxUp]   = 1; blindIdx[idxUp+1]   = 1; blindIdx[idxUp-1]   = 1;
               blindIdx[idxDown] = 1; blindIdx[idxDown+1] = 1; blindIdx[idxDown-1] = 1;
            }   
        }
*/




/*      

      // THIS IS A VERSION IN R,TIME,PHI SPACE. IT IS A BIT LESS EFFICIENT THAN THE PHI/TIME VERSION 
      // SINCE THE VECTOR TO STORE THE RESULTS ARE QUITE LARGE

      const float rmax_ = 650;
      const float rmin_ = 350;
      const float rbin_ = 50;
     
      const unsigned nTimeBins = unsigned((tmax_-tmin_)/tbin_);
      const unsigned nPhiBins  = unsigned(2*M_PI/pbin_)+1;
      const unsigned nRadBins  = unsigned((rmax_-rmin_)/rbin_);
      const unsigned nRTBins   = nTimeBins*nRadBins;
      const unsigned nTotBins  = nTimeBins*nRadBins*nPhiBins;
    
      std::vector<short unsigned> timePhiHist(nTotBins,0), blindIdx(nTotBins,0);

      for (unsigned ich=0; ich<nch;++ich)
      {
          const ComboHit& hit = chcol[ich];          
          if (testflag_ && (!hit.flag().hasAllProperties(sigmask_) || hit.flag().hasAnyProperty(bkgmask_))) continue;            
          
          unsigned pOffset = unsigned( (hit.phi()+M_PI)/pbin_);
          unsigned tOffset = unsigned( (hit.time()-tmin_)/tbin_);
          unsigned rOffset = unsigned( (sqrtf(hit.pos().perp2())-rmin_)/rbin_);
          unsigned idx     = tOffset + rOffset*nTimeBins + pOffset*nRTBins;
          timePhiHist[idx] += 1;
      }

      int numCluster(1);
      for (unsigned idx=0;idx<nTotBins-1;++idx){
          
          if (timePhiHist[idx]<minSeedHit_) continue;
          
          std::queue<unsigned> toProcess;
          toProcess.emplace(idx);
          
          while (!toProcess.empty()){          
             
             unsigned j = toProcess.front();                          
             if (timePhiHist[j] >= minPeakHit_) {
                 
                 blindIdx[j] = numCluster;

                 //look at neighboring time cells
                 unsigned tidx = j%nTimeBins;
                 if (tidx>0)           toProcess.emplace(j-1);
                 if (tidx+1<nTimeBins) toProcess.emplace(j+1);

                 //neighboring in radius cells
                 unsigned ridx = (j%nRTBins)/nTimeBins;
                 if (ridx>0)          toProcess.emplace(j-nTimeBins);
                 if (ridx+1<nRadBins) toProcess.emplace(j+nTimeBins);
                
                 //neighboring in phi cells with cyclic couner
                 toProcess.emplace((j+nRTBins+nTotBins)%nTotBins);
                 toProcess.emplace((j-nRTBins+nTotBins)%nTotBins);
            } 

             timePhiHist[j]=0;                     
             toProcess.pop();
          }
          ++numCluster;           
      }
*/






/*
       //ANOTHER VERSION OF THE CLUSTERING USING A LINEARIZED VECTOR, MIGHT BE BETTER FOR FPGA
       
       //adjust the time binning to index clusters in the clustering algo
       const unsigned vecLen(2048);
       float tbin   = std::max(mbtime/float(veclen);
       int   ditime = int(maxHitdt_/tbin);

       std::vector<int> cluIdx(vecLen,-1);
       for (size_t ihit=0; ihit<nch; ++ihit)
       {                        
           unsigned pOffset = unsigned( (chcol[ihit].phi()+M_PI)/pbin_);
           unsigned tOffset = unsigned( (chcol[ihit].time()-tmin_)/tbin_)+1;
           unsigned idx     =  pOffset*nTimeBins + tOffset;
           if (blindIdx[idx]>0) continue;

           // -- Find cluster closest to hit
           int minc(-1);                  
           float mindist(dseed_+1.0f);         
           int itime = int(chcol[ihit].time());

           for (int ic=0;ic<ditime;++ic)
           {
               if (cluIdx[itime+ic]>-1) {
                 float dist = distance(bkgccol[cluIdx[itime+ic]],chcol[ihit]);
                 if (dist < mindist) {mindist = dist; minc = itime+ic;}
                 if (mindist < dhit_) break;
               }               
               if (cluIdx[itime-ic]>-1) {
                 float dist = distance(bkgccol[cluIdx[itime-ic]],chcol[ihit]);
                 if (dist < mindist) {mindist = dist; minc = itime-ic;}
                 if (mindist < dhit_) break;
               }               
           }
           // -- Form new cluster, add hit to new cluster or do nothing
           if (mindist < dhit_) bkgccol[cluIdx[minc]].addHit(ihit);
           if (mindist > dseed_)
           {
               bkgccol.emplace_back(BkgCluster(chcol[ihit].pos(),chcol[ihit].time())); 
               bkgccol.back().addHit(ihit);         
               
	       int itimeClu  = int(chcol[ihit].time());
               if (cluIdx[itimeClu]==-1) cluIdx[itimeClu]=bkgccol.size()-1;
           } 
       }
*/
