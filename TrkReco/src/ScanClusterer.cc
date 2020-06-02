#include "TrkReco/inc/ScanClusterer.hh"
#include <vector>
#include <queue>
#include <algorithm>

namespace mu2e
{
  
   ScanClusterer::ScanClusterer(const Config& config) :
      tbin_       (config.tbin()),
      pbin_       (config.pbin()),     
      rbin_       (config.rbin()),     
      rmin_       (config.rmin()),     
      rmax_       (config.rmax()),     
      minPeakHit_ (config.minPeakHit()),     
      minSeedHit_ (config.minSeedHit()),     
      minRadHit_  (config.minRadHit()),     
      filterAlgo_ (config.filterAlgo()),
      bkgmask_    (config.bkgmsk()),
      sigmask_    (config.sigmsk()),
      testflag_   (config.testflag()),
      diag_       (config.diag())
   {
   }


        
   //---------------------------------------------------------------------------------------
   void ScanClusterer::init(){}


   //----------------------------------------------------------------------------------------------------------
   void ScanClusterer::findClusters(BkgClusterCollection& preFilterClusters, BkgClusterCollection& postFilterClusters,
                                      const ComboHitCollection& chcol, float mbtime, int iev)
   {             
       switch ( filterAlgo_ )
       {
         case 1:
            fastFilter1(preFilterClusters,chcol,mbtime);
            break;
         case 2:
            fastFilter2(preFilterClusters,chcol,mbtime);
            break;
        default:
            throw cet::exception("RECO")<< "[ScanClusterer] Unknown filter option " << filterAlgo_<< std::endl;
       }
   }
   
   //----------------------------------------------------------------------------------------------------------------------
   void ScanClusterer::fastFilter1(BkgClusterCollection& clusters, const ComboHitCollection& chcol, const float mbtime)
   {                            
       const unsigned nTimeBins = unsigned(mbtime/tbin_)+1;
       const unsigned nPhiBins  = unsigned(2*M_PI/pbin_)+1;
       const unsigned nRadBins  = unsigned((rmax_-rmin_)/rbin_)+1;
       const unsigned nRTBins   = nTimeBins*nRadBins;
       const unsigned nTotBins  = nTimeBins*nRadBins*nPhiBins;
    
    
       std::vector<short unsigned> timePhiHist(nTotBins,0);
       for (unsigned ich=0; ich<chcol.size();++ich)
       {
           const ComboHit& hit = chcol[ich];          
           if (testflag_ && (!hit.flag().hasAllProperties(sigmask_) || hit.flag().hasAnyProperty(bkgmask_))) continue;            

           unsigned pOffset = unsigned( (hit.phi()+M_PI)/pbin_);
           unsigned tOffset = unsigned( hit.time()/tbin_);
           unsigned rOffset = unsigned( (sqrtf(hit.pos().perp2())-rmin_)/rbin_);
           unsigned idx     = tOffset + rOffset*nTimeBins + pOffset*nRTBins;
           timePhiHist[idx] += 1;
       }
 

       std::vector<short unsigned> blindIdx(nTotBins,0);
       for (unsigned i=1;i<nTotBins-1;++i)
       {
           if (timePhiHist[i]<minPeakHit_) continue;
                  
           blindIdx[i]=1;

           unsigned tidx = i%nTimeBins;
           if (tidx>0)           blindIdx[i+1]=1;
           if (tidx+1<nTimeBins) blindIdx[i-1]=1;

           unsigned ridx = (i%nRTBins)/nTimeBins;
           if (ridx>0)          blindIdx[i-nTimeBins]=1;
           if (ridx+1<nRadBins) blindIdx[i+nTimeBins]=1;

           blindIdx[(i+nRTBins+nTotBins)%nTotBins]=1;
           blindIdx[(i-nRTBins+nTotBins)%nTotBins]=1;
       }


       //collect all preFiltered hits in a single cluster
       clusters.emplace_back(BkgCluster(XYZVec(0,0,0), 0));
       for (unsigned ich=0; ich<chcol.size();++ich)
       {
           const ComboHit& hit = chcol[ich];          
           unsigned pOffset = unsigned( (hit.phi()+M_PI)/pbin_);
           unsigned tOffset = unsigned( hit.time()/tbin_);
           unsigned rOffset = unsigned( (sqrtf(hit.pos().perp2())-rmin_)/rbin_);
           unsigned idx     = tOffset + rOffset*nTimeBins + pOffset*nRTBins;
           
//           if (timePhiHist[idx]>minPeakHit_) clusters.back().addHit(ich);
           if (blindIdx[idx]>0) clusters.back().addHit(ich);
       }        
   }


   //----------------------------------------------------------------------------------------------------------------------
   void ScanClusterer::fastFilter2(BkgClusterCollection& clusters, const ComboHitCollection& chcol, const float mbtime)
   {                            
       const unsigned nTimeBins = unsigned(mbtime/tbin_)+1;
       const unsigned nPhiBins  = unsigned(2*M_PI/pbin_)+1;
       const unsigned nRadBins  = unsigned((rmax_-rmin_)/rbin_)+1;
       const unsigned nRTBins   = nTimeBins*nRadBins;
       const unsigned nTotBins  = nTimeBins*nRadBins*nPhiBins;
    
       std::vector<short unsigned> timePhiHist(nTotBins,0), blindIdx(nTotBins,0);

       for (unsigned ich=0; ich<chcol.size();++ich)
       {
           const ComboHit& hit = chcol[ich];          
           if (testflag_ && (!hit.flag().hasAllProperties(sigmask_) || hit.flag().hasAnyProperty(bkgmask_))) continue;            

           unsigned pOffset = unsigned( (hit.phi()+M_PI)/pbin_);
           unsigned tOffset = unsigned( hit.time()/tbin_);
           unsigned rOffset = unsigned( (sqrtf(hit.pos().perp2())-rmin_)/rbin_);
           unsigned idx     = tOffset + rOffset*nTimeBins + pOffset*nRTBins;
           timePhiHist[idx] += 1;
       }
 
       int numCluster(1);
       for (unsigned idx=0;idx<nTotBins-1;++idx)
       {
           if (timePhiHist[idx]<minSeedHit_) continue;

           std::queue<unsigned> toProcess;
           toProcess.emplace(idx);

           while (!toProcess.empty())
           {          
              unsigned j = toProcess.front();                          
              if (timePhiHist[j] >= minPeakHit_)
              {
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
           unsigned tOffset = unsigned( hit.time()/tbin_);
           unsigned rOffset = unsigned( (sqrtf(hit.pos().perp2())-rmin_)/rbin_);
           unsigned idx     = tOffset + rOffset*nTimeBins + pOffset*nRTBins;
           
           if (blindIdx[idx]>0) clusters.back().addHit(ich);
       }        
   }


}



