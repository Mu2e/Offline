// $Id: FlagBkgHits2_module.cc, without diagnostics $
// $Author: brownd & mpettee $ 
// $Date: 2016/11/30 $
//

#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "GeometryService/inc/GeomHandle.hh"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"

#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StrawHitPosition.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/BkgCluster.hh"
#include "RecoDataProducts/inc/BkgQual.hh"
#include "RecoDataProducts/inc/StereoHit.hh"

#include "TrkReco/inc/TLTClusterer.hh"
#include "TrkReco/inc/TNTClusterer.hh"
#include "TrkReco/inc/TNTBClusterer.hh"
#include "Mu2eUtilities/inc/MVATools.hh"

#include "CLHEP/Units/PhysicalConstants.h"
#include "TMath.h"
#include "TH1F.h"
#include "TMVA/Reader.h"

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/weighted_variance.hpp> 

#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <functional>
#include <float.h>
#include <vector>
#include <set>
#include <map>
#include <numeric>

using namespace boost::accumulators;

namespace mu2e 
{

  class FlagBkgHits2 : public art::EDProducer
  {
    public:
       enum clusterer { TwoLevelThreshold=1, TwoNiveauThreshold=2, TwoNiveauThresholdB=3};
       explicit FlagBkgHits2(fhicl::ParameterSet const&);
       virtual ~FlagBkgHits2();
       virtual void beginJob();
       virtual void produce(art::Event& event ); 
    
    private:

      int _debug;
      int _printfreq;
      art::InputTag _shtag;
      art::InputTag _shptag;
      art::InputTag _shftag;
      art::InputTag _sthTag;
      bool          _flagall;
      unsigned _minnhits;
      unsigned _minnstereo;
      unsigned _minnp;
      unsigned _maxisolated;

      BkgClusterer* _clusterer;
      double        _cperr2;
      bool          _useMVA;
      double        _bkgMVAcut;
      MVATools      _bkgMVA;

      // use TMVA Reader until problems with MVATools is resolved
      mutable TMVA::Reader _reader;
      mutable std::vector<float> _mvavars;
      TH1F *_v1a,*_v2a,*_v3a,*_v4a,*_v5a,*_v6a,*_v7a,*_v8a,*_v9a;
      TH1F *_v1b,*_v2b,*_v3b,*_v4b,*_v5b,*_v6b,*_v7b,*_v8b,*_v9b;

      void classifyClusters(BkgClusterCollection& clusters,BkgQualCollection& cquals, 
                            const StrawHitCollection& shcol, const StrawHitPositionCollection& shpcol) const;
      void fillBkgQual(const BkgCluster& cluster,BkgQual& cqual, const TTracker& tt, 
                       const StrawHitCollection& shcol, const StrawHitPositionCollection& shpcol) const;
      void countHits(const BkgCluster& cluster, unsigned& nactive, unsigned& nstereo) const;
      void countPlanes(const BkgCluster& cluster,BkgQual& cqual, const TTracker& tt, const StrawHitCollection& shcol) const;
  };

  FlagBkgHits2::FlagBkgHits2(const fhicl::ParameterSet& pset) :
    _debug(pset.get<int>(           "debugLevel",0)),
    _printfreq(pset.get<int>(       "printFrequency",101)),
    _shtag(pset.get<art::InputTag>( "StrawHitCollectionLabel","makeSH")),
    _shptag(pset.get<art::InputTag>("StrawHitPositionCollectionLabel","MakeStereoHits")),
    _shftag(pset.get<art::InputTag>("StrawHitFlagCollectionLabel","MakeStereoHits")),
    _sthTag(pset.get<art::InputTag>("StereoHitCollection","MakeStereoHits")),
    _flagall(pset.get<bool>(        "FlagAllHits",true)), 
    _minnhits(pset.get<unsigned>(   "MinActiveHits",5)),
    _minnstereo(pset.get<unsigned>( "MinStereoHits",2)),
    _minnp(pset.get<unsigned>(      "MinNPlanes",4)),
    _maxisolated(pset.get<unsigned>("MaxIsolated",1)),
    _useMVA(pset.get<bool>(         "UseBkgMVA",true)),
    _bkgMVAcut(pset.get<double>(    "BkgMVACut",0.8)),
    _bkgMVA(pset.get<fhicl::ParameterSet>("BkgMVA",fhicl::ParameterSet()))
  {
      produces<StrawHitFlagCollection>();
      produces<BkgClusterCollection>();
      produces<BkgQualCollection>();
      
      clusterer ctype = static_cast<clusterer>(pset.get<int>("Clusterer",TwoLevelThreshold));
      switch ( ctype ) {
        case TwoLevelThreshold:
	  _clusterer = new TLTClusterer(pset.get<fhicl::ParameterSet>("TLTClusterer",fhicl::ParameterSet()));
	  break;
         case TwoNiveauThreshold:
	  _clusterer = new TNTClusterer(pset.get<fhicl::ParameterSet>("TNTClusterer",fhicl::ParameterSet()));
	  break;
         case TwoNiveauThresholdB:
	  _clusterer = new TNTBClusterer(pset.get<fhicl::ParameterSet>("TNTBClusterer",fhicl::ParameterSet()));
	  break;
       default:
	  throw cet::exception("RECO")<< "Unknown clusterer" << ctype << std::endl;
      }
    
      double cperr = pset.get<double>("ClusterPositionError",10.0); 
      _cperr2 = cperr*cperr;

      std::vector<std::string> varnames(pset.get<std::vector<std::string> >("MVANames"));
      _mvavars.reserve(varnames.size());
      for (const auto& vname : varnames)
      {
         _mvavars.push_back(0.0);
         _reader.AddVariable(vname,&_mvavars.back());
      }

      ConfigFileLookupPolicy configFile;
      std:: string weights = configFile(pset.get<std::string>("BkgMVA.MVAWeights"));
      _reader.BookMVA("MLP method", weights);
  }

  FlagBkgHits2::~FlagBkgHits2(){}


  void FlagBkgHits2::beginJob()
  {
     _clusterer->init();
     if(_useMVA) _bkgMVA.initMVA();
  }


  void FlagBkgHits2::produce(art::Event& event )
  {
     unsigned iev=event.id().event();
     if(_debug > 0 && (iev%_printfreq)==0) std::cout<<"FlagBkgHits2: event="<<iev<<std::endl;


     art::Handle<StrawHitCollection> strawHitsHandle;
     event.getByLabel(_shtag,strawHitsHandle);
     const StrawHitCollection& shcol(*strawHitsHandle);

     art::Handle<StrawHitPositionCollection> strawHitPosHandle;
     event.getByLabel(_shptag,strawHitPosHandle);
     const StrawHitPositionCollection& shpcol(*strawHitPosHandle);

     art::Handle<StrawHitFlagCollection> strawHitFlagsHandle;
     event.getByLabel(_shftag,strawHitFlagsHandle);
     const StrawHitFlagCollection& shfcol(*strawHitFlagsHandle);

     //art::Handle<StereoHitCollection> stereoHitsHandle;
     //event.getByLabel(_sthTag,stereoHitsHandle);
     //const StereoHitCollection& sthcol(*stereoHitsHandle);

     std::unique_ptr<StrawHitFlagCollection> bkgfcol(new StrawHitFlagCollection(shfcol));
     std::unique_ptr<BkgClusterCollection> bkgccol(new BkgClusterCollection);
     std::unique_ptr<BkgQualCollection> bkgqcol(new BkgQualCollection);


     _clusterer->findClusters(*bkgccol,shcol,shpcol,shfcol);

     bkgqcol->reserve(bkgfcol->size());
     classifyClusters(*bkgccol, *bkgqcol, shcol, shpcol);


     for (auto const& cluster : *bkgccol)
     {  
       for (auto const& chit : cluster.hits()) 
         bkgfcol->at(chit.index()).merge(StrawHitFlag::bkgclust);

       if (cluster.flag().hasAllProperties(BkgClusterFlag::bkg))
       {
          for (const auto& chit : cluster.hits()) 
          {
            if (_flagall || chit.flag().hasAllProperties(StrawHitFlag::active))
	      bkgfcol->at(chit.index()).merge(StrawHitFlag::bkg);	
          }
       }

       if (cluster.hits().size() <= _maxisolated)
         for (auto const& chit : cluster.hits()) 
	     bkgfcol->at(chit.index()).merge(StrawHitFlag::isolated);      
     }


     event.put(std::move(bkgccol));
     event.put(std::move(bkgqcol));
     event.put(std::move(bkgfcol));
  }





  void FlagBkgHits2::classifyClusters(BkgClusterCollection& clusters,BkgQualCollection& cquals, 
                                      const StrawHitCollection& shcol, const StrawHitPositionCollection& shpcol) const 
  {
      const TTracker& tracker = dynamic_cast<const TTracker&>(getTrackerOrThrow());
      for (auto& cluster : clusters)
      { 
         BkgQual cqual;
         fillBkgQual(cluster,cqual,tracker, shcol,shpcol);

         if(cqual.MVAOutput() > _bkgMVAcut) 
           cluster._flag.merge(BkgClusterFlag::bkg);

         // ALWAYS record a quality to keep the vectors in sync.
         cquals.push_back(cqual);
     }
  }




  void FlagBkgHits2::fillBkgQual(const BkgCluster& cluster, BkgQual& cqual, const TTracker& tracker, 
                                 const StrawHitCollection& shcol, const StrawHitPositionCollection& shpcol) const 
  {

    unsigned nactive, nstereo;
    countHits(cluster,nactive,nstereo);
    cqual.setMVAStatus(BkgQual::unset);

    if (nactive >= _minnhits && nstereo >= _minnstereo)
    {
        cqual[BkgQual::nhits] = nactive;
        cqual[BkgQual::sfrac] = static_cast<double>(nstereo)/nactive;
        cqual[BkgQual::crho] = cluster.pos().perp();

        countPlanes(cluster,cqual,tracker,shcol);
        if (cqual[BkgQual::np] >= _minnp)
        {
           accumulator_set<double, stats<tag::weighted_variance>, double> racc;
           accumulator_set<double, stats<tag::variance(lazy)> > tacc;
           std::vector<double> hz;
           for (const auto& chit : cluster.hits())
           {
	       //std::cout<<chit.index()<<" ";
               if (chit.flag().hasAllProperties(StrawHitFlag::active))
               {
	           const StrawHit& sh = shcol.at(chit.index());
	           const StrawHitPosition& shp = shpcol.at(chit.index());
	           hz.push_back(shp.pos().z());
	           double dt = sh.time() - cluster.time();
	           tacc(dt);
	           // compute the transverse radius of this hit WRT the cluster center
	           CLHEP::Hep3Vector psep = (shp.pos()-cluster.pos()).perpPart();
	           double rho = psep.mag();
	           // project the hit position error along the radial direction to compute the weight
	           CLHEP::Hep3Vector pdir = psep.unit();
	           CLHEP::Hep3Vector tdir(-shp.wdir().y(),shp.wdir().x(),0.0);
	           double rwerr = shp.posRes(StrawHitPosition::wire)*pdir.dot(shp.wdir());
	           double rterr = shp.posRes(StrawHitPosition::trans)*pdir.dot(tdir);
                   // include the cluster center position error when computing the weight
	           double rwt = 1.0/sqrt(rwerr*rwerr + rterr*rterr + _cperr2/nactive);
	           racc(rho,weight=rwt);
	       }
           }
           cqual[BkgQual::hrho] = extract_result<tag::weighted_mean>(racc);
           cqual[BkgQual::shrho] = sqrt(std::max(extract_result<tag::weighted_variance>(racc),0.0));
           cqual[BkgQual::sdt] = sqrt(std::max(extract_result<tag::variance>(tacc),0.0));

           // find the min, max and gap from the sorted Z positions
           std::sort(hz.begin(),hz.end());
           cqual[BkgQual::zmin] = hz.front();
           cqual[BkgQual::zmax] = hz.back();

           // find biggest Z gap
           double zgap = 0.0;
           for(unsigned iz=1;iz<hz.size();++iz)
	     if(hz[iz]-hz[iz-1] > zgap)zgap = hz[iz]-hz[iz-1]; 
           cqual[BkgQual::zgap] = zgap;

           // compute MVA
           cqual.setMVAStatus(BkgQual::filled);
           if (_useMVA)
           {
             // reduce values down to what's actually used in the MVA.  This functionality
               // should be in MVATool, FIXME!
	       _mvavars[0] = cqual.varValue(BkgQual::hrho);
	       _mvavars[1] = cqual.varValue(BkgQual::shrho);
	       _mvavars[2] = cqual.varValue(BkgQual::crho);
	       _mvavars[3] = cqual.varValue(BkgQual::zmin);
	       _mvavars[4] = cqual.varValue(BkgQual::zmax);
	       _mvavars[5] = cqual.varValue(BkgQual::zgap);
	       _mvavars[6] = cqual.varValue(BkgQual::np);
	       _mvavars[7] = cqual.varValue(BkgQual::npfrac);
	       _mvavars[8] = cqual.varValue(BkgQual::nhits);
	       double readerout = _reader.EvaluateMVA( "MLP method" );
               // std::vector<double> mvavars(9,0.0);
               // for(size_t ivar =0;ivar<_mvavars.size(); ++ivar)
               // mvavars[ivar] = _mvavars[ivar];
               // double toolout = _bkgMVA.evalMVA(mvavars);
               // if(fabs(toolout - readerout) > 0.01) 
               // std::cout << "mva disagreement: MVATool " << toolout << " TMVA::Reader: " << readerout << std::endl;
	       cqual.setMVAValue(readerout);
	       cqual.setMVAStatus(BkgQual::calculated);

           }
        } else {
            cqual[BkgQual::hrho] = -1.0;
            cqual[BkgQual::shrho] = -1.0;
            cqual[BkgQual::sdt] = -1.0;
            cqual[BkgQual::zmin] = -1.0;
            cqual[BkgQual::zmax] = -1.0;
            cqual[BkgQual::zgap] = -1.0;
        }
     }
  }

  void FlagBkgHits2::countPlanes(const BkgCluster& cluster, BkgQual& cqual, const TTracker& tracker, 
                                 const StrawHitCollection& shcol) const 
  {
      std::vector<int> hitplanes(tracker.nPlanes(),0);
      for (const auto& chit : cluster.hits()) 
      {
         if (chit.flag().hasAllProperties(StrawHitFlag::active))
         {
            const StrawHit& sh = shcol.at(chit.index());
	    unsigned iplane = (unsigned)(tracker.getStraw(sh.strawIndex()).id().getPlaneId());
	    ++hitplanes[iplane];
         }
      }

      unsigned ipmin = 0;
      unsigned ipmax = tracker.nPlanes()-1;
      while (hitplanes[ipmin]==0) ++ipmin;
      while (hitplanes[ipmax]==0) --ipmax;

      unsigned npexp(0);
      unsigned np(0);
      double nphits(0.0);
      const std::vector<Plane>& planes = tracker.getPlanes();

      for(unsigned ip = ipmin; ip <= ipmax; ++ip)
      {
         if (planes[ip].exists()) ++npexp;
         if (hitplanes[ip]> 0)++np;
         nphits += hitplanes[ip];
      }
      cqual[BkgQual::np] = np; 
      cqual[BkgQual::npexp] = npexp; 
      cqual[BkgQual::npfrac] = np/static_cast<float>(npexp); 
      cqual[BkgQual::nphits] = nphits/static_cast<float>(np);
  }

  void FlagBkgHits2::countHits(const BkgCluster& cluster, unsigned& nactive, unsigned& nstereo) const 
  {
      nactive = nstereo = 0;
      for (const auto& chit : cluster.hits()) 
      {
         if (chit.flag().hasAllProperties(StrawHitFlag::active))
         {
	    ++nactive;
	    if (chit.flag().hasAllProperties(StrawHitFlag::stereo))++nstereo;
         }
     }
  }

}

using mu2e::FlagBkgHits2;
DEFINE_ART_MODULE(FlagBkgHits2);
