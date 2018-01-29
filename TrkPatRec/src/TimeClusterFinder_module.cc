//
// TTracker time peak finder
//
// $Id: TimeClusterFinder_module.cc,v 1.3 2014/08/25 12:08:29 tassiell Exp $
// $Author: tassiell $
// $Date: 2014/08/25 12:08:29 $
//
// Original author D. Brown and G. Tassielli
//

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
// Mu2e
#include "GeneralUtilities/inc/Angles.hh"
#include "Mu2eUtilities/inc/MVATools.hh"
// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
// tracking
#include "TrkReco/inc/TrkUtilities.hh"
#include "TrkReco/inc/TrkTimeCalculator.hh"
// root
#include "TH1F.h"
// boost
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/weighted_median.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/weighted_mean.hpp>
#include <boost/accumulators/statistics/weighted_variance.hpp>
// C++
#include <memory>
#include <algorithm>
#include <utility>
using namespace std; 
using namespace boost::accumulators;
using CLHEP::Hep3Vector;



namespace {

   struct TimePeakMVA
   {
      vector<Double_t> _pars;
      Double_t& _dt;
      Double_t& _dphi;
      Double_t& _rho;
      TimePeakMVA() : _pars(3,0.0), _dt(_pars[0]), _dphi(_pars[1]), _rho(_pars[2]) {}
   };
   
   
   struct meanAccumulator 
   {
      meanAccumulator(): sum(0),weight(0) {};
      void add(double m, double w) {sum +=m*w; weight +=w;}
      void remove(double m, double w) {sum -=m*w; weight -=w;}
      double mean() {return weight>0 ? sum/weight : 0;}
      double sum;
      double weight;
   };
 
   
}



namespace mu2e {


  class TimeClusterFinder : public art::EDProducer {
    public:
      enum ClusterAlgorithm {peak=0,scan};
      explicit TimeClusterFinder(fhicl::ParameterSet const& pset);
      virtual ~TimeClusterFinder();

      virtual void beginJob();
      void produce(art::Event & e);

    private:
        int               _iev;
        int               _debug;
        int               _printfreq;
        ClusterAlgorithm  _algo;
        art::InputTag     _shTag;
        art::InputTag     _shpTag;
        art::InputTag     _shfTag;
        art::InputTag     _ccFastTag;
        StrawHitFlag      _hsel, _hbkg;
        double            _maxdt;
        unsigned          _minnhits;
        double            _minpeakmva, _maxpeakdt;
        double            _maxdPhi;
        double            _tmin, _tmax, _tbin;
        TH1F	          _timespec;
        double            _ymin;
        bool              _refine;
        bool              _preFilter;
        bool              _usecc;
        double            _ccmine, _ccwt;
        double	          _maxover;
        MVATools          _peakMVA; // MVA for peak cleaning
        TimePeakMVA       _pmva; // input variables to TMVA for peak cleaning
        TrkTimeCalculator _ttcalc;
        int               _deltaNbins; 

        typedef std::pair<Float_t,int> BinContent;

        void findClusters(TimeClusterCollection& tclusts, StrawHitFlagCollection& flags, 
                          const StrawHitCollection& shcol, const StrawHitPositionCollection& shpcol, 
                          const StrawHitFlagCollection& shfcol, 
                          const art::Handle<CaloClusterCollection> CaloClusterHandle, 
                          const CaloClusterCollection* caloClusters);
        void fillTimeSpectrum(const StrawHitCollection& shcol, const StrawHitPositionCollection& shpcol,
                              const std::vector<int>& goodHitFlag,const CaloClusterCollection* caloClusters);
        void initCluster(TimeCluster& tclust, const StrawHitCollection& shcol, const StrawHitPositionCollection& shpcol,
                         const std::vector<int>& goodHitFlag);
        void refineCluster(TimeCluster& tclust, const StrawHitCollection& shcol,const StrawHitPositionCollection& shpcol);
        void findPeaks(std::vector<double>& tctimes);
        void scanPeaks(std::vector<double>& tctimes);
        bool goodHit(const StrawHitFlag& flag) const; 
        double clusterWeight(const TimeCluster& tc) const;
        void addCaloCluster(TimeCluster& tc, const art::Handle<CaloClusterCollection> CaloClusterHandle, 
                            const CaloClusterCollection* caloClusters); 



    };

    TimeClusterFinder::~TimeClusterFinder() {
    }

    TimeClusterFinder::TimeClusterFinder(fhicl::ParameterSet const& pset) :
      _debug             (pset.get<int>("debugLevel",0)),
      _printfreq         (pset.get<int>("printFrequency",101)),
      _algo	         (static_cast<ClusterAlgorithm>(pset.get<int>("ClusterAlgorithm",peak))),
      _shTag	         (pset.get<art::InputTag>("StrawHitCollection","makeSH")),
      _shpTag	         (pset.get<art::InputTag>("StrawHitPositionCollection","MakeStereoHits")),
      _shfTag	         (pset.get<art::InputTag>("StrawHitFlagCollection","FlagBkgHits")),
      _ccFastTag         (pset.get<art::InputTag>("caloClusterModuleLabel","CaloClusterFast")),
      _hsel	         (pset.get<std::vector<std::string> >("HitSelectionBits",vector<string>{"EnergySelection","TimeSelection","RadiusSelection"})),
      _hbkg              (pset.get<vector<string> >("HitBackgroundBits",vector<string>{"Background"})),
      _maxdt             (pset.get<double>(  "DtMax",30.0)),
      _minnhits          (pset.get<unsigned>("MinNHits",10)),
      _minpeakmva        (pset.get<double>(  "MinTimePeakMVA",0.2)),
      _maxpeakdt         (pset.get<double>(  "MaxTimePeakDeltat",25.0)),
      _maxdPhi           (pset.get<double>(  "MaxdPhi",1.2)),
      _tmin              (pset.get<double>(  "tmin",500.0)),
      _tmax              (pset.get<double>(  "tmax",1700.0)),
      _tbin              (pset.get<double>(  "tbin",15.0)),
      _ymin              (pset.get<double>(  "ymin",5.0)),
      _refine	         (pset.get<bool>(    "RefineClusters",true)),    
      _preFilter         (pset.get<bool>(    "PrefilterCluster",true)),    
      _usecc	         (pset.get<bool>(    "UseCaloCluster",false)),
      _ccmine            (pset.get<double>(  "CaloClusterMinE",50.0)),   
      _ccwt              (pset.get<double>(  "CaloClusterWeight",10.0)), 
      _maxover           (pset.get<double>(  "MaxOverlap",0.3)),         // Maximum hit overlap to consider clusters as different
      _peakMVA           (pset.get<fhicl::ParameterSet>("PeakCleanMVA",fhicl::ParameterSet())),
      _ttcalc            (pset.get<fhicl::ParameterSet>("T0Calculator",fhicl::ParameterSet()))
    {    
        unsigned nbins = (unsigned)rint((_tmax-_tmin)/_tbin);
        _timespec = TH1F("timespec","time spectrum",nbins,_tmin,_tmax);
        _deltaNbins = int(_maxdt/_tbin)-1;

        produces<TimeClusterCollection>();
        produces<StrawHitFlagCollection>();
    }

    void TimeClusterFinder::beginJob()
    {
       _peakMVA.initMVA();
       if (_debug > 0)
       {
         std::cout << "TimeClusterFinder MVA : " << std::endl; 
         _peakMVA.showMVA();
       }
    }


    //--------------------------------------------------------------------------------------------------------------
    void TimeClusterFinder::produce(art::Event & event )
    {     
       if (_debug > 0 && (event.id().event()%_printfreq)==0) std::cout<<"TimeClusterFinder: event="<<event.id().event()<<std::endl;

       const CaloClusterCollection* caloClusters(0);
       art::Handle<CaloClusterCollection> CaloClusterHandle;
       if (event.getByLabel(_ccFastTag, CaloClusterHandle)) caloClusters = CaloClusterHandle.product();

       art::Handle<StrawHitCollection> strawHitsHandle;
       event.getByLabel(_shTag,strawHitsHandle);
       const StrawHitCollection& shcol(*strawHitsHandle);

       art::Handle<StrawHitPositionCollection> strawHitPosHandle;
       event.getByLabel(_shpTag,strawHitPosHandle);
       const StrawHitPositionCollection& shpcol(*strawHitPosHandle);

       art::Handle<StrawHitFlagCollection> strawHitFlagsHandle;
       event.getByLabel(_shfTag,strawHitFlagsHandle);
       const StrawHitFlagCollection& shfcol(*strawHitFlagsHandle);

       if (_usecc && caloClusters==nullptr)
          throw cet::exception("RECO")<<"mu2e::TimeClusterFinder: No caloCluster collection but useCalorimeter flag set to true" << std::endl; 

       std::unique_ptr<TimeClusterCollection> timeClusterColl(new TimeClusterCollection);
       std::unique_ptr<StrawHitFlagCollection> flagColl(new StrawHitFlagCollection(shfcol));

       
       findClusters(*timeClusterColl,*flagColl,shcol,shpcol, shfcol, CaloClusterHandle, caloClusters);
       

       if (_debug > 0) std::cout << "Found " << timeClusterColl->size() << " Time Clusters " << std::endl;
       for (auto tpc : *timeClusterColl)
       {
          if (_debug > 1)
	    std::cout << "Time Cluster time = " << tpc.t0().t0() << " +- " << tpc.t0().t0Err()
	    << " position = " << tpc._pos << std::endl;

          for (auto shi : tpc._strawHitIdxs )
          {
	     flagColl->at(shi).merge(StrawHitFlag::tclust);
	     if (_debug > 3) std::cout << "Time Cluster hit at index " << shi << std::endl;
          }
       }

       event.put(std::move(timeClusterColl));
       event.put(std::move(flagColl));
    }





    //--------------------------------------------------------------------------------------------------------------
    void TimeClusterFinder::findClusters(TimeClusterCollection& tclusts, StrawHitFlagCollection& flags, 
                                          const StrawHitCollection& shcol, const StrawHitPositionCollection& shpcol, 
                                          const StrawHitFlagCollection& shfcol, 
                                          const art::Handle<CaloClusterCollection> CaloClusterHandle, 
                                          const CaloClusterCollection* caloClusters) 
    {

      //Strangely, this saves a lot of time!!!
      std::vector<int> goodHitFlag(shcol.size(),0);       
      for (unsigned istr=0; istr<shcol.size();++istr) 
         if (shfcol.at(istr).hasAllProperties(_hsel) && !shfcol.at(istr).hasAnyProperty(_hbkg)) goodHitFlag[istr]=1;

      fillTimeSpectrum(shcol, shpcol, goodHitFlag,  caloClusters);


      vector<double> tctimes;
      switch (_algo )
      {
        case peak : default:
	  findPeaks(tctimes);
	  break;
        case scan :
	  scanPeaks(tctimes);
	  break;
      }


      for (auto tctime: tctimes)
      {
        if(_debug > 1) std::cout << "Peak Time " << tctime  << std::endl;

        TimeCluster tclust;
        tclust._t0 = TrkT0(tctime,1.0);


        for(size_t istr=0; istr<shcol.size(); ++istr)
        {
	   if (goodHitFlag[istr])
           {
	      double time = _ttcalc.strawHitTime(shcol.at(istr),shpcol.at(istr));
	      if (fabs(time-tctime) < _maxdt)tclust._strawHitIdxs.push_back(StrawHitIndex(istr));	  
	   }
        }

        if (_usecc) addCaloCluster(tclust,CaloClusterHandle, caloClusters);

        initCluster(tclust, shcol, shpcol, goodHitFlag);
        if (_refine) refineCluster(tclust,shcol,shpcol);


        if (tclust._strawHitIdxs.size() >= _minnhits)
        {
	    bool overl(false);
	    for (auto itclust = tclusts.begin(); itclust < tclusts.end(); ++itclust)
            {
	      if (TrkUtilities::overlap(tclust,*itclust) < _maxover) continue;          
	      if (tclust._strawHitIdxs.size() > itclust->_strawHitIdxs.size())
              { 
	         tclusts.erase(itclust);
	         break;
	      } else {	    
	         overl = true;
	         break;
	      }	  
	   }
           if (!overl) tclusts.push_back(tclust);
        }
        //std::cout<<"Collection size final"<<tclust._strawHitIdxs.size()<<std::endl;
      }

      for (auto tpc : tclusts)
         for (auto shi : tpc._strawHitIdxs ) flags.at(shi).merge(StrawHitFlag::tclust);


       // debug test of histogram
       if (_debug > 2)
       {
          art::ServiceHandle<art::TFileService> tfs;
          TH1F* tspec = tfs->make<TH1F>(_timespec);
          char name[40];
          char title[100];
          snprintf(name,40,"tspec_%i",_iev);
          snprintf(title,100,"time spectrum event %i;nsec",_iev);
          tspec->SetNameTitle(name,title);
       }
    }



    //--------------------------------------------------------------------------------------------------------------
    void TimeClusterFinder::fillTimeSpectrum(const StrawHitCollection& shcol, const StrawHitPositionCollection& shpcol,
                                              const std::vector<int>& goodHitFlag, const CaloClusterCollection* caloClusters)
    {
        _timespec.Reset();

        for (unsigned istr=0; istr<shcol.size();++istr)
        {
           if (!goodHitFlag[istr]) continue;
           double time = _ttcalc.strawHitTime(shcol.at(istr),shpcol.at(istr));
           _timespec.Fill(time);
        }

       // weight the cluster WRT hits by an ad-hoc value.  This is more about signal/noise than resolution
       if (_usecc)
       {
           for(auto& calo : *caloClusters)
           {
	      if (calo.energyDep() < _ccmine) return;
              double time = _ttcalc.caloClusterTime(calo);
              _timespec.Fill(time, _ccwt);	
           }
       }
    }

    //--------------------------------------------------------------------------------------------------------------
    void TimeClusterFinder::findPeaks(std::vector<double>& tctimes)
    {
       tctimes.clear();

       std::vector<BinContent> bcv;
       int nbins = _timespec.GetNbinsX()+1;
       std::vector<bool> alreadyUsed(nbins,false);

       for (int ibin=1;ibin < nbins; ++ibin)
         if (_timespec.GetBinContent(ibin) >= _ymin) bcv.push_back(make_pair(_timespec.GetBinContent(ibin),ibin));
       std::sort(bcv.begin(),bcv.end(),[](const BinContent& x, const BinContent& y){return x.first > y.first;});

       for (const auto& bc : bcv)
       {
           if (alreadyUsed[bc.second]) continue;

	   double tctime(0.0);
	   double norm(0.0);
	   for (int ibin = std::max(1,bc.second-_deltaNbins);ibin < std::min(nbins,bc.second+_deltaNbins+1); ++ibin)
           {
	       norm += _timespec.GetBinContent(ibin);
	       tctime += _timespec.GetBinCenter(ibin)*_timespec.GetBinContent(ibin);
               alreadyUsed[ibin] = true;
	   }
	   tctime /= norm;  
	   if (norm > _minnhits) tctimes.push_back(tctime);
       }
    }

    //--------------------------------------------------------------------------------------------------------------
    void TimeClusterFinder::scanPeaks(std::vector<double>& tctimes)
    {
       //implement this FIXME!!
    }


    //--------------------------------------------------------------------------------------------------------------
    void TimeClusterFinder::initCluster(TimeCluster& tclust, const StrawHitCollection& shcol, 
                                         const StrawHitPositionCollection& shpcol, const std::vector<int>& goodHitFlag)
    {
       // use medians to initialize robustly
       accumulator_set<double, stats<tag::min > > tmin;
       accumulator_set<double, stats<tag::max > > tmax;
       accumulator_set<double, stats<tag::weighted_median(with_p_square_quantile) >, double > tacc;
       accumulator_set<double, stats<tag::median(with_p_square_quantile) > > xacc;
       accumulator_set<double, stats<tag::median(with_p_square_quantile) > > yacc;
       accumulator_set<double, stats<tag::median(with_p_square_quantile) > > zacc;

       unsigned nstrs = tclust._strawHitIdxs.size();

       for(auto ish :tclust._strawHitIdxs)
       { 
           if (!goodHitFlag[ish]) continue;
	   const CLHEP::Hep3Vector& pos = shpcol.at(ish).pos();
	   double htime = _ttcalc.strawHitTime(shcol.at(ish),shpcol.at(ish));
	   double wt = std::pow(1.0/_ttcalc.strawHitTimeErr(),2);
	   tmin(htime);
	   tmax(htime);
	   tacc(htime,weight=wt);
	   xacc(pos.x());
	   yacc(pos.y());
	   zacc(pos.z());      
       }

       if (tclust._caloCluster.isNonnull())
       {
          double ctime = _ttcalc.caloClusterTime(*tclust._caloCluster);
          double wt = std::pow(1.0/_ttcalc.caloClusterTimeErr(tclust._caloCluster->diskId()),2);
          tacc(ctime,weight=wt);
       }

       static double invsqrt12(1.0/sqrt(12.0));
       tclust._t0._t0 = extract_result<tag::weighted_median>(tacc);
       tclust._t0._t0err = ( boost::accumulators::extract::max(tmax)-boost::accumulators::extract::min(tmin))*invsqrt12/sqrt(nstrs);
       tclust._pos = CLHEP::Hep3Vector(median(xacc),median(yacc),median(zacc));

       if (_debug > 0) std::cout<<"Init time peak "<<tclust._t0._t0<<std::endl;    
    }


    //--------------------------------------------------------------------------------------------------------------
    void TimeClusterFinder::refineCluster(TimeCluster& tclust, const StrawHitCollection& shcol, 
                                           const StrawHitPositionCollection& shpcol)
    {
       double pphi(tclust._pos.phi());
       double ptime(tclust._t0.t0());              
       meanAccumulator faccu,taccu;
       

       
       if (_preFilter)
       {
           /*
           //preFilter all in one shot
           std::vector<int> toremove;    
           for (size_t ips=0; ips<tclust._strawHitIdxs.size(); ++ips)
           {
	      unsigned ish = tclust._strawHitIdxs[ips];
	      double phi   = shpcol.at(ish).pos().phi(); 
	      double dphi  = fabs(Angles::deltaPhi(phi,pphi));                
              if (dphi > _maxdPhi) toremove.push_back(ips);
           }

           for (auto irm=toremove.rbegin();irm!=toremove.rend();++irm)
           {
              std::swap(tclust._strawHitIdxs[*irm],tclust._strawHitIdxs.back());
              tclust._strawHitIdxs.pop_back();
           }

           accumulator_set<double, stats<tag::mean > > facc;
           for (size_t ips=0;ips<tclust._strawHitIdxs.size();++ips)
           {
	      unsigned ish = tclust._strawHitIdxs[ips];
	      double   phi = shpcol.at(ish).pos().phi();
	      Angles::deltaPhi(phi,pphi);
	      facc(phi);
           }
           pphi  = extract_result<tag::mean>(facc); 
           */

           //prefilter one after another         
           while (tclust._strawHitIdxs.size() >= _minnhits)
           {                    
               int iworst(-1);
               double maxadPhi(0);
               double maxdPhi(0),sumphi(0);

               for (size_t ips=0; ips<tclust._strawHitIdxs.size(); ++ips)
               {
	          unsigned ish = tclust._strawHitIdxs[ips];
	          double phi   = shpcol.at(ish).phi(); 
	          double dphi  = Angles::deltaPhi(phi,pphi);
	          double adphi = std::abs(dphi);
                  sumphi += phi;

                  if (adphi > maxadPhi) {iworst=ips; maxadPhi=adphi; maxdPhi=phi;}             
               }

               if (maxadPhi<_maxdPhi) break;

               std::swap(tclust._strawHitIdxs[iworst],tclust._strawHitIdxs.back());
	       tclust._strawHitIdxs.pop_back();      
               pphi  = (sumphi-maxdPhi)/float(tclust._strawHitIdxs.size());
          }
       }



       while (tclust._strawHitIdxs.size() >= _minnhits)
       {
          size_t iworst(0);
          double worstmva(100.0);

          for (size_t ips=0;ips<tclust._strawHitIdxs.size();++ips)
          {
	     unsigned ish = tclust._strawHitIdxs[ips];
	     double dt = _ttcalc.strawHitTime(shcol.at(ish),shpcol.at(ish)) - ptime;

             double rho = shpcol.at(ish).pos().perp();
	     double phi = shpcol.at(ish).phi(); 
	     double dphi = Angles::deltaPhi(phi,pphi);

	     _pmva._dt = dt;
	     _pmva._dphi = dphi;
	     _pmva._rho = rho;
             
             double mvaout = _peakMVA.evalMVA(_pmva._pars);
             if (mvaout < worstmva)
             {
	        worstmva = mvaout;
	        iworst = ips;
	     }
          }
          
          if (worstmva > _minpeakmva) break;

          std::swap(tclust._strawHitIdxs[iworst],tclust._strawHitIdxs.back());
          tclust._strawHitIdxs.pop_back();      
	  
          // re-compute the average phi and range
          accumulator_set<double, stats<tag::mean > > facc;
          accumulator_set<double, stats<tag::weighted_mean >, double > tacc;
          for (size_t ips=0;ips<tclust._strawHitIdxs.size();++ips)
          {
            unsigned ish = tclust._strawHitIdxs[ips];
            double  time = _ttcalc.strawHitTime(shcol.at(ish),shpcol.at(ish));
            double    wt = std::pow(1.0/_ttcalc.strawHitTimeErr(),2);
            double   phi = shpcol.at(ish).phi();
            Angles::deltaPhi(phi,pphi);
            tacc(time,weight=wt);
            facc(phi);
          }
          pphi  = extract_result<tag::mean>(facc);
          ptime = extract_result<tag::weighted_mean>(tacc);
       }

       
       // final pass: hard cut on dt 
       std::vector<size_t> toremove;
       accumulator_set<double, stats<tag::mean > > facc;
       accumulator_set<double, stats<tag::weighted_variance(lazy)>, double > terr;
       accumulator_set<double, stats<tag::mean > > racc;
       accumulator_set<double, stats<tag::mean > > zacc;
       for(size_t ips=0;ips<tclust._strawHitIdxs.size();++ips)
       {
          unsigned ish = tclust._strawHitIdxs[ips];
          double dt = _ttcalc.strawHitTime(shcol.at(ish),shpcol.at(ish)) - ptime;
          double wt = std::pow(1.0/_ttcalc.strawHitTimeErr(),2);
          double phi = shpcol.at(ish).phi();
          double rho = shpcol.at(ish).pos().perp();
          Angles::deltaPhi(phi,pphi);
          if (fabs(dt) < _maxpeakdt)
          {
	     terr(_ttcalc.strawHitTime(shcol.at(ish),shpcol.at(ish)),weight=wt);
             facc(phi);
	     racc(rho);
	     zacc(shpcol.at(ish).pos().z());
          } else {
	     toremove.push_back(ips);
          }
       }

       for (auto irm=toremove.rbegin();irm!=toremove.rend();++irm)
       {
          std::swap(tclust._strawHitIdxs[*irm],tclust._strawHitIdxs.back());
          tclust._strawHitIdxs.pop_back();
       }

       tclust._t0._t0 = extract_result<tag::weighted_mean>(terr);
       tclust._t0._t0err = sqrt(std::max(0.0,extract_result<tag::weighted_variance(lazy)>(terr))/extract_result<tag::count>(terr));
       pphi = extract_result<tag::mean>(facc);
       double prho = extract_result<tag::mean>(racc);
       double zpos = extract_result<tag::mean>(zacc);
       tclust._pos = CLHEP::Hep3Vector(prho*cos(pphi),prho*sin(pphi),zpos);

       //if (_debug > 0) std::cout<<"final time "<<tclust._t0._t0<<std::endl;    
    }

    //--------------------------------------------------------------------------------------------------------------
    void TimeClusterFinder::addCaloCluster(TimeCluster& tc, const art::Handle<CaloClusterCollection> CaloClusterHandle, 
                                            const CaloClusterCollection* caloClusters) 
    {
        auto bestcc = caloClusters->end();
        for (auto icc = caloClusters->begin();icc != caloClusters->end(); ++icc)
        {
	    if (icc->energyDep() > _ccmine) continue;
	    double time = _ttcalc.caloClusterTime(*icc);
	    if (fabs(tc._t0._t0-time) > _maxdt) continue;
	    if (bestcc == caloClusters->end() || icc->energyDep() > bestcc->energyDep()) bestcc = icc;
        }

        if (bestcc != caloClusters->end())
        {
	   size_t index = std::distance(caloClusters->begin(),bestcc);
	   tc._caloCluster = art::Ptr<CaloCluster>(CaloClusterHandle,index);
        }
    }


    //--------------------------------------------------------------------------------------------------------------
    bool TimeClusterFinder::goodHit(const StrawHitFlag& flag) const 
    {
       return flag.hasAllProperties(_hsel) && !flag.hasAnyProperty(_hbkg);
    }

    //--------------------------------------------------------------------------------------------------------------
    double TimeClusterFinder::clusterWeight(const TimeCluster& tc) const
    {
       double wt = tc._strawHitIdxs.size();
       if (tc._caloCluster.isNonnull()) wt += _ccwt;
       return wt;
    }
  
}  

using mu2e::TimeClusterFinder;
DEFINE_ART_MODULE(TimeClusterFinder);

