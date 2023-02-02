///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art_root_io/TFileService.h"
#include "art/Utilities/make_tool.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "Offline/RecoDataProducts/inc/StrawHitPosition.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"

#include "Offline/TrkReco/inc/TrkUtilities.hh"
#include "Offline/TrkReco/inc/TrkTimeCalculator.hh"
#include "Offline/GeneralUtilities/inc/Angles.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "Offline/Mu2eUtilities/inc/polyAtan2.hh"
#include "Offline/CalPatRec/inc/PhiClusterFinder_types.hh"

#include "TH1.h"

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/weighted_median.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/weighted_mean.hpp>

#include <iostream>
#include <fstream>
#include <memory>
#include <algorithm>
#include <utility>
#include <functional>
#include <vector>

using namespace std;
using namespace boost::accumulators;

namespace mu2e {

  using namespace PhiClusterFinderTypes;

  class PhiClusterFinder : public art::EDProducer {
    public:

      struct Config
      {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;
        fhicl::Atom<int>                                    diag{                    Name("diag"),                 Comment("Diag"),0 };
        fhicl::Atom<int>                                    debug{                   Name("debug"),                Comment("Debug"),0 };
        fhicl::Atom<int>                                    minNSH{                  Name("minNSH"),               Comment("Minimum number of straw hits") };
        fhicl::Atom<int>                                    threshold{               Name("threshold"),            Comment("Lower limit of hits in a phi cluster") };
        fhicl::Atom<int>                                    nPhiClusters{            Name("nPhiClusters"),         Comment("Possible number of Phi Clusters in an event") };
        fhicl::Atom<float>                                  phiMin{                  Name("phiMin"),               Comment("Phi histogram start") };
        fhicl::Atom<float>                                  phiMax{                  Name("phiMax"),               Comment("Phi histogram end") };
        fhicl::Atom<int>                                    nPhiBins{                Name("nPhiBins"),             Comment("Phi histogram N(bins)") };
        fhicl::Atom<float>                                  minDPhi{                 Name("minDPhi"),              Comment("Minimum delta phi between two phi clusters") };
        fhicl::Atom<float>                                  pitch{                   Name("pitch"),                Comment("Average helix pitch (= dz/dflight, =sin(lambda)") };
        fhicl::Atom<float>                                  tMin{                    Name("tMin"),                 Comment("Time histogram start") };
        fhicl::Atom<float>                                  tMax{                    Name("tMax"),                 Comment("Time histogram end") };
        fhicl::Atom<float>                                  tBin{                    Name("tBin"),                 Comment("Time histogram bin width") };
        fhicl::Table<TrkTimeCalculator::Config>             t0Calc{                  Name("t0Calc"),               Comment("TimeTracker calculator config") };
        fhicl::Atom<float>                                  yMin{                    Name("yMin"),                 Comment("Minimum hit in time hist bin for peak") };
        fhicl::Atom<float>                                  maxDt{                   Name("maxDt"),                Comment("Maximum time difference between T0 and hit time") };
        fhicl::Atom<float>                                  minSigma{                Name("minSigma"),             Comment("Min sigma cut to remove compton particles") };
        fhicl::Atom<int>                                    peakWidth{               Name("peakWidth"),            Comment("Time Peak Width") };
        fhicl::Atom<art::InputTag>                          ComboHitCollection{      Name("ComboHitCollection"),   Comment("ComboHit collection name") };
        fhicl::Atom<art::InputTag>                          TimeClusterCollection{   Name("TimeClusterCollection"),Comment("TimeCluster collection name") };
        fhicl::Atom<bool>                                   useCC{                   Name("useCC"),                Comment("Use calo cluster") };
        fhicl::Table<PhiClusterFinderTypes::Config>         DiagPlugin{              Name("DiagPlugin"),           Comment("Diag plugin") };
      };

      explicit PhiClusterFinder(const art::EDProducer::Table<Config>& config);
      void beginJob();
      void produce(art::Event& event);

  private:
    typedef std::pair<Float_t,int> BinContent;
    int                                 _iev;
    int                                 _diag,_debug;
    int                                 _minnsh;
    int                                 _threshold;
    int                                 _nphiclusters;
    float                               _phimin,_phimax;
    int                                 _nphibins;
    float                               _mindphi,_pitch;
    float                               _tmin, _tmax, _tbin;
    TrkTimeCalculator                   _t0calc;
    float                               _ymin,_maxdt,_minsigma;
    int                                 _npeak;
    art::ProductToken<ComboHitCollection> const _chToken;
    art::ProductToken<TimeClusterCollection> const _tcToken;
    bool                                _useCC;
    const ComboHitCollection*           _chcol;
    TH1F*                               _hist1;
    TH1F                                _timespec;
    std::vector<int>                    _clustno;
    std::vector<float>                  _phitotal;
    std::vector<int>                    _nhitotal;
    //-----------------------------------------------------------------------------
    // diagnostics
    //-----------------------------------------------------------------------------
    Data_t                              _data;
    std::unique_ptr<ModuleHistToolBase> _hmanager;

    void findClusters (TimeClusterCollection& tccol1, const std::vector<StrawHitIndex>& ordchcol);
    void clusterminmax(float& cluphimin, float& cluphimax);
    void fillTimeSpectrum(int ClustNo, const std::vector<StrawHitIndex>& ordchcol);
    void fillTimeSpectrumAlt(const std::vector<StrawHitIndex>& ordchcol);
    void findPeaks  (TimeCluster& tc, const std::vector<StrawHitIndex>& ordchcol);
    void assignHits(TimeCluster& otc, const std::vector<StrawHitIndex>& ordchcol);
    void initCluster(TimeCluster& tc);
    float checkdelta(TimeCluster& tc, int ClustNo, const std::vector<StrawHitIndex>& ordchcol);
    void addCaloClusters(TimeClusterCollection& tccolnew, const TimeClusterCollection& tccol);
};

 PhiClusterFinder::PhiClusterFinder(const art::EDProducer::Table<Config>& config):
    art::EDProducer{config},
    _diag         (config().diag()),
    _debug        (config().debug()),
    _minnsh       (config().minNSH()),
    _threshold    (config().threshold()),
    _nphiclusters (config().nPhiClusters()),
    _phimin       (config().phiMin()),
    _phimax       (config().phiMax()),
    _nphibins     (config().nPhiBins()),
    _mindphi      (config().minDPhi()),
    _pitch        (config().pitch()),
    _tmin         (config().tMin()),
    _tmax         (config().tMax()),
    _tbin         (config().tBin()),
    _t0calc       (config().t0Calc()),
    _ymin         (config().yMin()),
    _maxdt        (config().maxDt()),
    _minsigma     (config().minSigma()),
    _npeak        (config().peakWidth()),
    _chToken      {consumes<ComboHitCollection>(config().ComboHitCollection()) },
    _tcToken      {consumes<TimeClusterCollection>(config().TimeClusterCollection()) },
    _useCC        (config().useCC())
    {
      art::ServiceHandle<art::TFileService> tfs;
      _hist1 = tfs->make<TH1F>( "hist1" , "phi spectrum",_nphibins,_phimin,_phimax);
      unsigned nbinst = (unsigned)rint((_tmax-_tmin)/_tbin);
      _timespec = TH1F("timespec","time spectrum",nbinst,_tmin,_tmax);
      produces<TimeClusterCollection>();

      if (_diag != 0) _hmanager = art::make_tool<ModuleHistToolBase>(config().DiagPlugin," ");
      else            _hmanager = std::make_unique<ModuleHistToolBase>();
    }

  void PhiClusterFinder::beginJob() {
    if (_diag > 0) {
      art::ServiceHandle<art::TFileService> tfs;
      _hmanager->bookHistograms(tfs);
    }
  }

  void PhiClusterFinder::produce(art::Event & event ){
    _iev = event.id().event();
    std::unique_ptr<TimeClusterCollection> tccol1(new TimeClusterCollection);
    auto const& chH = event.getValidHandle(_chToken);
    _chcol = chH.product();

    // Save the sum phi of a cluster to calculate delta phi between the clusters
    _phitotal.clear();
    _nhitotal.clear();
    _phitotal.resize(_nphiclusters);
    _nhitotal.resize(_nphiclusters);

    for(int m=0;m<_nphiclusters;m++){
      _phitotal[m] = 0.;
      _nhitotal[m] = 0.;
    }
    auto const& tcH = event.getValidHandle(_tcToken);
    const TimeClusterCollection& tccol(*tcH);
    // For the diagnostics
    _data._tccol = tcH.product();
    if (tccol.size() > 0) {
      for(size_t ipeak=0; ipeak<tccol.size(); ipeak++) {
        const auto& tc = tccol[ipeak];
        int nhits = tc.nhits();
        _clustno.resize(nhits);
        for(int k=0; k<nhits; k++) {
           _clustno[k] = -1;
        }
        // Find the time clusters separated in phi
        findClusters(*tccol1,tc.hits());
      }
    }
    if(_debug > 0 and tccol.size() > 0){
      int newsize = tccol1->size();
      int oldsize = tccol.size();
      int eventno = _iev;
      printf("Output tccol size = %i Input tccol size =  %i event = %i\n",newsize,oldsize,eventno);
    }

    TimeClusterCollection& tccolnew(*tccol1);
    // if(newsize == 2 and _debug>2) std::cout<<"Event with two phi clusters = "<<_iev<<std::endl;

    // Check for calo clusters
    if(_useCC and tccol.size() > 0){
      addCaloClusters(tccolnew,tccol);
    }
    // Save the output time cluster collection for the diagnostics
    _data._tccolnew = tccol1.get();
    // if(_debug>0) std::cout<<"Diag data size = "<<_data._tccolnew->size()<<" old = "<<_data._tccol->size()<<std::endl;
    if (_diag > 0) {
      _hmanager->fillHistograms(&_data);
    }
    // if (tccol1->size() == 0 and tccol.size()>0)
    // std::cout<<"tccol new size = "<<tccol1->size()<<" tccol size = "<<tccol.size()<<" event = "<<_iev<<std::endl;
    event.put(std::move(tccol1));
  }

//-----------------------------------------------------------------------------
// process input timecluster to find the cluster of hits separated in phi
//-----------------------------------------------------------------------------
  void PhiClusterFinder::findClusters(TimeClusterCollection& tccol1, const std::vector<StrawHitIndex>& ordchcol) {
    int nh = ordchcol.size();
    // Counts the no. of hits used to form the clusters
    int countusedhit(0);
    // Counts the no. of clusters
    int counter(0);
    // Count the no. of hits in a cluster when it is an event with two phi clusters
    int count1(0),count2(0);
    float dphi(0);
    // Mark true if a hit is used to form a cluster
    std::vector<bool> usedhit(ordchcol.size(),false);

    if (nh > (_minnsh/2)) {
      // Check if all the hits are used and keep trying to find peaks until < 5 hits are left in the collection
      while(nh-countusedhit > (_minnsh/2)) {
        int   finalcount(0);
        float phi1(0);
        _hist1->Reset("ICSEM");

        if (_debug > 1) printf("Cluster : %2i\n",counter);
        // Fill the phi histogram with all the unused hits
        for(int i=0; i<nh; i++){
          int ind = ordchcol[i];
          const mu2e::ComboHit* ch = &_chcol->at(ind);
          if(usedhit[i]==true) continue;
          float phi = ch->phi();
          // change the range to [0,2pi]
          if(phi < 0) phi += 2*M_PI;
          _hist1->Fill(phi,ch->nStrawHits());
          // if(_debug>3) std::cout<<"Phi "<<phi<<" ch phi "<< phi << "n straw hits = "<<ch->nStrawHits()<<std::endl;
        }

        // Find the min and max phi around the highest phi bin in the phi spectrum
        float cluphimin(0),cluphimax(0);
        clusterminmax(cluphimin, cluphimax);
        if (_debug > 2) printf("Min phi = %f Max phi = %f \n",cluphimin,cluphimax);
        float bin = _hist1->GetBinWidth(0);

        // Simple Case : When min phi numerically lower than max phi
        if(cluphimax > cluphimin and ((cluphimax-cluphimin) >= (bin+bin))){
          for(int i=0; i<nh; i++){
            const mu2e::ComboHit* ch = &_chcol->at(ordchcol[i]);
            float phi = ch->phi();
            if (phi < 0) phi += 2*M_PI;
            // if the phi is between min and max add it to the cluster
            if(phi>=cluphimin and phi<=cluphimax){
              finalcount++;
              usedhit[i] = true;
              countusedhit++;
              // Increment the sum phi of the cluster
              phi1 = phi1+phi;
              count1++;
              // Insert the cluster no. for each hit
              _clustno[i] = counter;
              // if(_debug>3) std::cout<<"Simple cluster case = "<<_clustno[i]<<" Cluster no. = "<<counter<<" hit = "<<i<<" Pos x = "<<ch->pos().X()<<" y = "<<ch->pos().Y()<<std::endl;
            }
          }
        }
        //  Alternate case when the hits are around 0 and 2PI
        else if(cluphimax < cluphimin){
          for(int i=0; i<nh; i++){
            const mu2e::ComboHit* ch = &_chcol->at(ordchcol[i]);
            float phi = ch->phi();
            if(phi < 0) phi += 2*M_PI;
            if((phi>=cluphimin and phi< _phimax) or (phi<=cluphimax and phi>=0)){
              finalcount++;
              if(phi > M_PI) phi = phi - 2*M_PI;
              usedhit[i] = true;
              countusedhit++;
              count2++;
              phi1 = phi1+phi;
              _clustno[i] = counter;
              // if(_debug>3)std::cout<<"Alternate cluster case = "<<_clustno[i]<<" Cluster no. = "<<counter<<" hit = "<<i<<" Pos x = "<<ch->pos().X()<<" y = "<<ch->pos().Y()<<std::endl;
           }
         }
        }
        else {
          // Hit associated to no cluster. Note : Need to investigate further
          // printf("Min phi = %f Max phi = %f \n",cluphimin,cluphimax);
          _clustno.clear();
          break;
        }
        //Total phi and no. of hits in a cluster
        if (phi1 !=0) {
          _phitotal[counter] = phi1;
          _nhitotal[counter] = finalcount;
        }
        // Increment the cluster number
        counter += 1;
      }
    }
    // phi separation between the clusters. Note : Only used for the events with two phi clusters at the moment
    if (counter == 2) {
      dphi = fabs(_phitotal[0]/_nhitotal[0]-_phitotal[1]/_nhitotal[1]);
      if (dphi > M_PI) dphi = 2*M_PI-dphi;
    }
    // Events where two phi clusters are found but they are separated < min delta phi.
    if (counter == 2 and dphi < _mindphi){
      TimeCluster otc;
      fillTimeSpectrumAlt(ordchcol);
      findPeaks(otc, ordchcol);
      assignHits(otc, ordchcol);
      initCluster(otc);
      if(otc._nsh > _minnsh) tccol1.push_back(otc);
      if (_debug > 1) {
        int onsh = otc._nsh;
        printf("dphi = %f nh = %i nsh = %i \n",dphi,nh,onsh);
      }
    }
    // Loop through the phi clusters and form time clusters
    for(int j=0;j<counter;j++){
       TimeCluster tc;
       fillTimeSpectrum(j ,ordchcol);
       findPeaks       (tc ,ordchcol);
       for(int ih=0; ih<nh; ih++) {
        //Fill the straw hit indices if the cluster number of the hit == j
        if (_clustno[ih]==j) tc._strawHitIdxs.push_back(ordchcol[ih]);
       }
       initCluster(tc);
       if (_debug > 1) {
         int nsh = tc._nsh;
         printf("Time cluster : %i No. straw hits : %i \n",counter,nsh);
       }
       // Sigma of the phi spectrum
       float sigma(0);
       if (tc._nsh > _minnsh) {
         if (tc._nsh > _nphiclusters) {
           sigma = checkdelta(tc, j, ordchcol);
           // if (_debug>2 and sigma>0) std::cout<<"Phi cluster sigma = "<<sigma<<" n straw hits = "<<tc._nsh<<" n combo hits = "<<tc._strawHitIdxs.size()<<std::endl;
         }
         if (sigma == 0 or sigma > _minsigma) {
           if (counter == 2){
             if (dphi >= _mindphi)
               tccol1.push_back(tc);
           }
          else tccol1.push_back(tc);
          }
          // if(_debug>1) std::cout<<"Delta phi = "<<dphi<<" T0 = "<<tc._t0._t0<<" n straw hits = "<<tc._nsh<<"n combo hits = "<<tc._strawHitIdxs.size()<<std::endl;
       }
       // if(_debug>2) std::cout<<"No. of time clusters = "<<counter<<" T0 = "<<tc._t0._t0<<" n straw hits = "<<tc._nsh<<" n combo hits = "<<tc._strawHitIdxs.size()<<std::endl;
       _timespec.Reset();
    }
  }

//------------------------------------------------------------------------------
// Function to find the min and max phi in the phi spectrum
//-----------------------------------------------------------------------------
  void PhiClusterFinder::clusterminmax(float& cluphimin, float& cluphimax) {
    int   nbx = _hist1->GetNbinsX();
    float bin = _hist1->GetBinWidth(0);

    int nsteps(0);
    int max_bin  = _hist1->GetMaximumBin();
    // Check to the left of the peak bin
    int bincheck = max_bin;
    while(_hist1->GetBinContent(bincheck) >=_threshold) {
      // if(_debug>3) std::cout<<"bincheck = "<<_hist1->GetBinContent(bincheck)<<" bincontent = "<<bincheck<<"   "<<_hist1->GetBinCenter(bincheck)<<std::endl;
      if (bincheck > _threshold) bincheck--;
      else bincheck = nbx;
      nsteps++;
      if (nsteps >= nbx) {
        cluphimin = 0;
        cluphimax = 2*M_PI;
        if (_debug>2) printf("Phi min = %10.3f max : %10.3f\n",cluphimin,cluphimax);
        return;
      }
    }
    // Check to the right of the highest bin
    int bincheckr=max_bin;
    while(_hist1->GetBinContent(bincheckr)>=_threshold){
      // if(_debug>3) std::cout<<"bincheckr = "<<_hist1->GetBinContent(bincheckr)<<"bincontent = "<<bincheckr<<"  "<<_hist1->GetBinCenter(bincheckr)<<std::endl;
      if(bincheckr<nbx) bincheckr++;
      else bincheckr = 1;
    }
    cluphimax = _hist1->GetXaxis()->GetBinCenter(bincheckr)+bin/2;
    cluphimin = _hist1->GetXaxis()->GetBinCenter(bincheck )-bin/2;
    if (_debug>3) printf("Phi min = %10.3f max : %10.3f\n",cluphimin,cluphimax);
  }

//------------------------------------------------------------------------------
// Function to fill the time spectrum
//-----------------------------------------------------------------------------
  void PhiClusterFinder::fillTimeSpectrum(int ClustNo, const std::vector<StrawHitIndex>& ordchcol) {
    _timespec.Reset();
    for(unsigned istr=0; istr<ordchcol.size();++istr) {
      int i = istr;
      if(_clustno[i]==ClustNo){
        ComboHit const& ch = (*_chcol)[istr];
        float time = _t0calc.comboHitTime((*_chcol)[istr],_pitch);
        _timespec.Fill(time,ch.nStrawHits());
      }
    }
  }

//-----------------------------------------------------------------------------------
// Function to fill the time spectrum when the hit is not associated to a phi cluster
//-----------------------------------------------------------------------------------
  void PhiClusterFinder::fillTimeSpectrumAlt(const std::vector<StrawHitIndex>& ordchcol) {
    _timespec.Reset();
    for(unsigned istr=0; istr<ordchcol.size();++istr) {
      ComboHit const& ch = (*_chcol)[istr];
      float time = _t0calc.comboHitTime((*_chcol)[istr],_pitch);
      _timespec.Fill(time,ch.nStrawHits());
    }
  }


//------------------------------------------------------------------------------
// Function to find the peak bin the time spectrum
//-----------------------------------------------------------------------------
  void PhiClusterFinder::findPeaks(TimeCluster& tc, const std::vector<StrawHitIndex>& ordchcol) {
    int nbins = _timespec.GetNbinsX()+1;
    // loop over spectrum to find peaks
    std::vector<BinContent> bcv;
    for (int ibin=1;ibin < nbins; ++ibin)
      if (_timespec.GetBinContent(ibin) >= _ymin) bcv.push_back(make_pair(_timespec.GetBinContent(ibin),ibin));
    std::sort(bcv.begin(),bcv.end(),[](const BinContent& x, const BinContent& y){return x.first > y.first;});
    for (const auto& bc : bcv) {
      float nsh(0.0);
      float t0(0.0);
      for (int ibin = std::max(1,bc.second-_npeak);ibin < std::min(nbins,bc.second+_npeak+1); ++ibin) {
        nsh += _timespec.GetBinContent(ibin);
        t0 += _timespec.GetBinCenter(ibin)*_timespec.GetBinContent(ibin);
      }
      t0 /= nsh;
      if (nsh >= _minnsh){
        tc._t0 = TrkT0(t0,_tbin*0.5);
        tc._nsh = nsh;
      }
    }
    // if(_debug>1) std::cout<<"Find peak T0 =   "<<tc._t0._t0<<" n straw hits = "<<tc._nsh<<std::endl;
  }

//------------------------------------------------------------------------------
// Function to assign hits to the time cluster
//-----------------------------------------------------------------------------
  void PhiClusterFinder::assignHits(TimeCluster& otc, const std::vector<StrawHitIndex>& ordchcol) {
    // assign hits to the closest time peak
    for(unsigned istr=0; istr<ordchcol.size();++istr) {
      float time = _t0calc.comboHitTime((*_chcol)[istr],_pitch);
      float dt = fabs(time - otc._t0._t0);
      if (dt < _maxdt+otc._t0._t0err){
        otc._strawHitIdxs.push_back(istr);
      }
    }
  }

//------------------------------------------------------------------------------
// Function to create the time clusters
//-----------------------------------------------------------------------------
  void PhiClusterFinder::initCluster(TimeCluster& tc){
    // use medians to initialize robustly
    accumulator_set<float, stats<tag::min > > tmin;
    accumulator_set<float, stats<tag::max > > tmax;
    accumulator_set<float, stats<tag::weighted_median(with_p_square_quantile) >, float > tacc, xacc, yacc, zacc;
    // No. of combo hits in the time cluster
    int nstrs = tc._strawHitIdxs.size();
    tc._nsh = 0;

    for (int i=0; i<nstrs; i++) {
      int loc = tc._strawHitIdxs[i];
      const ComboHit* ch = &_chcol->at(loc);
      float htime = _t0calc.comboHitTime(*ch,_pitch);
      unsigned nsh = ch->nStrawHits();
      tc._nsh += nsh;
      const XYZVectorF& pos = ch->pos();
      float hwt = ch->nStrawHits();
      tmin(htime);
      tmax(htime);
      tacc(htime,weight=hwt);
      xacc(pos.x(),weight=hwt);
      yacc(pos.y(),weight=hwt);
      zacc(pos.z(),weight=hwt);
    }

    static float invsqrt12(1.0/sqrt(12.0));
    tc._t0._t0 = extract_result<tag::weighted_median>(tacc);
    tc._t0._t0err = ( boost::accumulators::extract::max(tmax)-boost::accumulators::extract::min(tmin))*invsqrt12/sqrt(nstrs);
    tc._pos = XYZVectorF(extract_result<tag::weighted_median>(xacc),
    extract_result<tag::weighted_median>(yacc),
    extract_result<tag::weighted_median>(zacc));
    // if (_debug > 2) std::cout<<"Init Cluster T0 = "<<tc._t0._t0<<" n straw hits = "<<tc._nsh<<" n combo hits = "<<tc._strawHitIdxs.size()<<std::endl;
  }

//------------------------------------------------------------------------------
// Function to find the sigma of phi and time clusters. Note : Does not work yet for the clusters where some hits are close to 0 and some close to 2pi
//-----------------------------------------------------------------------------
  float PhiClusterFinder::checkdelta(TimeCluster& tc,int ClustNo, const std::vector<StrawHitIndex>& ordchcol) {
    float meanphi(0),sigphi(0);//,sig(0);
    for(auto ish : tc._strawHitIdxs) {
      const ComboHit* ch = &_chcol->at(ish);
      float phi = ch->phi();
      if(phi < 0) phi += 2*M_PI;
      meanphi = meanphi + phi;
    }
    meanphi = meanphi/tc._strawHitIdxs.size();
    for(auto ish :tc._strawHitIdxs) {
      const ComboHit* ch = &_chcol->at(ish);
      float phi = ch->phi();
      if(phi < 0) phi += 2*M_PI;
      sigphi = sigphi + fabs(meanphi-phi);
    }
    sigphi = sigphi/(tc._strawHitIdxs.size()-1);
    // if(_debug>3) std::cout<<"Check Delta time sigma = "<<sig<<" phi sigma = "<<sigphi<<" t0 = "<<tc._t0._t0<<" n straw hits = "<<tc._nsh<<" n combo hits = "<<tc._strawHitIdxs.size()<<std::endl;
    return sigphi;
  }

//------------------------------------------------------------------------------
// Function to add calo cluster to the time clusters
//-----------------------------------------------------------------------------
  void PhiClusterFinder::addCaloClusters(TimeClusterCollection& tccolnew, const TimeClusterCollection& tccol) {
    for(size_t ipeak=0; ipeak<tccol.size(); ipeak++) {
      auto& tc = tccol[ipeak];
      if (tc.hasCaloCluster()) {
       // if(_debug>1) std::cout<<"TC calo x = "<<tc._caloCluster->cog3Vector().x()<<" y = "<<tc._caloCluster->cog3Vector().y()<<" time = "<<tc._caloCluster->time()<<" nhits = "<<tc.nhits()<<std::endl;
        for(size_t inewpeak = 0; inewpeak<tccolnew.size();inewpeak++){
         auto& newtc = tccolnew[inewpeak];
         float dphi = fabs(polyAtan2( newtc._pos.y(), newtc._pos.x()) - polyAtan2( tc._caloCluster->cog3Vector().y(), tc._caloCluster->cog3Vector().x()));
         if (dphi > M_PI) dphi = 2*M_PI-dphi;
         if (dphi < _mindphi) {
          newtc._caloCluster = tc._caloCluster;
          // if(_debug>0) std::cout<<"New tc calo x = "<<newtc._caloCluster->cog3Vector().x()<<" time = "<<newtc._caloCluster->time()<<" nhits = "<<newtc.nhits()<<std::endl;
         }
        }
      }
    }
  }

}

DEFINE_ART_MODULE(mu2e::PhiClusterFinder)
