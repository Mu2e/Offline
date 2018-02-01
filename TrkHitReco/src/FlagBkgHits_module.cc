// $Id: FlagBkgHits_module.cc, without diagnostics $
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
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"

#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/BkgCluster.hh"
#include "RecoDataProducts/inc/BkgQual.hh"

#include "TrkReco/inc/TLTClusterer.hh"
#include "TrkReco/inc/TNTClusterer.hh"
#include "TrkReco/inc/TNTBClusterer.hh"
#include "Mu2eUtilities/inc/MVATools.hh"

#include "CLHEP/Units/PhysicalConstants.h"
#include "TMath.h"
#include "TH1F.h"
#include "TMVA/Reader.h"
#include "Math/VectorUtil.h"

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
using namespace ROOT::Math::VectorUtil;

namespace mu2e 
{

  class FlagBkgHits : public art::EDProducer
  {
    public:
      enum clusterer { TwoLevelThreshold=1, TwoNiveauThreshold=2, TwoNiveauThresholdB=3};
      explicit FlagBkgHits(fhicl::ParameterSet const&);
      virtual ~FlagBkgHits();
      virtual void beginJob();
      virtual void produce(art::Event& event ); 

    private:

      int _debug;
      int _printfreq;
      art::InputTag _chtag;
      bool          _flagall;
      unsigned _minnhits;
      unsigned _minnstereo;
      unsigned _minnp;
      unsigned _maxisolated;
      bool _savebkg;
      const ComboHitCollection* _chcol;
      BkgClusterer* _clusterer;
      float        _cperr2;
      bool          _useMVA;
      float        _bkgMVAcut;
      bool	    _filter;
      MVATools      _bkgMVA;

      // use TMVA Reader until problems with MVATools is resolved
      mutable TMVA::Reader _reader;
      mutable std::vector<float> _mvavars;

      void classifyClusters(BkgClusterCollection& clusters,BkgQualCollection& cquals) const; 
      void fillBkgQual(const BkgCluster& cluster,BkgQual& cqual) const;
      void countHits(const BkgCluster& cluster, unsigned& nactive, unsigned& nstereo) const;
      void countPlanes(const BkgCluster& cluster,BkgQual& cqual) const;
  };

  FlagBkgHits::FlagBkgHits(const fhicl::ParameterSet& pset) :
    _debug(pset.get<int>(           "debugLevel",0)),
    _printfreq(pset.get<int>(       "printFrequency",101)),
    _chtag(pset.get<art::InputTag>("ComboHitCollection","MakeStereoHits")),
    _flagall(pset.get<bool>(        "FlagAllHits",true)), 
    _minnhits(pset.get<unsigned>(   "MinActiveHits",5)),
    _minnstereo(pset.get<unsigned>( "MinStereoHits",2)),
    _minnp(pset.get<unsigned>(      "MinNPlanes",4)),
    _maxisolated(pset.get<unsigned>("MaxIsolated",0)),
    _savebkg(pset.get<bool>("SaveBkgClusters",false)),
    _useMVA(pset.get<bool>(         "UseBkgMVA",true)),
    _bkgMVAcut(pset.get<float>(    "BkgMVACut",0.8)),
    _filter(pset.get<bool>(  "FilterOutput",true)),
    _bkgMVA(pset.get<fhicl::ParameterSet>("BkgMVA",fhicl::ParameterSet()))
  {
    if(_filter)
      produces<ComboHitCollection>();
    else
      produces<StrawHitFlagCollection>();
    if(_savebkg) {
      produces<BkgClusterCollection>();
      produces<BkgQualCollection>();
    }

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

    float cperr = pset.get<float>("ClusterPositionError",10.0); 
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

  FlagBkgHits::~FlagBkgHits(){}


  void FlagBkgHits::beginJob()
  {
    _clusterer->init();
    if(_useMVA) _bkgMVA.initMVA();
  }


  void FlagBkgHits::produce(art::Event& event )
  {
    unsigned iev=event.id().event();
    if(_debug > 0 && (iev%_printfreq)==0) std::cout<<"FlagBkgHits: event="<<iev<<std::endl;

    art::Handle<ComboHitCollection> comboHitHandle;
    event.getByLabel(_chtag,comboHitHandle);
    _chcol = comboHitHandle.product();
// the primary output is either a deep copy of selected inputs or a flag collection on those
    std::unique_ptr<ComboHitCollection> chcol;    
    std::unique_ptr<StrawHitFlagCollection> bkgfcol;
// intermediate results: keep these on the heap unless requested for diagnostics later
    BkgClusterCollection bkgccol;
    BkgQualCollection bkgqcol;
 
    if(_filter){
      chcol = std::unique_ptr<ComboHitCollection>( new ComboHitCollection);
      chcol->reserve(_chcol->size());
    } else {
      bkgfcol = std::unique_ptr<StrawHitFlagCollection>(new StrawHitFlagCollection(_chcol->size()));
    }
    // guess on the size.
    bkgccol.reserve(_chcol->size()/2);
    // find clusters
    _clusterer->findClusters(bkgccol,*_chcol);
    if(_savebkg)bkgqcol.reserve(bkgccol.size());
    for (auto& cluster : bkgccol) {  
      BkgQual cqual;
      fillBkgQual(cluster,cqual);
      if(_savebkg)bkgqcol.push_back(cqual);
      if(cqual.MVAOutput() < _bkgMVAcut && cluster.hits().size() > _maxisolated ){
	if(_filter){
	  for (auto const& chit : cluster.hits())
	    chcol->push_back((*_chcol)[chit.index()]);
	}
      } else {
	StrawHitFlag flag;
	if(cluster.hits().size() <= _maxisolated){
	  cluster._flag.merge(BkgClusterFlag::iso);
	  flag.merge(StrawHitFlag::isolated);      
	}
	if(cqual.MVAOutput() > _bkgMVAcut){
	  cluster._flag.merge(BkgClusterFlag::bkg);
	  flag.merge(StrawHitFlag::bkgclust);
	}
	if(!_filter){
	  for (auto const& chit : cluster.hits())
	    (*bkgfcol)[chit.index()] = flag;
	}
      }
    }

    if(_filter)
      event.put(std::move(chcol));
    else
      event.put(std::move(bkgfcol));

    if(_savebkg){
      event.put(std::unique_ptr<BkgClusterCollection>(new BkgClusterCollection(bkgccol)));
      event.put(std::unique_ptr<BkgQualCollection>(new BkgQualCollection(bkgqcol)));
    }
  }

  void FlagBkgHits::fillBkgQual(const BkgCluster& cluster, BkgQual& cqual) const {

    unsigned nactive, nstereo;
    countHits(cluster,nactive,nstereo);
    cqual.setMVAStatus(BkgQual::unset);

    if (nactive >= _minnhits && nstereo >= _minnstereo)
    {
      cqual[BkgQual::nhits] = nactive;
      cqual[BkgQual::sfrac] = static_cast<float>(nstereo)/nactive;
      cqual[BkgQual::crho] = sqrtf(cluster.pos().perp2());

      countPlanes(cluster,cqual);
      if (cqual[BkgQual::np] >= _minnp)
      {
	accumulator_set<float, stats<tag::weighted_variance>, float> racc;
	accumulator_set<float, stats<tag::variance(lazy)> > tacc;
	std::vector<float> hz;
	for (const auto& chit : cluster.hits())
	{
	  if (chit.flag().hasAllProperties(StrawHitFlag::active))
	  {
	    const ComboHit& ch = (*_chcol)[chit.index()];
	    hz.push_back(ch.pos().z());
	    float dt = ch.time() - cluster.time();
	    tacc(dt);
	    // compute the transverse radius of this hit WRT the cluster center
	    XYZVec psep = PerpVector(ch.pos()-cluster.pos(),Geom::ZDir());
	    float rho = sqrtf(psep.mag2());
	    // project the hit position error along the radial direction to compute the weight
	    XYZVec pdir = psep.unit();
	    XYZVec tdir(-ch.wdir().y(),ch.wdir().x(),0.0);
	    float rwerr = ch.posRes(ComboHit::wire)*pdir.Dot(ch.wdir());
	    float rterr = ch.posRes(ComboHit::trans)*pdir.Dot(tdir);
	    // include the cluster center position error when computing the weight
	    float rwt = 1.0/sqrtf(rwerr*rwerr + rterr*rterr + _cperr2/nactive);
	    racc(rho,weight=rwt);
	  }
	}
	cqual[BkgQual::hrho] = extract_result<tag::weighted_mean>(racc);
	cqual[BkgQual::shrho] = sqrtf(std::max(extract_result<tag::weighted_variance>(racc),float(0.0)));
	cqual[BkgQual::sdt] = sqrtf(std::max(extract_result<tag::variance>(tacc),float(0.0)));

	// find the min, max and gap from the sorted Z positions
	std::sort(hz.begin(),hz.end());
	cqual[BkgQual::zmin] = hz.front();
	cqual[BkgQual::zmax] = hz.back();

	// find biggest Z gap
	float zgap = 0.0;
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
	  float readerout = _reader.EvaluateMVA( "MLP method" );
	  // std::vector<float> mvavars(9,0.0);
	  // for(size_t ivar =0;ivar<_mvavars.size(); ++ivar)
	  // mvavars[ivar] = _mvavars[ivar];
	  // float toolout = _bkgMVA.evalMVA(mvavars);
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

  void FlagBkgHits::countPlanes(const BkgCluster& cluster, BkgQual& cqual ) const 
  {
    std::array<int,StrawId::_nplanes> hitplanes{0};
    for (const auto& chit : cluster.hits()) 
    {
      if (chit.flag().hasAllProperties(StrawHitFlag::active))
      {
	const ComboHit& ch = (*_chcol)[chit.index()];
	hitplanes[ch.sid().plane()] += ch.nStrawHits();
      }
    }

    unsigned ipmin = 0;
    unsigned ipmax = StrawId::_nplanes-1; 
    while (hitplanes[ipmin]==0) ++ipmin;
    while (hitplanes[ipmax]==0) --ipmax;

    unsigned npexp(0);
    unsigned np(0);
    unsigned nphits(0);

    for(unsigned ip = ipmin; ip <= ipmax; ++ip)
    {
      npexp++; // should use TTracker to see if plane is physically present FIXME!
      if (hitplanes[ip]> 0)++np;
      nphits += hitplanes[ip];
    }
    cqual[BkgQual::np] = np; 
    cqual[BkgQual::npexp] = npexp; 
    cqual[BkgQual::npfrac] = np/static_cast<float>(npexp); 
    cqual[BkgQual::nphits] = nphits/static_cast<float>(np);
  }

  void FlagBkgHits::countHits(const BkgCluster& cluster, unsigned& nactive, unsigned& nstereo) const 
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

using mu2e::FlagBkgHits;
DEFINE_ART_MODULE(FlagBkgHits);
