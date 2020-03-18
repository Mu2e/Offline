// Author : R Bonventre
// Data : Jan 2020
// Purpose: Cosmic Track Finder with t0 and drift fit

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Persistency/Common/Ptr.h"

#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/LineSeed.hh"
#include "CosmicReco/inc/PDFFit.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "DataProducts/inc/XYZVec.hh"
#include "Mu2eUtilities/inc/ParametricFit.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
#include "ProditionsService/inc/ProditionsHandle.hh"

#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <utility>
#include <functional>
#include <float.h>
#include <vector>
#include <map>

#include <Minuit2/FCNBase.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnSimplex.h>
#include <Minuit2/MnMinos.h>
#include <Minuit2/MnStrategy.h>
#include <Minuit2/MnPrint.h>
#include <Minuit2/MnUserParameters.h>
#include <Minuit2/FunctionMinimum.h>

namespace mu2e{

  class CosmicDriftFit : public art::EDProducer {
    public:
      struct Config{
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;
        fhicl::Atom<int> diag{Name("diag"), Comment("Create diag histograms"),0};
        fhicl::Atom<float> lerror{Name("Lerror"), Comment("Starting intercept error for gaussian fit"), 40};
        fhicl::Atom<float> terror{Name("Terror"), Comment("Starting direction error for gaussian fit"), 2.5};
        fhicl::Atom<art::InputTag> shToken{Name("ComboHitCollection"),Comment("tag for combo hit collection")};
        fhicl::Atom<art::InputTag> lsToken{Name("LineSeedCollection"),Comment("tag for line seed collection")};
      };
      typedef art::EDProducer::Table<Config> Parameters;
      explicit CosmicDriftFit(const Parameters& conf);
      virtual ~CosmicDriftFit(){};
      virtual void produce(art::Event& event ) override;

    private:

      Config _conf;

      //config parameters:
      int _diag;
      float _lerror;
      float _terror;
      art::InputTag  _shToken;
      art::InputTag  _lsToken;

  };


 CosmicDriftFit::CosmicDriftFit(const Parameters& conf) :
   art::EDProducer(conf),
   	_diag (conf().diag()),
        _lerror (conf().lerror()),
        _terror (conf().terror()),
    	_shToken (conf().shToken()),
        _lsToken (conf().lsToken())
{
  consumes<ComboHitCollection>(_shToken);
  consumes<LineSeedCollection>(_lsToken);
  produces<LineSeedCollection>(); 
 }

void CosmicDriftFit::produce(art::Event& event ) {

  auto const& chH = event.getValidHandle<ComboHitCollection>(_shToken);
  const ComboHitCollection& shcol(*chH);
  auto  const& lsH = event.getValidHandle<LineSeedCollection>(_lsToken);
  const LineSeedCollection& lscol(*lsH);

  mu2e::GeomHandle<mu2e::Tracker> th;
  const Tracker* tracker = th.get();

  ProditionsHandle<StrawResponse> _strawResponse_h;
  auto _srep = _strawResponse_h.getPtr(event.id());
  StrawResponse const& srep = * _srep.get();

 

  std::unique_ptr<LineSeedCollection> seed_col(new LineSeedCollection());

  for (size_t index=0;index<lscol.size();++index) {
    const auto& lseed = lscol[index];
    auto tclust = (*lseed._timeCluster);

    const std::vector<StrawHitIndex>& shIndices = lseed._strawHitIdxs;
    ComboHitCollection tchits;

    for (size_t i=0; i<shIndices.size(); ++i) {
      int loc = shIndices[i];
      const ComboHit& sh  = shcol[loc];
      //FIXME add flag check
      //if(ch.flag().hasAnyProperty(_hsel) && !ch.flag().hasAnyProperty(_hbkg) ) {
      tchits.push_back(ComboHit(sh));
    }

    CLHEP::Hep3Vector dir = lseed._seedDir;
    CLHEP::Hep3Vector intercept = lseed._seedInt;
    dir /= -1*dir.y();
    intercept -= dir*intercept.y()/dir.y();
  
    // now gaussian fit, transverse distance only
    std::vector<double> errors(5,0);
    std::vector<double> seed(5,0);
    seed[0] = intercept.x();
    seed[1] = intercept.z();
    seed[2] = dir.x();
    seed[3] = dir.z();
    seed[4] = lseed._t0;

    errors[0] = _lerror;
    errors[1] = _lerror;
    errors[2] = _terror;
    errors[3] = _terror;
    errors[4] = 1;
	  
    //Define the PDF used by Minuit:
    GaussianDriftFit fit(tchits, srep, tracker);
	  
    //Initiate Minuit Fit:
    ROOT::Minuit2::MnStrategy mnStrategy(2); 
    ROOT::Minuit2::MnUserParameters params(seed,errors);
    ROOT::Minuit2::MnMigrad migrad(fit,params,mnStrategy);
//    migrad.Fix(2);
//    migrad.Fix(3);
//    migrad.Fix(4);
    //Define Minimization method as "MIGRAD" (see minuit documentation)
    ROOT::Minuit2::FunctionMinimum min = migrad(0, 0.1);
    if(_diag > 1){
      ROOT::Minuit2::MnPrint::SetLevel(3);
      ROOT::Minuit2::operator<<(cout, min);
    }else{
      ROOT::Minuit2::MnPrint::SetLevel(0);
    }
    //Will be the results of the fit routine:
    ROOT::Minuit2::MnUserParameters results = min.UserParameters();
    double minval = min.Fval();

    LineSeed newseed;
    newseed._seedSize = lseed._seedSize;
    newseed._timeCluster = lseed._timeCluster;
    newseed._strawHitIdxs = lseed._strawHitIdxs;

    if (min.IsValid()){
      std::cout << "CONVERGED " << minval << std::endl;
      newseed._converged = 1;
    }else{
      std::cout << "DID NOT CONVERGE " << minval << std::endl;
      newseed._converged = 0;
    }
    fit.setResults(newseed,results.Params());
        
    LineSeedCollection*  lcol  = seed_col.get();
    lcol->push_back(newseed);
  }

  event.put(std::move(seed_col));    
}

}//end mu2e namespace
using mu2e::CosmicDriftFit;
DEFINE_ART_MODULE(CosmicDriftFit);
