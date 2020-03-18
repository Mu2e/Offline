//
// $Id: LineFinder_module.cc,v 1.2 2014/08/30 12:19:38 tassiell Exp $
// $Author: tassiell $
// $Date: 2014/08/30 12:19:38 $
//
// Original author D. Brown and G. Tassielli
//

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Persistency/Common/Ptr.h"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"

#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/LineSeed.hh"

#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include "TH2F.h"

#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <utility>
#include <functional>
#include <float.h>
#include <vector>
#include <map>

namespace mu2e{

  class LineFinder : public art::EDProducer {
    public:
      struct Config{
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;
        fhicl::Atom<int> diag{Name("diag"), Comment("Create diag histograms"),0};
        fhicl::Atom<int> minpeak{Name("minPeak"), Comment("Minimum hits in accumulator peak"),3};
        fhicl::Atom<float> maxDOCA{Name("maxDOCA"), Comment("Largest DOCA that is considered on the track (mm)"),3};
        fhicl::Atom<float> t0offset{Name("t0offset"), Comment("T0 offset"), 0};
        fhicl::Atom<int> nsteps{Name("NSteps"), Comment("Number of steps per straw"), 8};
        fhicl::Atom<float> stepsize{Name("StepSize"), Comment("Size of each step in fraction of res"), 0.5};
        fhicl::Atom<art::InputTag> shToken{Name("ComboHitCollection"),Comment("tag for straw hit collection")};
        fhicl::Atom<art::InputTag> tcToken{Name("TimeClusterCollection"),Comment("tag for time cluster collection")};
      };
      typedef art::EDProducer::Table<Config> Parameters;
      explicit LineFinder(const Parameters& conf);
      virtual ~LineFinder(){};
      virtual void produce(art::Event& event ) override;

    private:

      Config _conf;

      //config parameters:
      int _diag;
      int  _minPeak;
      float _maxDOCA;
      float _t0offset;
      int _Nsteps;
      float _stepSize;
      art::InputTag  _shToken;
      art::InputTag  _tcToken;

      void findLine(const ComboHitCollection& shC, LineSeed &lseed);
  };


 LineFinder::LineFinder(const Parameters& conf) :
   art::EDProducer(conf),
   	_diag (conf().diag()),
        _minPeak (conf().minpeak()),
        _maxDOCA (conf().maxDOCA()),
        _t0offset (conf().t0offset()),
        _Nsteps (conf().nsteps()),
        _stepSize (conf().stepsize()),
    	_shToken (conf().shToken()),
	_tcToken (conf().tcToken())
{
  consumes<ComboHitCollection>(_shToken);
  consumes<TimeClusterCollection>(_tcToken);
  produces<LineSeedCollection>();
	    
 }

void LineFinder::produce(art::Event& event ) {

  auto const& shH = event.getValidHandle<ComboHitCollection>(_shToken);
  const ComboHitCollection& shcol(*shH);
  auto  const& tcH = event.getValidHandle<TimeClusterCollection>(_tcToken);
  const TimeClusterCollection& tccol(*tcH);



  std::unique_ptr<LineSeedCollection> seed_col(new LineSeedCollection());

  for (size_t index=0;index< tccol.size();++index) {
    const auto& tclust = tccol[index];

    const std::vector<StrawHitIndex>& shIndices = tclust.hits();
    ComboHitCollection tchits;
    
    int nhits = 0;
    for (size_t i=0; i<shIndices.size(); ++i) {
      int loc = shIndices[i];
      const ComboHit& sh  = shcol[loc];
      tchits.push_back(ComboHit(sh));
      nhits += sh.nStrawHits();
    }

    LineSeed lseed;
    lseed._timeCluster = art::Ptr<TimeCluster>(tcH,index);
    lseed._converged = 1;
    findLine(tchits, lseed);

    if (lseed._seedSize >= _minPeak){
      LineSeedCollection*  lcol  = seed_col.get();
      lcol->push_back(lseed);
    }
  }

  event.put(std::move(seed_col));    
}

void LineFinder::findLine(const ComboHitCollection& shC, LineSeed& lseed){

  mu2e::GeomHandle<mu2e::Tracker> th;
  auto tracker = th.get();
  int bestcount = 0;
  double bestll = 0;
  
 
  // lets get the best pairwise vector
  for (size_t i=0;i<shC.size();i++){
    Straw const& strawi = tracker->getStraw(shC[i].strawId());
    for (size_t j=i+1;j<shC.size();j++){
      Straw const& strawj = tracker->getStraw(shC[j].strawId());
      for (int is=-1*_Nsteps;is<_Nsteps+1;is++){
        CLHEP::Hep3Vector ipos = shC[i].posCLHEP() + strawi.getDirection()*shC[i].wireRes()*_stepSize*is;
        for (int js=-1*_Nsteps;js<_Nsteps+1;js++){
          CLHEP::Hep3Vector jpos = shC[j].posCLHEP() + strawj.getDirection()*shC[j].wireRes()*_stepSize*js;

          CLHEP::Hep3Vector newdir = (jpos-ipos).unit();
          // now loop over all hits and see how many are in this track
          int count = 0;
          double ll = 0;
          for (size_t k=0;k<shC.size();k++){
            Straw const& strawk = tracker->getStraw(shC[k].strawId());
            TwoLinePCA pca( strawk.getMidPoint(), strawk.getDirection(),
                ipos, newdir);
            double dist = (pca.point1()-strawk.getMidPoint()).mag();
            if (pca.dca() < _maxDOCA && dist < strawk.halfLength()){
              count += 1;
              ll += pow(dist-shC[k].wireDist(),2)/shC[k].wireErr2();
            }
          }
          if (count > bestcount || (count == bestcount && ll < bestll)){
            bestcount = count;
            bestll = ll;
            lseed._seedDir = newdir.unit();
            if (lseed._seedDir.y() > 0) lseed._seedDir *= -1;
            lseed._seedInt = ipos - newdir*ipos.y()/newdir.y();
          }
        }
      }
    }
  }

  double avg_t0 = 0;
  int good_hits = 0;
  for (size_t k=0;k<shC.size();k++){
    Straw const& strawk = tracker->getStraw(shC[k].strawId());
    TwoLinePCA pca( strawk.getMidPoint(), strawk.getDirection(),
        lseed._seedInt, lseed._seedDir);
    double dist = (pca.point1()-strawk.getMidPoint()).dot(strawk.getDirection());
    if (pca.dca() < _maxDOCA && fabs(dist) < strawk.halfLength()){
      double traj_time = ((pca.point2() - lseed._seedInt).dot(lseed._seedDir))/299.9;
      double hit_t0 = shC[k].time() - shC[k].driftTime() - shC[k].propTime() - traj_time;
      avg_t0 += hit_t0;
      good_hits++;
      lseed._strawHitIdxs.push_back(k);
    }
  }
  avg_t0 /= good_hits;


  lseed._seedSize = good_hits;
  lseed._t0 = avg_t0-_t0offset;
}

}//end mu2e namespace
using mu2e::LineFinder;
DEFINE_ART_MODULE(LineFinder);
