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
#include "RecoDataProducts/inc/CosmicTrackSeed.hh"

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
        fhicl::Atom<art::InputTag> chToken{Name("ComboHitCollection"),Comment("tag for straw hit collection")};
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
      art::InputTag  _chToken;
      art::InputTag  _tcToken;

      int findLine(const ComboHitCollection& shC, art::Event const& event, CosmicTrackSeed &tseed);
  };


 LineFinder::LineFinder(const Parameters& conf) :
   art::EDProducer(conf),
   	_diag (conf().diag()),
        _minPeak (conf().minpeak()),
        _maxDOCA (conf().maxDOCA()),
        _t0offset (conf().t0offset()),
        _Nsteps (conf().nsteps()),
        _stepSize (conf().stepsize()),
    	_chToken (conf().chToken()),
	_tcToken (conf().tcToken())
{
  consumes<ComboHitCollection>(_chToken);
  consumes<TimeClusterCollection>(_tcToken);
  produces<CosmicTrackSeedCollection>();
	    
 }

void LineFinder::produce(art::Event& event ) {

  auto const& chH = event.getValidHandle<ComboHitCollection>(_chToken);
  const ComboHitCollection& chcol(*chH);
  auto  const& tcH = event.getValidHandle<TimeClusterCollection>(_tcToken);
  const TimeClusterCollection& tccol(*tcH);



  std::unique_ptr<CosmicTrackSeedCollection> seed_col(new CosmicTrackSeedCollection());

  for (size_t index=0;index< tccol.size();++index) {
    const auto& tclust = tccol[index];

    ComboHitCollection tchits;
    
    int nhits = 0;
    std::vector<ComboHitCollection::const_iterator> chids;  
    chcol.fillComboHits(event, tclust.hits(), chids); 
    for (auto const& it : chids){
      tchits.push_back(it[0]);
      nhits += it[0].nStrawHits();
    }

    CosmicTrackSeed tseed;
    tseed._timeCluster = art::Ptr<TimeCluster>(tcH,index);
    tseed._track.converged = true;

    int seedSize = findLine(tchits, event, tseed);

    if (seedSize >= _minPeak){
      tseed._status.merge(TrkFitFlag::Straight);
      tseed._status.merge(TrkFitFlag::hitsOK);
      tseed._status.merge(TrkFitFlag::helixOK);
      tseed._status.merge(TrkFitFlag::helixConverged);

      CosmicTrackSeedCollection*  tcol  = seed_col.get();
      tcol->push_back(tseed);
    }
  }

  event.put(std::move(seed_col));    
}

int LineFinder::findLine(const ComboHitCollection& shC, art::Event const& event, CosmicTrackSeed& tseed){

  mu2e::GeomHandle<mu2e::Tracker> th;
  auto tracker = th.get();
  int bestcount = 0;
  double bestll = 0;
  
  CLHEP::Hep3Vector seedDir(0,0,0);
  CLHEP::Hep3Vector seedInt(0,0,0);
 
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
            seedDir = newdir.unit();
            if (seedDir.y() > 0) seedDir *= -1;
            seedInt = ipos - newdir*ipos.y()/newdir.y();
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
        seedInt, seedDir);
    double dist = (pca.point1()-strawk.getMidPoint()).dot(strawk.getDirection());
    if (pca.dca() < _maxDOCA && fabs(dist) < strawk.halfLength()){
      double traj_time = ((pca.point2() - seedInt).dot(seedDir))/299.9;
      double hit_t0 = shC[k].time() - shC[k].driftTime() - shC[k].propTime() - traj_time;
      avg_t0 += hit_t0;
      good_hits++;
      tseed._straw_chits.push_back(shC[k]);
    }
  }
  avg_t0 /= good_hits;

  tseed._t0._t0 = avg_t0-_t0offset;

  // get pos and direction into Z alignment
  if (seedDir.z() != 0){
    seedDir /= seedDir.z();
    seedInt -= seedDir*seedInt.z()/seedDir.z();
  }
  
  tseed._track.FitParams.A0 = seedInt.x(); 
  tseed._track.FitParams.B0 = seedInt.y(); 
  tseed._track.FitParams.A1 = seedDir.x(); 
  tseed._track.FitParams.B1 = seedDir.y(); 
  XYZVec X(1,0,0);
  XYZVec Y(0,1,0);
  XYZVec Z(0,0,1);
  TrackAxes XYZ(X,Y,Z);
  tseed._track.InitCoordSystem = XYZ; 
  tseed._track.FitCoordSystem = XYZ; 
  TrackEquation XYZTrack(Geom::toXYZVec(seedInt), Geom::toXYZVec(seedDir));
  tseed._track.SetFitEquation(XYZTrack);

  // For compatibility FIXME
  for(size_t ich= 0; ich<tseed._straw_chits.size(); ich++){  
    std::vector<StrawHitIndex> shitids;          	          		
    tseed._straw_chits.fillStrawHitIndices(event, ich, shitids);  

    for(auto const& ids : shitids){ 
      size_t    istraw   = (ids);
      TrkStrawHitSeed tshs;
      tshs._index  = istraw;
      tshs._t0 = tseed._t0;
      tseed._trkstrawhits.push_back(tshs); 
    }  
  }

  return good_hits;
}

}//end mu2e namespace
using mu2e::LineFinder;
DEFINE_ART_MODULE(LineFinder);
