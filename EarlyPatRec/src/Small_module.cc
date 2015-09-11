//=============================================================================
//
// Plugin to look at the straws in the TTracker.
// This is just a temporary tool to help learn how to write the
// PatRec geometry understander.
//
// $Id: Small_module.cc,v 1.14 2013/10/21 20:34:14 gandr Exp $
// $Author: gandr $
// $Date: 2013/10/21 20:34:14 $
//
// Original author: Mark Fischler
//
//=============================================================================
// C++ includes
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <cassert>

#include "CLHEP/Vector/TwoVector.h"

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Provenance.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Root includes.
#include "TFile.h"
#include "TH1F.h"

// Mu2e includes.
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "GeneralUtilities/inc/LineSegmentPCA.hh"
#include "Mu2eUtilities/inc/SimParticlesWithHits.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "RecoDataProducts/inc/StrawClusterCollection.hh"
#include "TrackerGeom/inc/Tracker.hh"

using namespace std;

namespace mu2e {;
  enum PrintLevel { quiet  =-1,
                    normal = 0,
                    verbose= 1};

 class pstraw{
   //
   // pseudo straw class
   //
 public:
   Int_t   lay;
   Int_t   did;
   Int_t   sec;
   Float_t hl;
   Float_t mpx;
   Float_t mpy;
   Float_t mpz;
   Float_t dirx;
   Float_t diry;
   Float_t dirz; // should always be 0
   /*
    bool operator>(const pstraw other) const {
      if (id > other.id) {
        return true;
      }
      else{
        return false;
      }
    }
   bool operator<(const pstraw other) const {
      if (id < other.id) {
        return true;
      }
      else{
        return false;
      }
   }
   bool operator==(const straw other) const {
      if (id == other.id) {
        return true;
      }
      else{
        return false;
      }
   }
   */
    void Print()
    {
      cout<< "Straw:  " << endl;
      cout<< "======  " << endl;
      cout<< "Layer:  " << lay  <<endl;
      cout<< "DID:    " << did  <<endl;
      cout<< "Sector: " << sec  <<endl;
      cout<< "hl:     " << hl   <<endl;
      cout<< "mpx:    " << mpx  <<endl;
      cout<< "mpy:    " << mpy  <<endl;
      cout<< "mpz:    " << mpz  <<endl;
      cout<< "dirx:   " << dirx <<endl;
      cout<< "diry:   " << diry <<endl;
      cout<< "dirz:   " << dirz <<endl;
    }
 };

  //--------------------------------------------------------------------
  //
  //
  class GrokGeometry : public art::EDAnalyzer {
  public:
    explicit GrokGeometry(fhicl::ParameterSet const& pset):
      art::EDAnalyzer(pset),
      _diagLevel(pset.get<int>("diagLevel",0)),
      _maxFullPrint(pset.get<int>("maxFullPrint",5)),
      _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
      _g4ModuleLabel(pset.get<string>("g4ModuleLabel")),
      _trackerStepPoints(pset.get<string>("trackerStepPoints","tracker")),
      _makerModuleLabel(pset.get<std::string>("makerModuleLabel")),
      _clmakerModuleLabel(pset.get<std::string>("clmakerModuleLabel")),
       _EnergyDep_s(0),
      _EnergyDepX_s(0),
      _beta_c(0),
      _tanTau_c(0),
      _tanTheta_c(0)
    {
    }
    virtual ~GrokGeometry() { }

    virtual void beginJob();

    void analyze( art::Event const& e);
  private:

    // Diagnostics level.
    int _diagLevel;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;
    // Module label of the geerator module.
    std::string _generatorModuleLabel;
    // Module label of the g4 module that made the hits.
    std::string _g4ModuleLabel;

    // Name of the tracker StepPoint collection
    std::string _trackerStepPoints;
    // Label of the module that made the StrawHits.
    std::string _makerModuleLabel;
    // Label of the module that made the Clusters.
    std::string _clmakerModuleLabel;

    // Some diagnostic histograms.
    TH1F* _EnergyDep_s;
    TH1F* _EnergyDepX_s;
    TH1F* _beta_c;
    TH1F* _tanTau_c;
    TH1F* _tanTheta_c;
    //
  }; // end of GrokGeometry class definition


  void GrokGeometry::beginJob(){
    cout << "Diaglevel: "
         << _diagLevel << " "
         << _maxFullPrint<<endl;

    art::ServiceHandle<art::TFileService> tfs;
    _EnergyDep_s     = tfs->make<TH1F>( "EnergyDep_s", "s Energy Deposited in Straw (KeV)", 100, 0.0, 10.0);
    _EnergyDepX_s    = tfs->make<TH1F>( "EnergyDepX_s","s Energy Deposited in Straw (KeV)", 100, 0.0, 10.0);
    _beta_c          = tfs->make<TH1F>( "beta_c",      "s straw incidence angle beta", 45, 0., 90. );
    _tanTau_c        = tfs->make<TH1F>( "tanTau_c",    "s panel attack angle tan tau", 50, 0., 2.0 );
    _tanTheta_c      = tfs->make<TH1F>( "tanTheta_c",  "s pitch tan theta", 50, 0., 2. );

   } // end of GrokGeometry beginJob

  void GrokGeometry::analyze(art::Event const& evt)
  {
    if ( _diagLevel > 2 ) cout << "GrokGeometry: analyze() begin"<<endl;
    // Geometry info for the TTracker.
    // Get a reference to one of the T trackers.
    // Throw exception if not successful.
    const Tracker& tracker = getTrackerOrThrow();
    //
    //
    // Position of the center of the tracker in mu2e coordinates.
    //
    //CLHEP::Hep2Vector tt =CLHEP::Hep2Vector( 0.0, 10200.);
    //CLHEP::Hep3Vector point =CLHEP::Hep3Vector( -3904, 0.0, 10200.);
    //CLHEP::Hep3Vector bf = bfMgr->getBField(point);
    //cout << " B-field: (getBField(center of tracker)) " <<bf<<endl;
    static int ncalls(0);
    ++ncalls;
    StrawId nsid;
    Straw str;
    StrawId sid;
    LayerId lid;
    DeviceId did;
    SectorId secid;
    int sector;

    vector<CLHEP::Hep3Vector> momentum_cluster(36);
    vector<CLHEP::Hep3Vector> deviceStrawDirections(36);
    vector<double> beta_cluster;            // angle against projection of wire
    vector<double> tanTau_cluster;          // attack angle of path to panel
    vector<double> tanTheta_cluster;        // helix pitch

    //double  edep[36] ;
    int nhitdev[36];
    CLHEP::Hep3Vector  MCPoint[36];

    tanTheta_cluster.clear();     // helix pitch
    beta_cluster.clear();         // angle against projection of wire
    tanTau_cluster.clear();       // attack angle of path to panel

    multimap<int,pstraw> mpstraws;
    mpstraws.clear();

    //
    // Get the straw clusters from the event
    art::Handle<StrawClusterCollection> clustersHandle;
    evt.getByLabel(_clmakerModuleLabel,clustersHandle);
    StrawClusterCollection const& clusters = *clustersHandle;

    // Get handles to the generated and simulated particles.
    art::Handle<GenParticleCollection> genParticles;
    evt.getByLabel(_generatorModuleLabel, genParticles);

    art::Handle<SimParticleCollection> simParticles;
    evt.getByLabel(_g4ModuleLabel, simParticles);

    // Handle to information about G4 physical volumes.
    art::Handle<PhysicalVolumeInfoCollection> volumes;
    evt.getRun().getByLabel(_g4ModuleLabel, volumes);

    //Some files might not have the SimParticle and volume information.
    bool haveSimPart = ( simParticles.isValid() && volumes.isValid() );

    // Other files might have empty collections.
    if ( haveSimPart ){
      haveSimPart = !(simParticles->empty() || volumes->empty());
    }

//cout << "[[  3 ]]\n";
    // Construct an object that ties together all of the simulated particle and hit info.
    SimParticlesWithHits sims( evt,
                               _g4ModuleLabel,
                               _makerModuleLabel,
                               "tracker",
                               0.001,
                               5 );
    if (sims.size()<1) return;  // no sim particles found
    bool foundcele=false;
    typedef SimParticlesWithHits::map_type map_type;
    for ( map_type::const_iterator i=sims.begin();
          i != sims.end(); ++i )                      // loop over simparticles
      {
        // All information about this SimParticle
        SimParticleInfo const& simInfo = i->second;
        // Information about StrawHits that belong on this SimParticle.
        vector<StrawHitMCInfo> const& infos = simInfo.strawHitInfos();
        did = -1;
        if (simInfo.simParticle().generatorIndex()>=0)
          {
            const GenParticle genpar  =genParticles->at(simInfo.simParticle().generatorIndex());

            //cout<< genpar.generatorId()<<endl;
            if (genpar.generatorId()== GenId::conversionGun)
              {
                foundcele=true;
                /*
                cout << "SimParticle associated to conversion electron: "
                     << " Event: " << evt.id().event()
                     << " Track: " << i->first
                     << " PdgId: " << simInfo.simParticle().pdgId()
                     << " |p|: "   << simInfo.simParticle().startMomentum().vect().mag()
                     << " Hits: "  << infos.size()
                     << " CC:   "  << simInfo.simParticle().creationCode()
                     << " GI:   "  << simInfo.simParticle().generatorIndex()
                     << endl;
                cout << "Polar: "<<simInfo.simParticle().startMomentum().vect().getTheta ()<<endl;
                */
                double _timetodist=149.8962;

                // calculate the average hit position of track at a plane
                for (int idev = 0; idev < 36 ; idev++) {
                  nhitdev[idev] = 0 ;
                  //edep[idev] = 0.0 ;
                  MCPoint[idev] = CLHEP::Hep3Vector(0.,0.,0.);
                }
                for ( size_t associatedHit=0; associatedHit<infos.size(); ++associatedHit) // Loop over associated Hits
                  {
                    StrawHitMCInfo const& info = infos.at(associatedHit);
                    StrawHit const& hit        = info.hit();
                    Straw const& str           = tracker.getStraw(hit.strawIndex());
                    sid = str.id();
                    did = sid.getDeviceId();
                    std::vector<StepPointMC const *> const& steps = info.steps();
                    // mf study 1
                    sector = sid.getSector();
                    const CLHEP::Hep3Vector straw_direction = str.getDirection();
                    if ((straw_direction.z() > .000001) || (straw_direction.mag() > 1.000001)) {
                      cout << "????? unexpected straw direction: \n      " << straw_direction << "\n";
                    }
                    if ((sector < 0) || (sector > 5)) {
                      cout << "????? unexpected straw sector: \n      " << sector << "\n";
                    }
                    // --- mf
                    deviceStrawDirections[did] = straw_direction;
                    //cout << "Device " << did << " Straw Direction " << straw_direction << "\n";
                    double energyAH  = 0.;
                    double energyAHX = 0.;
                    for ( size_t ks=0; ks<steps.size(); ++ks){
                      StepPointMC const& step = *(steps[ks]);
                      if (step.momentum().mag()>5)
                        {
                          MCPoint[did] =  MCPoint[did]+step.position();
                          // mf study 2
                          // cout << "Energy Step =  " << step.totalEDep() << "\n";
                          energyAH  += step.totalEDep();
                          energyAHX += step.totalEDep();
                          // mf study 1
                          momentum_cluster[did] +=step.momentum();
                          // --- mf
                          nhitdev[did]++;
                        }
                      else // step momentum is less than 5
                        {
                          energyAH  += step.totalEDep();
                          //                      cout << "delta: " << step.momentum().mag()<<endl;
                        } // end of if/else for momentum > 5
                    } // end of loops over steps in this hit
                    // mf study 2
                    // cout << "Energy for associated hit " << associatedHit << " =  " << energyAH << "\n";
                    _EnergyDep_s->Fill(1000.0*energyAH);
                    _EnergyDepX_s->Fill(1000.0*energyAHX);
                    // --- mf
                  } // end of loop over associated hits

                // Device quantity normalization loop
                for (int idev = 0; idev < 36 ; idev++) {
                  if (nhitdev[idev] <= 0) continue;
                  double a = 1.0/double(nhitdev[idev]);
                  MCPoint[idev] = MCPoint[idev]*a;
                  // mf study 1
                  momentum_cluster[idev] *= a;
                  CLHEP::Hep3Vector projectedMomentum = momentum_cluster[idev];
                  projectedMomentum.setZ(0);
                  double beta = std::acos(projectedMomentum.dot(deviceStrawDirections[idev])/projectedMomentum.mag());
                  // double beta = 0;
                  _beta_c->Fill(beta*180.0/3.141592653589793);
                  double tanTheta = momentum_cluster[idev].perp()/momentum_cluster[idev].z();
                  _tanTheta_c->Fill(tanTheta);
                  double tanTau = tanTheta*std::sin(beta);
                  _tanTau_c->Fill(tanTau);
                  // --- mf
                } // End of Device quantity normalization loop
//cout << "[[  4 ]]\n";

                for ( size_t jhit=0; jhit<infos.size(); ++jhit) // Loop over associated Hits
                  {
                    StrawHitMCInfo const& info = infos.at(jhit);
                    StrawHit const& hit        = info.hit();
                    Straw const& str           = tracker.getStraw(hit.strawIndex());
                    const CLHEP::Hep3Vector mpvec  = str.getMidPoint();
                    const CLHEP::Hep3Vector dirvec = str.getDirection();
                    double dt =hit.dt();
                    double disttomid =dt* _timetodist;
                    CLHEP::Hep3Vector hitpos = mpvec+disttomid*dirvec;
                    CLHEP::Hep3Vector smcpos = CLHEP::Hep3Vector( 0.0, 0.0, 0.0);
                    std::vector<StepPointMC const *> const& steps = info.steps();
                    for ( size_t k=0; k<steps.size(); ++k){
                      StepPointMC const& step = *(steps[k]);
                      smcpos= smcpos+step.position();
                    }
                    smcpos=smcpos/steps.size();
                  }                        // end loop over hits

              }   // end code done if the simparticle is conversionGun
        }         // end if on generatorINdex >= 0
      }           // end loop over simparticles
    if (!foundcele) return;       // no conversion electron found
    Int_t totalHits=0;
    CLHEP::Hep3Vector dvec;

    for ( size_t i=0; i< clusters.size(); ++i ) { // Apparently a loop over clusters
      double hlen=9999999.;
      StrawCluster      const& cluster(clusters.at(i));
      StrawHitPtrVector const& strawHits(cluster.strawHits());

      CLHEP::Hep3Vector pvec = CLHEP::Hep3Vector(0.,0.,0.);
      totalHits = totalHits+strawHits.size();
      Double_t totalEnergy = 0.0;
      for( size_t j=0; j<strawHits.size(); j++ ) {  // Loop over straws in the cluster
        StrawHit const& strawhit = *strawHits[j];
        Double_t Energy = strawhit.energyDep();
        //Double_t Time   = strawhit.time();
        //Double_t deltaT = strawhit.dt();
        totalEnergy=totalEnergy+Energy;
        str = tracker.getStraw(strawhit.strawIndex());
        sid = str.id();
        lid = sid.getLayerId();
        did = sid.getDeviceId();
        secid = sid.getSectorId();
        const CLHEP::Hep3Vector mpvec  = str.getMidPoint();
        const CLHEP::Hep3Vector dirvec = str.getDirection();
        dvec = CLHEP::Hep3Vector(dirvec.getX(),dirvec.getY(),dirvec.getZ());
        //pvec = pvec + Energy * mpvec;     // weight straw by energy deposition
        pvec = pvec + mpvec;
        if (str.getHalfLength()<hlen)
          {
            hlen=str.getHalfLength();
          }
      } // end of loop over straws in the cluster
      double a = 1./double(strawHits.size());
      pvec = pvec*a;
      //pvec = pvec/totalEnergy;             // mean weighted by energy deposition
      pstraw pstr;
      pstr.lay=lid.getLayer();
      assert (did != -1);
      pstr.did=did;
      pstr.sec=secid.getSector();
      pstr.hl=hlen;
      pstr.mpx=pvec.getX();
      pstr.mpy=pvec.getY();
      pstr.mpz=pvec.getZ();
      pstr.dirx=dvec.getX();
      pstr.diry=dvec.getY();
      pstr.dirz=dvec.getZ();
      mpstraws.insert(pair<int,pstraw>(did,pstr));
    }
    //cout << " size of pseudo straw map: " <<mpstraws.size()<<endl;

    // Loop over PAIRS of pseudostraws in the same plane, but only if they intersect
    // This will create doublets, and fill:
    //
    for (int i = 0;i<36;i++)
      {
        if (mpstraws.count(i)>1)
          {
            pair<multimap<int,pstraw>::iterator, multimap<int,pstraw>::iterator> ppp1;
            ppp1 = mpstraws.equal_range(i);
            multimap<int,pstraw>::iterator first11 = ppp1.first;
            multimap<int,pstraw>::iterator first22 = ppp1.first;
            multimap<int,pstraw>::iterator last1 = ppp1.second;
            last1--;
            multimap<int,pstraw>::iterator last2 = ppp1.second;
            for ( multimap<int,pstraw>::iterator first1=first11;first1 != last1;first1++)
              {
                first22=first1;
                first22++;
                for ( multimap<int,pstraw>::iterator first2=first22;first2 != last2;++first2)
                  {
                    pstraw junk  = (*first1).second;
                    pstraw pjunk = (*first2).second;
                    const CLHEP::Hep2Vector p0 =
                      CLHEP::Hep2Vector(junk.mpx-junk.hl*junk.dirx,junk.mpy-junk.hl*junk.diry);
                    const CLHEP::Hep2Vector p1 =
                      CLHEP::Hep2Vector(junk.mpx+junk.hl*junk.dirx,junk.mpy+junk.hl*junk.diry);
                    const CLHEP::Hep2Vector p2 =
                      CLHEP::Hep2Vector(pjunk.mpx-pjunk.hl*pjunk.dirx,pjunk.mpy-pjunk.hl*pjunk.diry);
                    const CLHEP::Hep2Vector p3 =
                      CLHEP::Hep2Vector(pjunk.mpx+pjunk.hl*pjunk.dirx,pjunk.mpy+pjunk.hl*pjunk.diry);
                    LineSegmentPCA linesegment0(p0, p1);
                    LineSegmentPCA linesegment1(p2, p3);
                  } // end for first2
              }// end for first1
          }// end count >1
      }   ///endloop over all devices


//cout << "[[ 10 ]]\n";
    int nclusters=0;
    double _timetodist=149.8962;
    for ( size_t i=0; i< clusters.size(); ++i )// Loop over Clusters
      {
        StrawCluster      const& cluster(clusters.at(i));
        StrawHitPtrVector const& strawHits(cluster.strawHits());

        CLHEP::Hep3Vector pvec = CLHEP::Hep3Vector(0.,0.,0.);
        CLHEP::Hep3Vector clusterpos =  CLHEP::Hep3Vector(0.,0.,0.);
        //totalHits = totalHits+strawHits.size();
        double totalEnergy = 0.0;
        for( size_t j=0; j<strawHits.size(); j++ ) // Loop over Straws in Cluster
          {
            StrawHit const& strawhit = *strawHits[j];
            double Energy = strawhit.energyDep();
            //double Time   = strawhit.time();
            double deltaT = strawhit.dt();
            StrawIndex si   = strawhit.strawIndex();
            Straw str       = tracker.getStraw(si);
            const CLHEP::Hep3Vector mpvec  = str.getMidPoint();
            const CLHEP::Hep3Vector dirvec = str.getDirection();
            double disttomid = deltaT* _timetodist;   // convert delta T into delta x along the wire
            CLHEP::Hep3Vector hitpos = mpvec+disttomid*dirvec;
            clusterpos=clusterpos+hitpos;
            totalEnergy=totalEnergy+Energy;
          } // end loop over straws in cluster
        double a = 1./double(strawHits.size());
        clusterpos=clusterpos*a;
        nclusters++;

      } //  end Loop over Clusters
  } // end of ::analyze.

}

using mu2e::GrokGeometry;
DEFINE_ART_MODULE(GrokGeometry);
