///////////////////////////////////////////////////////////////////////////////
//
// Original author G. Pezzullo
//
// 2014-06-03 P.Murat: will no longer work with the vanes-based geometry
// 2015-09-17 P.Murat: use Bertran's extrapolator
///////////////////////////////////////////////////////////////////////////////

// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Selector.h"
#include "art_root_io/TFileDirectory.h"

// From the art tool-chain
#include "fhiclcpp/ParameterSet.h"

// conditions
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"

#include "Offline/RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "Offline/RecoDataProducts/inc/TrkFitDirection.hh"

#include "Offline/RecoDataProducts/inc/TrkCaloIntersect.hh"
#include "Offline/RecoDataProducts/inc/TrackClusterMatch.hh"

#include "Offline/RecoDataProducts/inc/CaloCluster.hh"

#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"

#include "BTrk/TrkBase/HelixParams.hh"

// Other includes.
#include "cetlib_except/exception.h"

//root includes
#include "TFile.h"
#include "TDirectory.h"
#include "TMath.h"
#include "TVector2.h"

// From the art tool-chain
#include <cmath>
#include <deque>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <functional>
#include "cetlib/pow.h"
#include "BTrk/BaBar/Constants.hh"


using namespace std;
using cet::square;
using cet::sum_of_squares;
using CLHEP::Hep3Vector;

namespace mu2e {

  class TrackCaloMatching : public art::EDProducer {
  private:

    art::ProductToken<KalRepPtrCollection>        const _fitterToken;
    art::ProductToken<CaloClusterCollection>      const _caloClusterToken;
    art::ProductToken<TrkCaloIntersectCollection> const _trkToCaloExtrapolToken;

    int             _debugLevel;

    double          _minClusterEnergy;      //
    double          _maxDeltaT;             // time preselection for track-cluster matching
    double          _meanInteractionDepth;  // path length correction
    double          _dtOffset;              // shift of the Delta(T) distribution

                                            // resolutions used to calculate chi2's
    double          _sigmaE;
    double          _sigmaT;
    double          _sigmaU;
    double          _sigmaV;

    double          _chi2e, _chi2t, _chi2u, _chi2v;

                                            // no offset in Y ?
    double          _solenoidOffSetX;
    double          _solenoidOffSetZ;

    bool            _skipEvent;
    bool            _firstEvent;

  public:

    explicit TrackCaloMatching(fhicl::ParameterSet const& pset):
      art::EDProducer{pset},
      _fitterToken           {consumes<KalRepPtrCollection>       (pset.get<string>     ("fitterModuleLabel"           ))},
      _caloClusterToken      {consumes<CaloClusterCollection>     (pset.get<std::string>("caloClusterModuleLabel"      ))},
      _trkToCaloExtrapolToken{consumes<TrkCaloIntersectCollection>(pset.get<std::string>("trkToCaloExtrapolModuleLabel"))},

      _debugLevel            (pset.get<int>   ("debugLevel"          )),
      _minClusterEnergy      (pset.get<double>("minClusterEnergy"    )),
      _maxDeltaT             (pset.get<double>("maxDeltaT"           )),
      _meanInteractionDepth  (pset.get<double>("meanInteractionDepth")),
      _dtOffset              (pset.get<double>("dtOffset"            )),
      _sigmaE                (pset.get<double>("sigmaE"              )),
      _sigmaT                (pset.get<double>("sigmaT"              )),
      _sigmaU                (pset.get<double>("sigmaU"              )),
      _sigmaV                (pset.get<double>("sigmaV"              )),
      _firstEvent            (true)
    {
//-----------------------------------------------------------------------------
// Tell the framework what we make.
//-----------------------------------------------------------------------------
      produces<TrackClusterMatchCollection>();
    }

    void beginJob() override;
    void beginRun(art::Run& aRun) override;
    void produce (art::Event& e) override;

    void fromTrkToMu2eFrame(HepPoint& vec, CLHEP::Hep3Vector& res);

    void doMatching(art::Event & evt, bool skip);
  };


  //-----------------------------------------------------------------------------
  void TrackCaloMatching::fromTrkToMu2eFrame(HepPoint& vec, Hep3Vector& res) {
    res.setX(vec.x() - _solenoidOffSetX);
    res.setZ(vec.z() - _solenoidOffSetZ);
    res.setY(vec.y());
  }

  //-----------------------------------------------------------------------------
  void TrackCaloMatching::beginJob() {

    if (_debugLevel > 0) {
      printf("---- TrackCaloMatching::beginJob constants used: \n"    );
      printf("  minClusterEnergy     : %10.3f\n",_minClusterEnergy    );
      printf("  maxDeltaT            : %10.3f\n",_maxDeltaT           );
      printf("  meanInteractionDepth : %10.3f\n",_meanInteractionDepth);
      printf("  dtOffset             : %10.3f\n",_dtOffset            );
    }
  }

  //-----------------------------------------------------------------------------
  void TrackCaloMatching::beginRun(art::Run& aRun) {
    art::ServiceHandle<GeometryService> geom;

    // calculate DS offsets

    _solenoidOffSetX =  geom->config().getDouble("mu2e.solenoidOffset");
    _solenoidOffSetZ = -geom->config().getDouble("mu2e.detectorSystemZ0");
  }

  //-----------------------------------------------------------------------------
  void TrackCaloMatching::produce(art::Event & evt) {
    doMatching(evt, _skipEvent);
  }

  //-----------------------------------------------------------------------------
  void TrackCaloMatching::doMatching(art::Event & evt, bool skip) {
    constexpr const char* oname = "TrackCaloMatching::doMatching";

    double     cl_x, cl_y, cl_time, cl_energy;
    double     trk_x, trk_y, trk_mom, trk_time;
    double     s1, s2, sint, ds;
    double     nx, ny, dx, dy, du, dv, dt, xu, xv, xt, xe;
    double     chi2;
    constexpr double chi2_max(1.e12);

    int        idisk, iex, icl, ltrk;
    constexpr int  ndisks(2);

    TrackClusterMatch::Data_t tcm_data[100][ndisks];

    const TrkCaloIntersect     *extrk;
    const KalRep               *krep;
    const CaloCluster          *cl;
    // cp_mu2e: cluster position in the MU2E frame

    CLHEP::Hep3Vector           tp_disk, tp_mu2e, cp_mu2e, cp_st;
    CLHEP::Hep3Vector           mom, pos;
    HepPoint                    point, p1, p2, p12, p_closest;

    const KalRep               *kkrep[100];
    int                        nkkrep(0);
    //-----------------------------------------------------------------------------
    // Get handle to calorimeter
    //-----------------------------------------------------------------------------
    art::ServiceHandle<GeometryService> geom;
    GeomHandle<Calorimeter> cg;

    auto const& trjExtrapols        = evt.getValidHandle(_trkToCaloExtrapolToken);
    int const   nex                 = trjExtrapols->size();

    auto const& trksHandle          = evt.getValidHandle(_fitterToken);
    const KalRepPtrCollection* trks = trksHandle.product();
    int  const  ntracks             = trks->size();

    auto const& caloClusters        = evt.getValidHandle(_caloClusterToken);
    int  const  nclusters           = caloClusters->size();

    if (_debugLevel > 2) {
      printf("%s: Event Number : %8i ntracks: %2i nclusters: %4i nex: %2i\n",
             oname,evt.event(), ntracks, nclusters, nex);
    }

    unique_ptr<TrackClusterMatchCollection> tcmcoll(new TrackClusterMatchCollection);
    tcmcoll->reserve(nex);
    //-----------------------------------------------------------------------------
    // here the execution really starts
    //-----------------------------------------------------------------------------
    if (ntracks == 0)                                         goto END;

    for (int it=0; it<ntracks; it++) {
      for (int iv=0; iv<ndisks; iv++) {
        tcm_data[it][iv].chi2 = chi2_max+1;
        tcm_data[it][iv].iex  = -1;
        tcm_data[it][iv].icl  = -1;
      }
    }
    //-----------------------------------------------------------------------------
    // loop over the extrapolated tracks (TrkToCaloExtrapol's),
    // each of them corresponds to a track-to-disk intersection, so one track may
    // have two corresponding TrkToCaloExtrapol things
    //-----------------------------------------------------------------------------
    for (int jex=0; jex<nex; jex++) {

      if(_debugLevel > 2) {
        printf("%s: jex = %5i\n", oname, jex);
      }

      extrk = &trjExtrapols->at(jex);
      krep  = extrk->trk().get();
      //-----------------------------------------------------------------------------
      // track index, important: store one, the best, intersection per track per vane
      // the absolute ltrk number doesn't really matter
      //-----------------------------------------------------------------------------
      ltrk = nkkrep;
      for (int ik=0; ik<nkkrep; ik++) {
        if (krep == kkrep[ik]) {
          ltrk = ik;
        }
      }

      if (ltrk == nkkrep) {
        kkrep[nkkrep] = krep;
        nkkrep       += 1;
      }

      if (ltrk > 100) {
        printf(">>> ERROR in %s: ltrk = %i, skip the rest\n",oname,ltrk);
        goto NEXT_INTERSECTION;
      }
      else if (ltrk < 0) {
        printf(">>> ERROR in %s: ltrk = %i, skip the rest\n",oname,ltrk);
        goto NEXT_INTERSECTION;
      }
      //-----------------------------------------------------------------------------
      // apparently, ntupling stuff was mixed in, almost removed
      //-----------------------------------------------------------------------------
      idisk  = extrk->diskId();
      // assume Z(middle of the disk)

      s1       = extrk->pathLengthEntrance();
      s2       = extrk->pathLengthExit    ();
      ds       = s2-s1;
      // 'sint' - extrapolation length to the interaction point

      if (ds > _meanInteractionDepth) sint = s1+_meanInteractionDepth;
      else                            sint = s2;

      // shower starts developing when the particle
      // reaches the disk
      trk_time = krep->arrivalTime(sint);
      mom      = krep->momentum(sint);

      point    = krep->position(sint);
      mom      = krep->momentum(sint);
      trk_mom  = mom.mag();

      if (_debugLevel > 2) {
        printf("%s: idisk: %5i trk_mom: %10.5f [MeV/c]\n",oname,idisk,trk_mom);
      }
      // track direction cosines in XY plane

      nx      = mom.x()/sqrt(mom.x()*mom.x()+mom.y()*mom.y());
      ny      = mom.y()/sqrt(mom.x()*mom.x()+mom.y()*mom.y());

      // transform track position to Mu2e detector frame

      fromTrkToMu2eFrame(point, tp_mu2e);

      if (_debugLevel > 2){
      }

      tp_disk = cg->geomUtil().mu2eToDisk(idisk,tp_mu2e);
      trk_x   = tp_disk.x();
      trk_y   = tp_disk.y();

      if(_debugLevel > 2){
        printf("%s: tp_mu2e: %10.4f %10.4f %10.4f  tp_disk: %10.4f %10.4f %10.4f\n",
               oname,
               tp_mu2e.x(),tp_mu2e.y(),tp_mu2e.z(),
               tp_disk.x(),tp_disk.y(),tp_disk.z());
      }
      //-----------------------------------------------------------------------------
      // save track-only information
      //-----------------------------------------------------------------------------
      tcm_data[ltrk][idisk].xtrk      = tp_disk.x();
      tcm_data[ltrk][idisk].ytrk      = tp_disk.y();
      tcm_data[ltrk][idisk].ztrk      = tp_disk.z();
      tcm_data[ltrk][idisk].ttrk      = trk_time;
      tcm_data[ltrk][idisk].nx        = mom.x()/trk_mom;
      tcm_data[ltrk][idisk].ny        = mom.y()/trk_mom;
      tcm_data[ltrk][idisk].nz        = mom.z()/trk_mom;
      tcm_data[ltrk][idisk].int_depth = _meanInteractionDepth;
      tcm_data[ltrk][idisk].ds        = ds;

      //-----------------------------------------------------------------------------
      // track helix at Z = Z(middle of the disk) in the tracker frame
      //-----------------------------------------------------------------------------
      double     trk_d0, trk_om, trk_r, trk_phi0, trk_phi1, trk_x0, trk_y0, trk_tandip;
      double     cp_dx, cp_dy, cp_phi, cp_dphi, delta_x, delta_y, s12, s_cl;
      double     dds, dz, dr, sint;


      p1       = krep->position(s1);
      p2       = krep->position(s2);

      s12      = (s1+s2)/2;

      p12       = krep->position(s12);

      trk_d0     = krep->helix(s12).d0();
      trk_om     = krep->helix(s12).omega();
      trk_r      = fabs(1./trk_om);
      trk_phi0   = krep->helix(s12).phi0();
      trk_tandip = krep->helix(s12).tanDip();
      trk_x0     =  -(1/trk_om+trk_d0)*sin(trk_phi0);
      trk_y0     =   (1/trk_om+trk_d0)*cos(trk_phi0);
      trk_phi1   = atan2(p12.y()-trk_y0,p12.x()-trk_x0);
      //-----------------------------------------------------------------------------
      // loop over clusters
      //-----------------------------------------------------------------------------
      for (int icl=0; icl<nclusters; icl++) {
        cl      = &(*caloClusters).at(icl);
        cl_time = cl->time();
        // move peak to zero
        dt      = trk_time-cl_time-_dtOffset;

        if (cl->diskID() != idisk           )            goto NEXT_CLUSTER;
        if (cl->energyDep() < _minClusterEnergy)            goto NEXT_CLUSTER;
        if (std::fabs(dt)   > _maxDeltaT       )            goto NEXT_CLUSTER;
        //------------------------------------------------------------------------------
        // 2015-03-26: as of now, the cluster coordinates are defined in the local disk
        //             coordinate system
        //-----------------------------------------------------------------------------
        cl_x         = cl->cog3Vector().x();
        cl_y         = cl->cog3Vector().y();
        cl_energy    = cl->energyDep();
        //-----------------------------------------------------------------------------
        // X and Y - coordinates in the disk system
        // rotate them in the direction perpendicular to the track
        //-----------------------------------------------------------------------------
        dx  = trk_x-cl_x;
        dy  = trk_y-cl_y;

        du = dx*nx+dy*ny;
        dv = dx*ny-dy*nx;
        //-----------------------------------------------------------------------------
        // to study pattern recognition accuracies, calculate track parameters in the
        // point closest to the cluster in 2D (XY) - want to understand coordinate and
        // angle resolutions
        //-----------------------------------------------------------------------------
        cp_mu2e = cg->geomUtil().diskToMu2e(idisk,cl->cog3Vector());
        cp_st   = cg->geomUtil().mu2eToTracker(cp_mu2e);

        cp_dx   = cp_st.x()-trk_x0;
        cp_dy   = cp_st.y()-trk_y0;
        // cluster phi wrt helix center
        cp_phi  = atan2(cp_dy,cp_dx);
        cp_dphi = TVector2::Phi_mpi_pi(cp_phi-trk_phi1);

        // radial distance

        dr      = sqrt(cp_dx*cp_dx+cp_dy*cp_dy)-trk_r;
        delta_x = dr*cos(cp_dphi);
        delta_y = dr*sin(cp_dphi);
        // last: track Z in the closest point

        dds     = cp_dphi*trk_r*sqrt(1+trk_tandip*trk_tandip);
        if (trk_om < 0) dds = -dds;

        s_cl    = s12+dds;

        p_closest = krep->position(s_cl);
        dz        = p_closest.z()-cp_st.z();
        sint      = s_cl-s1;

        if (_debugLevel > 2) {
          printf("%s: s1     : %9.3f p1     : %9.3f %9.3f %9.3f\n",oname,s1, p1.x(),p1.y(),p1.z());
          printf("%s: s2     : %9.3f p2     : %9.3f %9.3f %9.3f\n",oname,s2, p2.x(),p2.y(),p2.z());
          printf("%s: s12    : %9.3f p12    : %9.3f %9.3f %9.3f\n",oname,s12,p12.x(),p12.y(),p12.z());
          printf("%s: mom    : %9.3f %9.3f %9.3f tandip: %9.3f\n",oname,mom.x(),mom.y(),mom.z(),trk_tandip);
          printf("%s: dds    : %9.3f\n",oname,dds);
          printf("%s: s_cl   : %9.3f\n",oname,s_cl);
          printf("%s: p_cl   : %9.3f %9.3f %9.3f\n",oname,p_closest.x(),p_closest.y(),p_closest.z());

          printf("%s: ",oname);
          printf("trk_x0 : %9.3f ",trk_x0);
          printf("trk_y0 : %9.3f ",trk_y0);
          printf("trk_ph0: %9.3f ",trk_phi0);
          printf("trk_ph1: %9.3f ",trk_phi1);
          printf("trk_r  : %9.3f ",trk_r );
          printf("\n");

          printf("%s: ",oname);
          printf("cl_x   : %9.3f ",cl_x  );
          printf("cl_y   : %9.3f ",cl_y  );
          printf("cp_phi : %9.3f ",cp_phi);
          printf("cp_dphi: %9.3f ",cp_dphi);
          printf("\n");

          printf("%s: ",oname);
          printf("delta_x: %9.3f ",delta_x);
          printf("delta_y: %9.3f ",delta_y);
          printf("dr     : %9.3f ",dr    );
          printf("dz     : %9.3f ",dz    );
          printf("sint   : %9.3f ",sint  );
          printf("\n");
        }
        //-----------------------------------------------------------------------------
        // 2014-05-14 P.Murat: use 10 MeV as the matching resolution,
        //            for a 10 MeV cluster the energy term in the chi2 would be (90/10)^2
        //            also set coordinate resolution to 5cm , need to try using dx only
        //-----------------------------------------------------------------------------
        xu       = du/_sigmaU;
        xv       = dv/_sigmaV;
        xt       = dt /_sigmaT;
        xe       = (trk_mom-cl_energy)/_sigmaE;

        _chi2u   = xu*xu;
        _chi2v   = xv*xv;
        _chi2t   = xt*xt;
        _chi2e   = xe*xe;

        chi2 = _chi2u + _chi2v;

        if (_debugLevel > 2){
          printf("%s: ",oname);
          printf("trk_x  : %9.3f "  ,trk_x  );
          printf("cl_x   : %9.3f "  ,cl_x   );
          printf("trk_y  : %9.3f "  ,trk_y  );
          printf("cl_y   : %9.3f "  ,cl_y   );
          printf("cl_time: %9.3f "  ,cl_time);
          printf("\n");
          printf("%s: ",oname);
          printf("du     : %9.3f "  ,du     );
          printf("dv     : %9.3f "  ,dv     );
          printf("xu     : %9.3f "  ,xu     );
          printf("xv     : %9.3f "  ,xv     );
          printf("chi2   : %9.3f "  ,chi2   );
          printf("\n");
        }

        if (chi2 < tcm_data[ltrk][idisk].chi2) {
          //-----------------------------------------------------------------------------
          // new best match
          //-----------------------------------------------------------------------------
          tcm_data[ltrk][idisk].icl       = icl;
          tcm_data[ltrk][idisk].iex       = jex;
          tcm_data[ltrk][idisk].dx        = dx;
          tcm_data[ltrk][idisk].dy        = dy;
          tcm_data[ltrk][idisk].dz        = -1e6;
          tcm_data[ltrk][idisk].du        = du;
          tcm_data[ltrk][idisk].dv        = dv;
          tcm_data[ltrk][idisk].dt        = dt;
          tcm_data[ltrk][idisk].ep        = cl_energy/trk_mom;
          tcm_data[ltrk][idisk].chi2      = chi2;
          tcm_data[ltrk][idisk].chi2_time = _chi2t;
          tcm_data[ltrk][idisk].dr        = dr;
          tcm_data[ltrk][idisk].sint      = sint;
        }
      NEXT_CLUSTER:;
      }
    NEXT_INTERSECTION:;
    }

    //-----------------------------------------------------------------------------
    // form output list of matches
    //-----------------------------------------------------------------------------
    for (int it=0; it<ntracks; it++) {
      for (int iv=0; iv<ndisks; iv++) {
        if (tcm_data[it][iv].chi2 < chi2_max) {
          iex = tcm_data[it][iv].iex;
          icl = tcm_data[it][iv].icl;

          art::Ptr<TrkCaloIntersect>  artPtrTex    (trjExtrapols,iex);
          art::Ptr<CaloCluster>       artPtrCluster(caloClusters,icl);
          TrackClusterMatch           tcm(artPtrTex,artPtrCluster,&tcm_data[it][iv]);
          tcmcoll->push_back(tcm);
        }
      }
    }

  END:;

    evt.put(std::move(tcmcoll));
    //-----------------------------------------------------------------------------
    // diagnostics printout
    //-----------------------------------------------------------------------------
    if (_debugLevel > 0) {
      if (evt.id().event() %100 == 0) printf("Event %d TrackCaloMatching done.\n",evt.id().event());
    }
  }

}

DEFINE_ART_MODULE(mu2e::TrackCaloMatching)
