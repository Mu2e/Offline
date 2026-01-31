//
// Original author G. Pezzullo
//

// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Selector.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"

//CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "BTrk/BbrGeom/HepPoint.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include "Offline/RecoDataProducts/inc/TrkFitDirection.hh"
#include "Offline/RecoDataProducts/inc/KalRepPtrCollection.hh"

// From the art tool-chain
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//tracker includes
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BaBar/Constants.hh"
#include "BTrk/TrkBase/HelixParams.hh"
#include "BTrk/TrkBase/TrkRep.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
// conditions
#include "Offline/TrackerGeom/inc/Tracker.hh"
// data
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/MCDataProducts/inc/PtrStepPointMCVector.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/RecoDataProducts/inc/TrkToCaloExtrapol.hh"

//calorimeter includes
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"

// Other includes.
#include "cetlib_except/exception.h"

// Mu2e includes.
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/MCDataProducts/inc/PtrStepPointMCVector.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"

//root includes
#include "TFile.h"
#include "TDirectory.h"
#include "TNtuple.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TGraph.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TMath.h"

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

using namespace std;

namespace mu2e {

  struct IntersectData_t {
      int    fSection;
      int    fRC;                        // return code, 0=success, <0: failure, details TBD
      double fSEntr;
      double fSExit;
    };

  class TrkExtrapol : public art::EDProducer {
  public:

    explicit TrkExtrapol(fhicl::ParameterSet const& pset):
      art::EDProducer{pset},
      _fitterModuleLabel(pset.get<string>("fitterModuleLabel")),
      _tpart((TrkParticle::type)(pset.get<int>("fitparticle",TrkParticle::e_minus))),
      _fdir((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection",TrkFitDirection::downstream))),
      _diagLevel(pset.get<int>("diagLevel",0)),
      _outPutNtup(pset.get<int>("outPutNtup",0)),
      _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel",
                                                  "generate")),
      _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run")),
      _caloReadoutModuleLabel(pset.get<std::string>("caloReadoutModuleLabel",
                                                    "CaloReadoutHitsMaker")),
      _caloCrystalModuleLabel(pset.get<std::string>("caloCrystalModuleLabel",
                                                    "CaloHitsMaker")),
      _directory(0),
      _firstEvent(true),
      _trkdiag(0){
      // Tell the framework what we make.
      produces<TrkToCaloExtrapolCollection>();

      // construct the data product instance name
      _fitDir = _fdir.fitDirection();

    }

    virtual ~TrkExtrapol() {}

    void beginJob();
    void endJob() {}

    void produce(art::Event & e );

    void caloExtrapol(int&             diagLevel,
                      int              evtNumber,
                      TrkFitDirection  fdir,
                      KalRep*          Krep,
                      double&          lowrange,
                      double&          highrange,
                      HelixTraj        &trkHel,
                      int              &res0,
                      int&             NIntersections,
                      IntersectData_t* Intersections);

    double ZfrontFaceCalo() const{ return _ZfrontFaceCalo;}

    double ZbackFaceCalo() const{ return _ZbackFaceCalo;}

  private:

    void doExtrapolation(art::Event & evt, bool skip);
    // Module label of the module that performed the fits.
    std::string _fitterModuleLabel;

    TrkParticle _tpart;

    TrkFitDirection _fdir;

    TrkFitDirection::FitDirection _fitDir;

    // diagnostic of Kalman fit

    // Diagnostic level
    int _diagLevel;
    int _outPutNtup;

    // Label of the generator.
    std::string _generatorModuleLabel;

    // Label of the G4 module
    std::string _g4ModuleLabel;

    // Label of the calo readout hits maker
    std::string _caloReadoutModuleLabel;

    // Label of the calo crystal hists maker
    std::string _caloCrystalModuleLabel;

    bool _skipEvent;

    // Save directory from beginJob so that we can go there in endJob. See note 3.
    TDirectory* _directory;
    bool _firstEvent;

    TTree* _trkdiag;

    Int_t _trkid
      ,_trkint;
    Float_t _trksection[1000]
      , _trkpath[1000]
      , _trktof[1000]
      , _trkx[1000]
      , _trky[1000]
      , _trkz[1000]
      , _trkmomx[1000]
      , _trkmomy[1000]
      , _trkmomz[1000]
      , _trkmom[1000];

    double _solenoidOffSetX;
    double _solenoidOffSetZ;
    double _ZfrontFaceCalo;
    double _ZbackFaceCalo;

    CLHEP::Hep3Vector fromTrkToMu2eFrame(CLHEP::Hep3Vector  &vec);

    void filltrkdiag(int itrk, IntersectData_t *intersec,
                     int size, KalRep const* kalrep);

  };

  CLHEP::Hep3Vector TrkExtrapol::fromTrkToMu2eFrame(CLHEP::Hep3Vector  &vec){
    art::ServiceHandle<GeometryService> geom;
    double solenoidOffSetX = geom->config().getDouble("mu2e.solenoidOffset");
    double solenoidOffSetZ = -geom->config().getDouble("mu2e.detectorSystemZ0");
    CLHEP::Hep3Vector res;

    res.setX(vec.x() - solenoidOffSetX);
    res.setZ(vec.z() - solenoidOffSetZ);
    res.setY(vec.y());
    return res;
  }

  void TrkExtrapol::caloExtrapol(int&             diagLevel,
                                 int              evtNumber,
                                 TrkFitDirection  fdir,
                                 KalRep*          Krep,
                                 double&          lowrange,
                                 double&          highrange,
                                 HelixTraj        &trkHel,
                                 int              &res0,
                                 int&              NIntersections,
                                 IntersectData_t*  Intersection  ) {
    art::ServiceHandle<GeometryService> geom;
    GeomHandle<DiskCalorimeter> cg;
    static const char* oname = "TrkExtrapol::caloExtrapol";

    if(diagLevel>2){

      cout<<"start caloExtrapol, lowrange = "<<lowrange<<
        ", highrange = "<<highrange<<endl;
      cout<<"point of traj at lowrange : "<<Krep->traj().position(lowrange)<<endl;
      cout<<"point of traj at highrange : "<<Krep->traj().position(highrange)<<endl;
      cout<<"fltLMin = "<<Krep->startValidRange()<<
        ", fltLMax = "<<Krep->endValidRange()<<endl;
    }

    TrkErrCode rc;
    rc = Krep->extendThrough(lowrange);

    if (rc.success() != 1 && rc.success() !=13) {
      printf(" %s ERROR: could not extend to lowrange = %10.3f, rc = %i, BAIL OUT\n",oname,lowrange,rc.success());
      return;
    }

    if(diagLevel>2){
      cout<<", after extention..."<<
        ", lowrange = "<<lowrange<<
        ", highrange = "<<highrange<<endl;
      cout<<"point of traj at lowrange : "<<Krep->traj().position(lowrange)<<endl;
      cout<<"point of traj at highrange : "<<Krep->traj().position(highrange)<<endl;
      cout<<"fltLMin = "<<Krep->startValidRange()<<
        ", fltLMax = "<<Krep->endValidRange()<<endl;
    }

    TrkDifTraj const &traj = Krep->traj();

    double circleRadius; // P.Murat: has to be always positive
    double startLowrange = lowrange;
    if(fdir.dzdt() == -1.0) startLowrange = highrange;
    circleRadius = fabs(1.0/trkHel.omega());

    int nSections = cg->nDisks();

    double *entr   = new double[nSections];
    double *ex     = new double[nSections];
    bool *isInside = new bool[nSections];

    for(int jSection=0; jSection<nSections; ++jSection){
      isInside[jSection] = false;
      entr[jSection] = 0.0;
      ex[jSection] = 0.0;
    }
    int nAngleSteps = 500;

    double pathStepSize = Constants::twoPi / (double) nAngleSteps;
    nAngleSteps *= 2.0;

    pathStepSize *= circleRadius/fabs(trkHel.cosDip());

    if (diagLevel>2){
      cout<< "circle radius = "<< circleRadius<<
        ", pathStepSize = "<<pathStepSize<<endl;
    }

    double tmpRange = startLowrange;

    NIntersections = 0;
    CLHEP::Hep3Vector trjVec;
    HepPoint trjPoint;

    for(int iStep = 0; iStep< nAngleSteps; ++iStep){
      for(int jSection=0; jSection<nSections; ++jSection){
        trjPoint = traj.position(tmpRange);
        if(diagLevel>4){
          cout<<" tmpRange = "<< tmpRange<<
            ", trj.position(tmpRange) = "<<trjPoint<<endl;
        }

        trjVec.setX(trjPoint.x());
        trjVec.setY(trjPoint.y());
        trjVec.setZ(trjPoint.z());

        trjVec = fromTrkToMu2eFrame(trjVec);

        if( cg->geomUtil().isInsideSection(jSection,trjVec ) ){
          if(!isInside[jSection]){
            if(diagLevel>4){
              cout<<"Event Number : "<< evtNumber<< endl;
              cout<<" vane "<<jSection<<
                "isInside : true"<<
                "pathLength entrance = "<<tmpRange<<endl;
            }
            isInside[jSection] = true;
            if(fdir.dzdt() == 1.0){
              entr[jSection] = tmpRange - pathStepSize;
            }else if(fdir.dzdt() == -1.0){
              entr[jSection] = tmpRange + pathStepSize;
            }
          }
        }else if(isInside[jSection]){
          ex[jSection] = tmpRange + pathStepSize;
          if(diagLevel>4){
            cout<<"Event Number : "<< evtNumber<< endl;
            cout<<" vane "<<jSection<<
              "isInside : true"<<
              "hasExit : true"<<
              "pathLength entrance = "<<entr[jSection]<<
              "pathLength exit = "<<tmpRange<<endl;
          }
          isInside[jSection] = false;

          if (NIntersections < 100) {
            Intersection[NIntersections].fSection  = jSection;
            Intersection[NIntersections].fRC    = 0;
            Intersection[NIntersections].fSEntr = entr[jSection];
            Intersection[NIntersections].fSExit = ex  [jSection];
            NIntersections++;
          }
          else {
            printf("%s ERROR: NIntersections > 100, TRUNCATE LIST\n",oname);
          }
        }
      }
      if(fdir.dzdt() == 1.0){
        tmpRange += pathStepSize;
      }else if(fdir.dzdt() == -1.0){
        tmpRange -= pathStepSize;
      }
    }
    //  }

    if (diagLevel>2) {
      cout<<"end search behindSection(), position is : "<<traj.position(tmpRange)<<endl;
    }

    double     lrange;
    TrkErrCode trk_rc;

    for (int i=0; i<NIntersections; i++) {
      lrange = Intersection[i].fSEntr;
      trk_rc = Krep->extendThrough(lrange);
      if (trk_rc.success() != 1) {
        //-----------------------------------------------------------------------------
        // failed to extend
        //-----------------------------------------------------------------------------
        Intersection[i].fRC = -1;
        if (diagLevel>2) {
          printf("%s ERROR vane = %2i FAILED to EXTEND TRAJECTORY, rc = %i\n",
                 oname,Intersection[i].fSection,trk_rc.success());
        }
      }

    }

    delete [] isInside ;
    delete [] entr ;
    delete [] ex;

  }//end proce_dUre

  void TrkExtrapol::beginJob() {

    if (_outPutNtup == 1) {
      art::ServiceHandle<art::TFileService> tfs;
      _trkdiag       = tfs->make<TTree>("trk", "trk extrapolated info");
      _trkdiag->Branch("trkid", &_trkid  ,"trkid/I");
      _trkdiag->Branch("trkint", &_trkint  ,"trkint/I");
      _trkdiag->Branch("trksection[trkint]", _trksection, "trksection[trkint]/F");
      _trkdiag->Branch("trkpath[trkint]", _trkpath, "trkpath[trkint]/F");
      _trkdiag->Branch("trktof[trkint]", _trktof, "trktof[trkint]/F");
      _trkdiag->Branch("trkx[trkint]", _trkx, "trkx[trkint]/F");
      _trkdiag->Branch("trky[trkint]", _trky, "trky[trkint]/F");
      _trkdiag->Branch("trkz[trkint]", _trkz, "trkz[trkint]/F");
      _trkdiag->Branch("trkmomx[trkint]", _trkmomx, "trkmomx[trkint]/F");
      _trkdiag->Branch("trkmomy[trkint]", _trkmomy, "trkmomy[trkint]/F");
      _trkdiag->Branch("trkmomz[trkint]", _trkmomz, "trkmomz[trkint]/F");
      _trkdiag->Branch("trkmom[trkint]", _trkmom, "trkmom[trkint]/F");
    }

  }

  void TrkExtrapol::filltrkdiag(int itrk, IntersectData_t *intersec, int size, KalRep const* kalrep){
    _trkid = itrk;
    double lenght(0.0);
    _trkint = size;
    TrkDifTraj const &traj = kalrep->traj();
    for(int i=0; i<size; ++i){
      _trksection[i] = intersec[i].fSection;
      lenght = intersec[i].fSEntr;
      _trkpath[i] = lenght;
      _trktof[i] =  kalrep->arrivalTime(lenght);
      _trkx[i] = traj.position(lenght).x();
      _trky[i] = traj.position(lenght).y();
      _trkz[i] = traj.position(lenght).z();
      _trkmomx[i] = kalrep->momentum(lenght).x();
      _trkmomy[i] = kalrep->momentum(lenght).y();
      _trkmomz[i] = kalrep->momentum(lenght).z();
      _trkmom[i]  = kalrep->momentum(lenght).mag();

    }
    _trkdiag->Fill();
  }

//-----------------------------------------------------------------------------
  void TrkExtrapol::produce(art::Event & evt ) {

    doExtrapolation(evt, _skipEvent);
  }

//-----------------------------------------------------------------------------
  void TrkExtrapol::doExtrapolation(art::Event & evt, bool skip){

    art::ServiceHandle<GeometryService> geom;
    GeomHandle<Calorimeter> cg;
    _solenoidOffSetX = geom->config().getDouble("mu2e.solenoidOffset");//3904.;//[mm]
    _solenoidOffSetZ = -geom->config().getDouble("mu2e.detectorSystemZ0");//-10200.;

    _ZfrontFaceCalo = cg->geomUtil().origin().z() + _solenoidOffSetZ;
    _ZbackFaceCalo = cg->geomUtil().origin().z() + _solenoidOffSetZ;

    const char* oname = "TrkExtrapol::doExtrapolation";
    double      lowrange, highrange, zmin, zmax;
    HepPoint    point;
    int         ntrk, res0;

    //create output
    unique_ptr<TrkToCaloExtrapolCollection> extrapolatedTracks(new TrkToCaloExtrapolCollection );
    TrkToCaloExtrapolCollection tmpExtrapolatedTracks;

    art::Handle<KalRepPtrCollection> trksHandle;
    evt.getByLabel(_fitterModuleLabel,trksHandle);
    const KalRepPtrCollection* trks = trksHandle.product();

    ntrk = trks->size();

    if(_diagLevel>2){
      cout<<endl<<"Event Number : "<< evt.event()<< endl;
      cout<<"\n start TrkExtrapol..."<<endl;
      cout<<"ntrk = "<< ntrk <<endl;
    }

    double circleRadius = 0.0, centerCircleX=0.0, centerCircleY = 0.0, angle = 0.0;

    for (int itrk=0; itrk< ntrk; ++itrk ){
      res0 = -1;
                                        // extrapolation extends the track and thus changes it...

      KalRep* krep = (KalRep*) trks->at(itrk).get();
      if ( !krep ) continue;
      TrkDifTraj const& traj = krep->traj();
      double pos = 0.0;

      double endTrk = krep->endFoundRange();
      // starting from the end of the tracker!!!FIXME
      HelixTraj trkHel(krep->helix(endTrk).params(),krep->helix(endTrk).covariance());

      angle = Constants::pi*0.5 + trkHel.phi0();

      circleRadius  = 1.0/trkHel.omega();
      centerCircleX = trkHel.d0() + circleRadius;
      centerCircleY = centerCircleX*sin(angle);
      centerCircleX *= cos(angle);

      // 2013-05-23 gianipez : add 10cm tolerance
      zmin = _ZfrontFaceCalo - 100.;
      zmax = _ZbackFaceCalo  + 100.;

      if(_fitDir ==  TrkFitDirection::downstream){
        lowrange  = trkHel.zFlight(zmin);  /*1740*/
        highrange = trkHel.zFlight(zmax); /*3500*/
      }else if(_fitDir ==  TrkFitDirection::upstream ){
        lowrange  = trkHel.zFlight(zmax);
        highrange = trkHel.zFlight(zmin);
      }

      if (_diagLevel>2) {

        cout<<endl<<"Event Number : "<< evt.event()<< endl;
        cout<<"------ trk number : "<<itrk<<" ------"<<endl;
        cout<<"found traj, point of traj at "<<traj.position(pos)<<endl;
        cout<<"*************** lowRange =  "<< lowrange <<", highRange = "<< highrange << endl;
        cout<< " is the particle in the 0 quadrant?"<<endl<<
          ", traj.position(lowrange).x() = "<< traj.position(lowrange).x()<<
          ", traj.position(lowrange).y() = "<< traj.position(lowrange).y()<<
          ", traj.position(lowrange).z() = "<< traj.position(lowrange).z()<<endl;

        printf("circle: R = %10.3f X0 = %10.3f  Y0 = %10.3f phi0 = %10.3f\n",
               circleRadius,centerCircleX,centerCircleY,angle);
      }

      IntersectData_t  intersection[100];
      int                                     nint(0);

      caloExtrapol(_diagLevel,
                   (int) evt.event(),
                   _fitDir, krep, lowrange, highrange,
                   trkHel,
                   res0,
                   nint,
                   intersection);

      if (nint == 0) {
        printf("\n%s , run / event : %d / %d, \nERROR: intersection not found : res0 = %i\nfitdirection = %s \n",
               oname,
                evt.id().run(), evt.id().event(),
               res0,
               _fdir.name().c_str());
        point = krep->traj().position(lowrange);
        printf("point of trj at lowrange(%10.3f)  : ( %10.3f, %10.3f, %10.3f )\n",
               lowrange,
               point.x(), point.y(), point.z());

        point = krep->traj().position(highrange);
        printf("point of trj at highrange(%10.3f) : ( %10.3f, %10.3f, %10.3f )\n",
               highrange,
               point.x(), point.y(), point.z());
      }

      if(_outPutNtup ==1){
        filltrkdiag(int(itrk), intersection, nint, krep);
      }

      for (int i=0; i<nint; i++) {
        KalRepPtr tmpRecTrk = trksHandle->at(itrk);
        tmpExtrapolatedTracks.push_back(
                                        TrkToCaloExtrapol(intersection[i].fSection,
                                                          itrk,
                                                          tmpRecTrk,
                                                          intersection[i].fSEntr,
                                                          intersection[i].fSExit)
                                        );
      }

      // P.Murat: why would one need to sort at this point?

      //      std::sort(tmpExtrapolatedTracks.begin(), tmpExtrapolatedTracks.end());
      for(TrkToCaloExtrapolCollection::iterator it = tmpExtrapolatedTracks.begin(); it != tmpExtrapolatedTracks.end(); ++it){
        extrapolatedTracks->push_back(*it);
      }

    }//end loop on recoTrj

    evt.put(std::move(extrapolatedTracks));
    if( evt.id().event() %100 ==0){
      printf("\nEvent %d %s done...\n",evt.id().event(),oname );
    }
  }

}

using mu2e::TrkExtrapol;
DEFINE_ART_MODULE(TrkExtrapol)
