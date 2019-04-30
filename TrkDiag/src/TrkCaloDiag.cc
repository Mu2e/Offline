//
// Diagnostics for track-calo matching
// $Author: brownd $ 
// $Date: 2014/09/20 14:34:22 $
//

#include "TrkDiag/inc/TrkCaloDiag.hh"
namespace mu2e 
{  

  TrkCaloDiag::TrkCaloDiag( TrkParticle const& tpart, TrkFitDirection const& fdir, fhicl::ParameterSet const& pset) : 
    _ncalo(0),
    _pid_dt(pset.get<fhicl::ParameterSet>("PIDdt",fhicl::ParameterSet())),
    _pid_ep(pset.get<fhicl::ParameterSet>("PIDEp",fhicl::ParameterSet())) 
  {
// construct the calo matching module name.  Convention is 1st letter of direction, 1st letter of particle name + 1st letter of charge
// This code will break if the fhicl prolog conventions change, FIXME!!!
    std::string chargename = tpart.charge() > 0.0 ? "p" : "m";
    std::string caloMatchingRoot = pset.get<std::string>("caloMatchingRoot","TrackCaloMatching");
    _caloMatchModule = caloMatchingRoot + fdir.name().substr(0,1) + tpart.name().substr(0,1) + chargename;
  }

  void TrkCaloDiag::addCaloInfo(KalRep const* krep) {
    if(_caloMatchHandle.isValid()){
      _caloinfo.clear();
      _ncalo = 0;
      for( auto const& tcm : *_caloMatchHandle ) {
	if(tcm.textrapol()->trk().get() == krep){
	  TrkCaloInfo tcinfo;
	  fillCaloInfo(tcm,tcinfo);
	  _caloinfo.push_back(tcinfo);
	  _ncalo++;
	}
      }
    }
  }

  void TrkCaloDiag::fillCaloInfo(TrackClusterMatch const& tcm, TrkCaloInfo& tcinfo) {
// matching info
    tcinfo._dt = tcm.dt();
    tcinfo._du = tcm.du();
    tcinfo._dv = tcm.dv();
    tcinfo._ds = tcm.ds();
    tcinfo._ep = tcm.ep();
    tcinfo._uvChisq = tcm.chi2();
    tcinfo._tChisq = tcm.chi2_time();
// PID information
    tcinfo._dtllr = _pid_dt.value(tcm.dt());
    tcinfo._epllr = _pid_ep.value(tcm.ep(),tcm.ds());
// cluster info
    const CaloCluster* cluster = tcm.caloCluster();
    tcinfo._eclust = cluster->energyDep();
    tcinfo._tclust = cluster->time();
    tcinfo._section = cluster->diskId();
    tcinfo._cpos = Geom::toXYZVec(cluster->cog3Vector());
// track information at intersection point.  Don't use this as there's an
// additional fltlen added for the depth (59mm).
//  KalRep const* krep = tcm.textrapol()->trk()  
//    double ipath = tcinfo.textrapol()->pathLengthEntrance();	       
//    tcinfo._tpos = Geom::toXYZVec(krep.position(ipath); 
//    tcinfo._tdir = Geom::toXYZVec(krep.direction(ipath);
//    tcinfo._ttrk = krel.arrivalTime(ipath);
    tcinfo._tpos = Geom::toXYZVec(CLHEP::Hep3Vector(tcm.xtrk(),tcm.ytrk(),tcm.ztrk()));
    tcinfo._tdir = Geom::toXYZVec(CLHEP::Hep3Vector(tcm.nx(),tcm.ny(),tcm.nz()));
    tcinfo._ttrk = tcm.ttrk();
  }

  void TrkCaloDiag::addBranches(TTree* tree,const char* suffix) {
    std::string ssuf(suffix);
    std::string bname = std::string("ncalo") +ssuf;
    tree->Branch(bname.c_str(),&_ncalo,"ncalo/I");
    bname = std::string("calo") +ssuf;
    tree->Branch(bname.c_str(),&_caloinfo);
  }

  void TrkCaloDiag::findData(const art::Event& event) {
    event.getByLabel(_caloMatchModule,_caloMatchHandle);
    // clear the branch in case there are no calo objects found
    _ncalo = 0;
    _caloinfo.clear();
  }

}

