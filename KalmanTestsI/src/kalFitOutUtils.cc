//
// output utilities for reco modules
//
// $Id: kalFitOutUtils.cc,v 1.2 2012/12/04 00:51:27 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/12/04 00:51:27 $
//

#include "KalmanTestsI/inc/kalFitOutUtils.hh"
#include "BaBar/Constants.hh"
#include "KalmanTests/inc/TrkDef.hh"
#include "TrkBase/TrkMomCalculator.hh"
#include "ConditionsService/inc/TrackerCalibrationsI.hh"

using namespace std; 

namespace mu2e 
{

void kalFitOutUtils::bookHitos()
  {
    art::ServiceHandle<art::TFileService> tfs;
    
    hpull  =tfs->make<TH2D>("hpull","pull;param;#delta parameter/#sigma",5,0,5,10000,-500,500);
    hdev   =tfs->make<TH2D>("hdev","dev",5,0,5,200000,-1,1);
    hchi2  =tfs->make<TH1D>("hchi2","chi2/ndf;#chi^{2}/ndf",100000,0,1000);
    hprob  =tfs->make<TH1D>("hprob","prob",10000,0,1);
    hchi2_0=tfs->make<TH1D>("hchi2_0","chi2",100000,0,1000);
    hmom   =tfs->make<TH1D>("hmom",";#delta P, MeV/c",10000,-100,100);

    hpars[0]=tfs->make<TH1D>("d0",";d0, mm",400,-30,30);
    hpars[1]=tfs->make<TH1D>("phi0",";#phi0, rad",400,-0.05,0.05);
    hpars[2]=tfs->make<TH1D>("omega",";#omega, 1/mm",400,-2e-4,2e-4);
    hpars[3]=tfs->make<TH1D>("z0",";z0, mm",400,-30,30);
    hpars[4]=tfs->make<TH1D>("tan",";tan",400,-0.05,0.05);

    hnhits     =tfs->make<TH1D>("hnhits",";number of hits per turn",500,0,500);
    hnhitstotal=tfs->make<TH1D>("hnhitstotal",";total number of hits",500,0,500);
    hnturns    =tfs->make<TH1D>("hnturns",";number of turns",10,0,10);
    heloss     =tfs->make<TH1D>("heloss",";E loss from first to last cell, MeV",10000,0,10);
    helossIWall=tfs->make<TH1D>("helossIWall",";E loss on inner wall, MeV",10000,0,10);
    helossIWallNorm=tfs->make<TH1D>("helossIWallNorm",";E loss on inner wall normilized to perp, MeV",10000,0,10);

    hdedx      =tfs->make<TH1D>("hdedx",";dEdX , MeV/cm",500000,0,10);
    hdedxturn  =tfs->make<TH1D>("hdedxturn",";dEdX for full turn, MeV/cm",500000,0,10);

    hdtheta     = tfs->make<TH1D>("hdtheta",";dtheta, rad/#sqrt{cm}",500000,-5e-2,5e-2);
    hdthetaturn = tfs->make<TH1D>("hdthetaturn",";dtheta for full turn, rad/#sqrt{cm}",500000,-5e-2,5e-2);
    hdthetaIWall= tfs->make<TH1D>("hdthetaIWall",";dtheta on inner wall, rad",100000,-1.,1.);
    hdthetaIWallNorm= tfs->make<TH1D>("hdthetaIWallNorm",";dtheta on inner wall normilized to perp as square, rad",100000,-1.,1.);

    hlenturn = tfs->make<TH1D>("hlenturn",";turn length, cm",1000,0.,500.);

    hdistres   =tfs->make<TH1D>("hdistres",";hit dist resol, mm",10000,-10,10);
    htargeteloss     =tfs->make<TH1D>("htargeteloss",";E loss at target, MeV",10000,0,10);

    treefit=tfs->make<TTree>("trfit","track fit params");
    treefit->Branch("startpar","HelixParams",&startpar);
    treefit->Branch("recopar","HelixParams",&recopar);
    treefit->Branch("hitinfo",&_tcellinfo);
    treefit->Branch("mchitinfo",&_mchitinfo);
    std::string fitStat ="fit/I:chi2/F:chi2in/F:chi2out/F:fitcon/F:radlen/F:ndof/I:niter/I:ninter/I:nactive/I:nsites/I:nhits/I:nhitstot/I:nturn/I";
    std::string fitData =":genmom/F:wgt/F:momtrackerin/F:momin/F:momout/F:t0/F:fitmom/F:fitmomerr/F:seedmom/F:fitmomout/F:fitmombeam/F:fitnturn/I:t0init/F:t0fit/F:errt0/F:iseed/I:iev/I";
    std::string recoData =":trkfrstHitPosX/F:trkfrstHitPosY/F:trkfrstHitPosZ/F:trk1Hitflgh/F:trk1SimHitflgh/F:trkto1SimHitdist/F:atSimHit/I";
    std::string genInfo =":genID/I:fromTargets/I:lowFitRange/F:hiFitRange/F:firsthitPLength/F:frstHitPosX/F:frstHitPosY/F:frstHitPosZ/F";
    std::string matData =":nHitSwire/I:nhitFwire/I:pathInSwire/F:pathInFwire/F:pathInGas/F:simPathInGas/F:simPathInGasFrstLoop/F:simTotPath/F";
    std::string momLossData =":simMomLossInVol/F:simMomLossInVolFrstLoop/F:simMomLossInIWFrstloop/F";
    //treefit->Branch("fitinfo",&recoinfo,"fit/I:chi2/F:chi2out/F:nhits/I:nhitstot/I:nturn/I:genmom/F:momin/F:momout/F:t0/F:fitmom/F:fitmomerr/F:fitnturn/I:t0init/F:t0fit/F:errt0/F:iseed/I:iev/I:genID/I:fromTargets/I:firsthitPLength/F:frstHitPosX/F:frstHitPosY/F:frstHitPosZ/F:nHitSwire/I:nhitFwire/I:pathInSwire/F:pathInFwire/F:pathInGas/F:simPathInGas/F:simPathInGasFrstLoop/F");
    treefit->Branch("fitinfo",&recoinfo,(fitStat+fitData+recoData+genInfo+matData+momLossData).c_str());
    treefit->Branch("eventinfo",&eventinfo,"elosstgt/F:ntgt/I");
    
  }


  bool kalFitOutUtils::FillMCInfo(art::Event const& event, std::vector<hitIndex> &strawhits, HelixTraj &seed){
    int iev=event.id().event();
    recoinfo.clear();
    findMCData(event);

    GeomHandle<ITracker> itr;
    CellGeometryHandle *itwp = itr->getCellGeometryHandle();
    //const Tracker& tracker = getTrackerOrThrow();
    ConditionsHandle<TrackerCalibrations> tcal("ignored");

    startpar=HelixParams(CLHEP::HepVector(5,0),CLHEP::HepSymMatrix(5,1));
    recopar=HelixParams(CLHEP::HepVector(5,0),CLHEP::HepSymMatrix(5,1));
    eventinfo.clear();

    art::Handle<GenParticleCollection> genParticles;
    event.getByLabel(_generatorModuleLabel, genParticles);

    // Steps in the stopping target
    art::Handle<StepPointMCCollection> stHitsHandle;
    event.getByLabel(_g4ModuleLabel,_targetStepPoints,stHitsHandle);
    StepPointMCCollection const& sthits(*stHitsHandle);

    art::Handle<StepPointMCCollection> vdHits;
    event.getByLabel(_g4ModuleLabel,_vdStepPoints,vdHits);

    art::Handle<SimParticleCollection> simParticles;
    event.getByLabel(_g4ModuleLabel, simParticles);

    art::Handle<mu2e::PtrStepPointMCVectorCollection> cellptr;
    event.getByLabel(_hitMakerModuleLabel,"StrawHitMCPtr",cellptr);
    PtrStepPointMCVectorCollection const* hits_mcptr = cellptr.product();
    if(hits_mcptr->size()==0) {
      cout<<"Empty ev = "<<event.id().event()<<endl;
      return false;
    }

    art::Handle<mu2e::StrawHitMCTruthCollection> cellmchits;
    art::Handle<mu2e::StrawHitCollection> cellhits;
    if(!(event.getByLabel(_hitMakerModuleLabel,cellhits)&&event.getByLabel(_hitMakerModuleLabel,cellmchits))){
      cout<<"no dc hits "
	  <<event.getByLabel(_hitMakerModuleLabel,cellhits)<<" "
	  <<event.getByLabel(_hitMakerModuleLabel,cellmchits)<<endl;
      return false;
    }

    //select first mc
    int firstmcstep=-1;
    double time0_min=1e10;
    time0_max=-1e10;
    for(unsigned int i=0;i<hits_mcptr->size();i++){
      const PtrStepPointMCVector&  mcptr=hits_mcptr->at(i);
      if(mcptr.size()){
	double tcell=mcptr.at(0)->time();
	GenId genid;
	if(mcptr.at(0)->simParticle()->isPrimary()&&
	   mcptr.at(0)->simParticle()->genParticle()) {
	  genid=mcptr.at(0)->simParticle()->genParticle()->generatorId();
	}
	if( genid==_genidcompare)
// 	  (genid==GenId::conversionGun
//            || genid==GenId::dioCzarnecki || genid==GenId::dioFlat || genid==GenId::dioE5 || genid==GenId::dioShankerWanatabe
//            || genid==GenId::cosmic || genid==GenId::cosmicDYB || genid==GenId::cosmicToy
//            || genid==GenId::PiCaptureCombined )
	{
	  if(tcell<time0_min){time0_min=tcell;firstmcstep=i;}
	  if(tcell>time0_max){time0_max=tcell;}
	}
      }
    }
    if(firstmcstep<0){
      cout<<"ev = "<<event.id().event()
          <<" Not electron in all cells, ncell="<<hits_mcptr->size();
      if(hits_mcptr->size()&&hits_mcptr->at(0).size())
	cout<<" id= "<< hits_mcptr->at(0).at(0)->trackId().asInt();
      cout<<endl;
      return false;
    }
    const PtrStepPointMCVector& mcptr=hits_mcptr->at(firstmcstep);
    
    double sumEDep=0.;
    if (true/*_fromOrigin*/) {
      gen_pos = mcptr.at(0)->simParticle()->genParticle()->position();
      //      gen_pos.setX(gen_pos.getX());
      gen_pos.setZ(gen_pos.getZ()+(12000-itr->z0()));
      gen_mom = mcptr.at(0)->simParticle()->genParticle()->momentum().vect();
      recoinfo.genmom=mcptr.at(0)->simParticle()->genParticle()->momentum().rho();
      if (_genidcompare==GenId::dioCzarnecki || 
          _genidcompare==GenId::dioFlat || 
          _genidcompare==GenId::dioE5 || 
          _genidcompare==GenId::dioShankerWanatabe) { recoinfo.wgt = DIOspectrum(recoinfo.genmom); }
      
      recoinfo.genID=_genidcompare.id();
      
      set<int> hitFoils;
      for ( StepPointMCCollection:: const_iterator i=sthits.begin();
	    i !=sthits.end(); ++i ){
	StepPointMC const& hit = *i;
	//	cout<<hit.trackId().asInt()<<" "<<hit.simParticle()->generatorIndex()<<" "<<hit.eDep()<<endl;
	if ( hit.trackId().asInt()!=1){
	  const SimParticle* sp=hit.simParticle().get();
	  while(sp&&sp->id().asInt()!=1) sp = sp->parent().get();
	  if(!sp||sp->id().asInt()!=1) continue;
	}
	sumEDep += hit.eDep();
	hitFoils.insert(hit.volumeId());
      }
      htargeteloss->Fill(sumEDep);
      eventinfo.elosstgt=sumEDep;
      eventinfo.ntgt=hitFoils.size();


      //Loop over vd hits
//      bool foundin = false;
//      bool foundout = false;
//      int vdIn(0), vdOut(0);
      double pIn(0);//, pOut(0);
      double intime = 1e15; 
//      double outtime = -1e15;
//      size_t inIndex(0), outIndex(0);
      if (vdHits.isValid()) {
        for(size_t i=0; i<vdHits->size(); ++i) {
          const StepPointMC& hit = (*vdHits)[i];
          //          std::cout<<i<<" "<<simParticles->at(hit.trackId()).fromGenerator()<<" vlmid "<<hit.volumeId()<<" must "<<VirtualDetectorId::IT_VD_InSurf<<" tim "<<hit.time()<<" mom "<<hit.momentum().mag()<<endl;
          if (simParticles->at(hit.trackId()).fromGenerator()) {
            int id = hit.volumeId();
            if (id == VirtualDetectorId::IT_VD_InSurf||
                id==VirtualDetectorId::IT_VD_EndCap_Front||
                id==VirtualDetectorId::IT_VD_EndCap_Back) {
              if (hit.time() < intime) {
                pIn = hit.momentum().mag();
                //foundin = true;
                //vdIn = id;
                //inIndex = i;
                intime = hit.time();
              }
              /*if (hit.time() > outtime) {
                pOut = hit.momentum().mag();
                foundout = true;
                vdOut = id;
                outIndex = i;
                outtime = hit.time();
              }*/
            }
          }
        }
      }

      recoinfo.momtrackerin=pIn;
      
    }


    time0=mcptr.at(0)->time();
    CLHEP::Hep3Vector pos=mcptr.at(0)->position();
    CLHEP::Hep3Vector mom=mcptr.at(0)->momentum();
    firstcell_pos=pos;

    double mom0=mom.mag();

    HepVector startPar(5);
    s0=0;
    double charge(-1.);

    if(_debugLvl>0) {
      cout<<"pos "<<pos<<" "<<mom<<" |p| "<<mom.mag()<<" targetloss "<<sumEDep<<endl;
    }

    TrkHelixUtils::helixFromMom( startPar, s0, 
				 HepPoint(pos.x(),pos.y(),pos.z()),
				 mom,
				 charge,
				 *_bfield);
    HepSymMatrix dummy(5,1); 
    dummy(1,1)=1.; dummy(2,2)=0.1*0.1;dummy(3,3)=1e-2*1e-2;
    dummy(4,4)=1.; dummy(5,5)=0.1*0.1;
    startpar=HelixParams(startPar,dummy);

    double helixdz=fabs(2*TMath::Pi()/startpar.omega()*startpar.tanDip());                                                                                                 
    int ndz=floor((pos.z()-startpar.z0())/helixdz);
    //std::cout<<ndz<<" "<<pos.z()<<" start "<<startpar.z0()<<" helixdz "<<helixdz<<" s0 "<<s0<<" ds0 "<<ndz*fabs(2*TMath::Pi()/startpar.omega()*sqrt(1+startpar.tanDip()*startpar.tanDip()))*(startpar.tanDip()>0?1.:-1.)<<std::endl;
    startpar.setZ0(startpar.z0()+helixdz*ndz);
    s0-=ndz*fabs(2*TMath::Pi()/startpar.omega()*sqrt(1+startpar.tanDip()*startpar.tanDip()))*(startpar.tanDip()>0?1.:-1.);
    seed=HelixTraj(startpar);

    double dt0atz0=seed.zFlight(0.0)/Constants::c;
    if(_debugLvl>1) {
      std::cout<<"s0 "<<s0<<" time0 "<<time0<<std::endl;
      seed.printAll(std::cout);
    }
    double mommin=mom0;
    

    unsigned int nn=cellhits->size();
    strawhits.clear();

    if(_debugLvl>1)std::cout<<"create hotlist"<<std::endl;
    std::vector<std::pair<double,unsigned int> > sorthits;
    hitflt.clear();
    for(unsigned int i=0;i<nn;++i){
        const mu2e::StepPointMC& stepmc=*cellptr->at(i).at(0);

	double time=stepmc.time();
	double ll=(time-time0)*Constants::c;
	sorthits.push_back(std::pair<double,unsigned int>(ll,i));
	hitflt.push_back(ll+s0-dt0atz0*Constants::c);//defined as relative to Z=0
    }

    std::sort(sorthits.begin(),sorthits.end());
    //sorted array
    double cosphidirprev=0.;
    double thetadirprev=0;
    double mmprev=0;
    double llprev=-1;
    double thetadirprevturn=0;
    double mmprevturn=0;
    double llprevturn=-1;
    int nhitstotal=0,nhitsperturn=0,nturn=1,nadd=0;
    int prevstep=-1;

    bool printevent=false;
    bool usehits=true;
    recoinfo.simTotPath=0.0;

    bool frstLoopNotstored = true;

    for(unsigned int i=0;i<nn;++i){
      int istr=sorthits[i].second;
      const mu2e::StepPointMC& stepmc=*cellptr->at(istr).at(0);
      GenId genid;
//      bool isOfInter = false;
      if(stepmc.simParticle()->isPrimary()&&
	 stepmc.simParticle()->genParticle()) {
	genid=stepmc.simParticle()->genParticle()->generatorId();
      }

      if(genid==_genidcompare){

          double mm=stepmc.momentum().mag();
	  double phidir=stepmc.momentum().phi()-stepmc.position().phi();
	  double thetadir=stepmc.momentum().theta();
	  double cosphidir=cos(phidir);
	  double dist=cellmchits->at(istr).driftDistance();
	  double time=stepmc.time();
	  double ll=(time-time0)*Constants::c;
	  if (ll>recoinfo.simTotPath) { recoinfo.simTotPath = ll; }
          if (llprev>0) {
            //                  std::cout<<"ihit "<<i<<" cellid "<<cellhits->at(istr).strawIndex()<<" abs path "<<ll<<" delta path "<<(ll-llprev)<<std::endl;
                  recoinfo.simPathInGas+=(ll-llprev);
                  recoinfo.simMomLossInVol+=(mmprev-mm);
          }

	  if((_oneturn&&cosphidirprev<0&&cosphidir>0)){
	    usehits=false;
	  }

	  if(prevstep>=0&&cosphidirprev<0&&cosphidir>0){
	    recoinfo.simPathInGas-=(ll-llprev);//(llprev-llprevturn);
	    recoinfo.simMomLossInVol-=(mmprev-mm);
	    if(mm>90){
	      helossIWall->Fill(mmprev-mm);
	      double perpnorm=(fabs(cosphidir)*sin(thetadir)+fabs(cosphidirprev)*sin(thetadirprev))/2;
	      helossIWallNorm->Fill((mmprev-mm)*perpnorm);
	      hdthetaIWall->Fill(thetadir-thetadirprev);
	      hdthetaIWallNorm->Fill((thetadir-thetadirprev)*sqrt(perpnorm));

	      if(nhitsperturn>2){
		hdedxturn->Fill((mmprevturn-mmprev)/((llprev-llprevturn)/cm));
		hdthetaturn->Fill((thetadirprev-thetadirprevturn)/sqrt((llprev-llprevturn)/cm));
		hlenturn->Fill((llprev-llprevturn)/cm);
		if (frstLoopNotstored) {
		        recoinfo.simPathInGasFrstLoop=(llprev-llprevturn);
		        recoinfo.simMomLossInVolFrstLoop=(mmprevturn-mmprev);
		        recoinfo.simMomLossInIWFrstloop=(mmprev-mm);
		        frstLoopNotstored = false;
		}
	      }
	    }

	    if(usehits){
	      hnhits->Fill(nhitsperturn);
	      nhitsperturn=0;
	      nturn++;
	    }	    
	    thetadirprevturn=thetadir;
	    llprevturn=ll;
	    mmprevturn=mm;
	  }

	  if(prevstep>=0&&mm>90&&fabs((ll-llprev)/cm)<10){
	    hdedx->Fill((mmprev-mm)/((ll-llprev)/cm));	  
	    hdtheta->Fill((thetadir-thetadirprev)/sqrt((ll-llprev)/cm));
	  }
	  cosphidirprev=cosphidir;
	  thetadirprev=thetadir;
	  llprev=ll;
	  mmprev=mm;
	  if(prevstep<0){
	    thetadirprevturn=thetadir;
	    llprevturn=ll;
	    mmprevturn=mm;
	  }
	  prevstep=i;

	  //	  std::cout<<nhitstotal<<" "<<_oneturn<<" "<<cosphidirprev<<" "<<cosphidir<<std::endl;
	  if(usehits){
	    if(mm<mommin) mommin=mm;
	    nhitsperturn++;
	    nhitstotal++;
	    nadd++;
	    strawhits.push_back(hitIndex(istr,dist<0?1:-1));
	  }
      }
    }
    hnhits->Fill(nhitsperturn);
    hnhitstotal->Fill(nhitstotal);
    hnturns->Fill(nturn);
    heloss->Fill(mom0-mommin);    

    for(unsigned int ii=0;ii<nn;++ii){
      int istr=sorthits[ii].second;
      const mu2e::StepPointMC& stepmc=*cellptr->at(istr).at(0);
      double time=stepmc.time();
      double mm=stepmc.momentum().mag();
      double ll=(time-time0)*Constants::c;
      double dist=cellmchits->at(istr).driftDistance();
      unsigned int dcell_id=cellhits->at(istr).strawIndex().asUint();
      itwp->SelectCellDet(dcell_id);
      const StrawHit& sh = cellhits->at(istr);
      //const Straw& straw = tracker.getStraw(sh.strawIndex());

      int layer=itwp->GetSuperLayer();//dcell_id/10000;
      int cell=itwp->GetWire();//dcell_id%10000;
      if(_debugLvl>1||printevent){
	GenId genid;
	if(stepmc.simParticle()->isPrimary()&&
	   stepmc.simParticle()->genParticle())
	  genid=stepmc.simParticle()->genParticle()->generatorId();

	std::cout<<iev<<" "<<istr<<" lay,side,cell "<<layer<<" "<<itwp->isUpStream()<<" "<<cell
		 <<" time "<<time<<" hittime "<<sh.time()<<" dist "<<dist
		 <<" flt "<<ll<<" id "<<stepmc.trackId().asInt()<<" eloss "<<mom0-mm
		 <<" xyz "<<stepmc.position().x()<<" "<<stepmc.position().y()<<" "<<stepmc.position().z()
		 <<" genid "<<genid.name()<<std::endl;
      }
      T2D t2d;
      double onWirePropTime = (itwp->GetCellHalfLength()-stepmc.position().z()/cos(itwp->GetWireEpsilon()) )/tcal->SignalVelocity(cellhits->at(istr).strawIndex());
      tcal->TimeToDistance(cellhits->at(istr).strawIndex(),(sh.time()-(time + onWirePropTime) ), stepmc.momentum(), t2d);
      hdistres->Fill(t2d._rdrift-fabs(dist));
      //StrawHit sh2(sh.strawIndex(),distup/0.035+time/*+1.37459e-01/0.035 (wire signal prop)*/,sh.dt(),sh.energyDep());

    }

    //recoinfo.clear();
    recoinfo.iev=iev;
    recoinfo.nhits=-1;
    recoinfo.nhitstot=nhitstotal;
    recoinfo.nturn=nturn;
    recoinfo.momin=mom0;
    recoinfo.momout=mommin;    

    recoinfo.t0=time0-s0/Constants::c+dt0atz0;
    recoinfo.firsthitPLength=s0;
    recoinfo.frstHitPosX=firstcell_pos.x();
    recoinfo.frstHitPosY=firstcell_pos.y();
    recoinfo.frstHitPosZ=firstcell_pos.z();

    return true;
  }

  // find the MC truth objects in the event and set the local cache
    bool
    kalFitOutUtils::findMCData(const art::Event& evt) {
            _mcdata.clear();
            //_mchitsums.clear();
            art::Handle<StrawHitMCTruthCollection> truthHandle;
            if(evt.getByLabel(_hitMakerModuleLabel,truthHandle))
                    _mcdata._mcstrawhits = truthHandle.product();
            // Get the persistent data about pointers to StepPointMCs
            art::Handle<PtrStepPointMCVectorCollection> mchitptrHandle;
            if(evt.getByLabel(_hitMakerModuleLabel,"StrawHitMCPtr",mchitptrHandle))
                    _mcdata._mchitptr = mchitptrHandle.product();
            // Get the persistent data about the StepPointMCs, from the tracker and the virtual detectors
            art::Handle<StepPointMCCollection> mctrackerstepsHandle;
            if(evt.getByLabel(_g4ModuleLabel,"tracker",mctrackerstepsHandle))
                    _mcdata._mcsteps = mctrackerstepsHandle.product();
            art::Handle<StepPointMCCollection> mcVDstepsHandle;
            if(evt.getByLabel(_g4ModuleLabel,"virtualdetector",mcVDstepsHandle))
                    _mcdata._mcvdsteps = mcVDstepsHandle.product();
            art::Handle<SimParticleCollection> simParticlesHandle;
            if(evt.getByLabel(_g4ModuleLabel,simParticlesHandle))
                    _mcdata._simparts = simParticlesHandle.product();
            /*if( _mcdata.good()){
                    // mc track patermeter info
                    mcTrkInfo();
                    // fill hit summary
                    fillMCHitSummary();
                    // count # of conversion straw hits
                    _strawhits = 0;
                    _nchits = 0;
                    art::Handle<mu2e::StrawHitCollection> strawhitsH;
                    if(evt.getByLabel(_hitMakerModuleLabel,strawhitsH))
                            _strawhits = strawhitsH.product();
                    if(_strawhits != 0){
                            unsigned nstrs = _strawhits->size();
                            for(unsigned istr=0; istr<nstrs;++istr){
                                    const std::vector<TrkSum>& mcsum = mcHitSummary(istr);
                                    if(mcsum.size()>0){
                                            bool conversion = (mcsum[0]._pdgid == 11 && mcsum[0]._gid == 2);
                                            if(conversion){
                                                    ++_nchits;
                                            }
                                    }
                            }
                    }
                    return true;
            }*/
            return false;
    }

  void kalFitOutUtils::hitsDiag(std::vector<const TrkCellHit*> const& hits) {
    _tcellinfo.clear();
//    _tainfo.clear();
//    _ncactive = 0;
//// find the arcs
//    std::vector<TrkArc> arcs;
//    findArcs(hits,arcs);
//    _narcs = arcs.size();
//// loop over arcs
//    for(size_t iarc=0;iarc < arcs.size(); ++iarc){
//      TrkArcInfo tainfo;
//      TrkArc const& arc = arcs[iarc];
//      tainfo._narctsh = arc._ntsh;
//      tainfo._narcactive = arc._nactive;
//      tainfo._arctshlen = hits[arc._end]->fltLen() - hits[arc._begin]->fltLen();
//      tainfo._arcactivelen = hits[arc._endactive]->fltLen() - hits[arc._beginactive]->fltLen();
//      _tainfo.push_back(tainfo);
//    }

    if (_debugLvl>0) {_mchitinfo.clear();}
 // loop over hits
    for(size_t itsh=0;itsh<hits.size();++itsh){
      const TrkCellHit* tsh = hits[itsh];
      if(tsh != 0){
        HitInfo mchitinfo;

        TrkCellHitInfo tcellinfo;
        tcellinfo._active = tsh->isActive();
        tcellinfo._usable = tsh->usability();
        tcellinfo._Slayer = tsh->GetSuperLayer();
        tcellinfo._layer = tsh->GetCelRing();
        tcellinfo._cell = tsh->GetCell();
        static HepPoint origin(0.0,0.0,0.0);
        CLHEP::Hep3Vector hpos = tsh->hitTraj()->position(tsh->hitLen()) - origin;
        tcellinfo._z = hpos.z();
        tcellinfo._phi = hpos.phi();
        tcellinfo._rho = hpos.perp();
        double resid,residerr;
        if(tsh->resid(resid,residerr,true)){
          tcellinfo._resid = resid;
          tcellinfo._residerr = residerr;
        } else {
          tcellinfo._resid = tcellinfo._residerr = -1000.;
        }
        tcellinfo._rdrift = tsh->driftRadius();
        tcellinfo._rdrifterr = tsh->driftRadiusErr();
        tcellinfo._trklen = tsh->fltLen();
        tcellinfo._hlen = tsh->hitLen();
        tcellinfo._ht = tsh->hitT0()._t0;
        tcellinfo._tddist = tsh->timeDiffDist();
        tcellinfo._tdderr = tsh->timeDiffDistErr();
        tcellinfo._ambig = tsh->ambig();
        if(tsh->poca() != 0)
          tcellinfo._doca = tsh->poca()->doca();
        else
          tcellinfo._doca = -100.0;
        tcellinfo._exerr = tsh->extErr();
        tcellinfo._penerr = tsh->penaltyErr();
        tcellinfo._t0err = tsh->t0Err();
//// arc information
//        int iarc = findArc(itsh,arcs);
//        tcellinfo._iarc = iarc;
//        if(iarc >= 0){
//          TrkArc const& arc = arcs[iarc];
//          tcellinfo._iarchit = itsh-arc._begin;
//          tcellinfo._architlen = tsh->fltLen() - hits[arc._beginactive]->fltLen();
//          if(itsh > arc._begin)
//            tcellinfo._gaplow = tsh->fltLen() - hits[itsh-1]->fltLen();
//          else
//            tcellinfo._gaplow = -1;
//          if(itsh < arc._end)
//            tcellinfo._gaphi = hits[itsh+1]->fltLen()-tsh->fltLen();
//          else
//            tcellinfo._gaphi = -1;
//        } else {
//          tcellinfo._iarchit = -1;
//          tcellinfo._architlen = -1;
//          tcellinfo._gaplow = -1;
//          tcellinfo._gaphi = -1;
//        }
        // MC information
        PtrStepPointMCVector const& mcptr(_mcdata._mchitptr->at(tsh->index()));
        tcellinfo._mcn = mcptr.size();
        if(_mcdata._mcsteps != 0){
          //std::vector<TrkSum> mcsum;
          //KalFitMC::fillMCHitSum(mcptr,mcsum);
          //tcellinfo._mcnunique = mcsum.size();
          //tcellinfo._mcppdg = mcsum[0]._pdgid;
          //tcellinfo._mcpgen = mcsum[0]._gid;
          //tcellinfo._mcpproc = mcsum[0]._pid;
// first hit is the one that set t0
          tcellinfo._mcht = mcptr[0]->time();
          art::Ptr<SimParticle> sp = mcptr[0]->simParticle();
          tcellinfo._mcpdg = sp->pdgId();
          tcellinfo._mcproc = sp->creationCode();
          tcellinfo._mcgen = -1;
          if(sp->genParticle().isNonnull())
            tcellinfo._mcgen = sp->genParticle()->generatorId().id();
// find the step midpoint
          Hep3Vector start = mcptr[0]->position();
          Hep3Vector dir = mcptr[0]->momentum().unit();
          Hep3Vector end = start + dir*mcptr[0]->stepLength();
          Hep3Vector mid = 0.5*(start+end);
          Hep3Vector mcsep = mid-tsh->straw().getMidPoint();
          Hep3Vector mcperp = (dir.cross(tsh->straw().getDirection())).unit();
          double dperp = mcperp.dot(mcsep);
          tcellinfo._mcdist = fabs(dperp);
          tcellinfo._mcambig = dperp > 0 ? -1 : 1; // follow TrkPoca convention
          tcellinfo._mclen = mcsep.dot(tsh->straw().getDirection());

          if (_debugLvl>0) {
                  StrawHitMCTruth const& mctrth(_mcdata._mcstrawhits->at(tsh->index()));
                  mchitinfo._mcpdg = sp->pdgId();
                  mchitinfo._mcgen = -1;
                  mchitinfo._mcproc = sp->creationCode();
                  mchitinfo._mcedep = 0.0;
                  for (size_t istep=0; istep<mcptr.size(); ++istep) {
                          mchitinfo._mcedep += mcptr.at(istep)->eDep();
                  }
                  mchitinfo._mcdisttomid = mctrth.distanceToMid();
                  mchitinfo._mcpos = tsh->cellHandle()->GetCellCenter()+mctrth.distanceToMid()*tsh->cellHandle()->GetWireDirection();

                  mchitinfo._edep = tsh->strawHit().energyDep();
                  mchitinfo._disttomid = tsh->timeDiffDist();
                  mchitinfo._disttomid_err = tsh->timeDiffDistErr();
                  mchitinfo._pos = tsh->cellHandle()->GetCellCenter()+tsh->timeDiffDist()*tsh->cellHandle()->GetWireDirection();
                  mchitinfo._Slayer = tsh->GetSuperLayer();
                  mchitinfo._layer = tsh->GetCelRing();
                  mchitinfo._cell = tsh->GetCell();

          }
        }
        _tcellinfo.push_back(tcellinfo);
        if (_debugLvl>0) { _mchitinfo.push_back(mchitinfo); }
// count active conversion hits
        //if(tcellinfo._mcgen==2&&tsh->isActive())++_ncactive;
      }
    }
  }

//  void
//  kalFitOutUtils::findArcs(std::vector<const TrkCellHit*> const& cells, std::vector<TrkArc>& arcs) const {
//    arcs.clear();
//// define an initial arc
//    size_t icell(0);
//    while(icell < cells.size()){
//      int igap(0);
//      TrkArc arc(icell);
//      bool firstactive(false);
//      do {
//        arc._end = icell;
//        ++arc._ntsh;
//        if(cells[icell]->isActive()){
//          arc._endactive = icell;
//          ++arc._nactive;
//          if(!firstactive){
//            arc._beginactive = icell;
//            firstactive = true;
//          }
//        }
//// advance to the next straw and compute the gap
//        ++icell;
//        if(icell< cells.size())igap = cells[icell]->straw().id().getDevice()-cells[arc._end]->straw().id().getDevice();
//// loop while the gap is small
//      } while(icell < cells.size() && igap <= _maxarcgap);
//// end of an arc: record it
//      arcs.push_back(arc);
//// update the arc to start with this new straw
//      arc = TrkArc(icell);
//    }
//  }
//
//  int
//  kalFitOutUtils::findArc(size_t itsh,std::vector<TrkArc>& arcs ) {
//    for(size_t iarc=0;iarc<arcs.size();++iarc){
//      if(itsh >= arcs[iarc]._begin && itsh <= arcs[iarc]._end)
//        return iarc;
//    }
//    return -1;
//  }

  double kalFitOutUtils::DIOspectrum(double ee) {
    double mal(25133);
//    double mmu(105.654);
    double emu(105.194);
//    double emue(104.973);
//    double me(0.511);
    double a5(8.6434e-17);
    double a6(1.16874e-17);
    double a7(-1.87828e-19);
    double a8(9.16327e-20);
    double delta = emu - ee - ee*ee/(2*mal);
    double deltaPO5 = delta*delta;
    deltaPO5*=deltaPO5;
    deltaPO5*=delta;
    double deltaPO6 =  deltaPO5*delta;
    double deltaPO7 =  deltaPO6*delta;
    double deltaPO8 =  deltaPO7*delta;
    return a5*deltaPO5 + a6*deltaPO6 + a7*deltaPO7 + a8*deltaPO8;
  }
  
  void kalFitOutUtils::FillHistos(KalFitResult& kres, HelixTraj &seed, int iseed){

    GeomHandle<ITracker> itr;                                                                                                                                                                       

    KalRep *kalrep = kres._krep;
    TrkErrCode err= kres._fit;

    recoinfo.clearrec();
    recoinfo.nhits=kres._krep?kres._krep->hotList()->nActive():-1;
    recoinfo.iseed=iseed;
    recoinfo.t0init=kres._tdef.t0()._t0;

    if(!err.success()) {
      treefit->Fill();
      //      delete kalrep;
      return;
    }

    recoinfo.nHitSwire=0; recoinfo.nhitFwire=0;
    recoinfo.pathInSwire=0.0; recoinfo.pathInFwire=0.0; recoinfo.pathInGas=0.0;
    const TrkHotList* hots = kalrep->hotList();
    std::vector<const TrkCellHit*> hits;
    hits.reserve(hots->nHit());
    for(TrkHotList::hot_iterator ihot=hots->begin();ihot != hots->end();++ihot){
            const TrkCellHit* hit = dynamic_cast<const TrkCellHit*>(ihot.get());
            if(hit != 0){
                    hits.push_back(hit);
                    if (hit->hitSwire()) { recoinfo.nHitSwire++; recoinfo.pathInSwire+=hit->getSwirePath(); }
                    if (hit->hitFwireBot()) { recoinfo.nhitFwire++; recoinfo.pathInFwire+=hit->getFwirePathBot(); }
                    if (hit->hitFwireSide()) { recoinfo.nhitFwire++; recoinfo.pathInFwire+=hit->getFwirePathSide(); }
                    if (hit->hitFwireTop()) { recoinfo.nhitFwire++; recoinfo.pathInFwire+=hit->getFwirePathTop(); }
                    recoinfo.pathInGas+=(hit->exitDeltaPath()-hit->entryDeltaPath());
            }
    }
    //std::sort(hits.begin(),hits.end(),devicecomp());
    hitsDiag(hits);

    if(_debugLvl>1){
      const TrkHotList* hots = kalrep->hotList();
      for(TrkHotList::hot_iterator ihot=hots->begin();ihot != hots->end();++ihot){
	const TrkCellHit* hit = dynamic_cast<const TrkCellHit*>(ihot.get());
	if(hit != 0){
	  double resid=-1000,residerr=-1000;
	  hit->resid(resid,residerr,true);
	  static HepPoint origin(0.0,0.0,0.0);
	  CLHEP::Hep3Vector hpos = hit->hitTraj()->position(hit->hitLen())-origin;
	  cout<<"hitleng "<<hit->hitLen()<<" active "<<hit->isActive()<<" usable "<<hit->usability()
	      <<" ambig "<<hit->ambig()<<" resid "<<resid<<" residErr "<<residerr
	      <<" driftR "<<hit->driftRadius()<<" driftRErr "<<hit->driftRadiusErr()
	      <<" xyz "<<hpos<<" hitindex "<<hit->strawHit().strawIndex().asUint()<<endl;
	}
      }
    }
    
    if(_debugLvl>2){
      std::cout<<"fit kalrep"<<std::endl;
      kalrep->printAll(std::cout);
    }

    recoinfo.ndof = kalrep->nDof();
    recoinfo.nactive=kalrep->nActive();
    recoinfo.niter=kalrep->iterations();
    recoinfo.ninter=kalrep->intersections();
    recoinfo.nsites=kalrep->siteList().size();

    double chi2=kalrep->chisq();
    double chi2in=kalrep->chisquared(trkIn);
    double chi2out=kalrep->chisquared(trkOut);
    if(_debugLvl) {
      cout<<"chi2 in "<<chi2in<<" out "<<chi2out<<endl;
      //      cout<<"hello err="<<err<<endl;
      // seed.printAll(std::cout);
      cout<<"true t0 helix apprx at z=0 "<<recoinfo.t0
	  <<" (t0 at first hit "<<time0<<" )fitted at z=0 of helix niters "<<kres._nt0iter
	  <<" initt0 "<<kres._tdef.t0()._t0<<" +- "<<kres._tdef.t0()._t0err
	  <<" after fit "<<kalrep->t0()._t0<<" +- "<<kalrep->t0()._t0err<<endl;
    }
    recoinfo.t0fit=kalrep->t0()._t0;
    recoinfo.errt0=kalrep->t0()._t0err;

    double Bval = _bfield->bFieldNominal();
    std::cout<<"Nominal field "<<Bval<<" "<<std::endl;

    double fltlen=kalrep->firstHit()->globalLength();
    double fltlenEnd=kalrep->lastHit()->globalLength();
    recopar = kalrep->helix(fltlen);

    const TrkDifPieceTraj& recotraj = kalrep->pieceTraj();


    double dotin  = TMath::Cos(kalrep->momentum(fltlen).phi()-recotraj.position(fltlen).phi());
    double dotout = TMath::Cos(kalrep->momentum(fltlenEnd).phi()-recotraj.position(fltlenEnd).phi());
    double fnturns=(fltlenEnd-fltlen)/sqrt(1+recopar.tanDip()*recopar.tanDip())*fabs(recopar.omega())/(2*TMath::Pi());
    recoinfo.fitnturn=int(1+fnturns+int(dotin<0&&dotout>0)*1.+(dotin*dotout>0)*0.5);

    TrkLineTraj firstp(HepPoint(firstcell_pos.x(),firstcell_pos.y(),firstcell_pos.z()), 
		       Hep3Vector(0, 0, 1), 1.e-3);
    double simhitfltlen=fltlen+(firstcell_pos.z()-recotraj.position(fltlen).z())*sqrt(1+1./recopar.tanDip()*1./recopar.tanDip())*(recopar.tanDip()>0?1.:-1.);
    TrkPoca recpoca1(recotraj,simhitfltlen,false,firstp,0,true, 1e-12);
    double firstplen=fltlen;
    if(recpoca1.status().success()){
      firstplen = recpoca1.flt1();
    }

    recoinfo.lowFitRange=kalrep->lowFitRange();
    recoinfo.hiFitRange=kalrep->hiFitRange();
    recoinfo.trk1Hitflght=fltlen;
    recoinfo.trk1SimHitflght=recpoca1.flt1();
    recoinfo.trkto1SimHitdist=recpoca1.doca();
    recoinfo.trkfrstHitPosX=recotraj.position(firstplen).x();
    recoinfo.trkfrstHitPosY=recotraj.position(firstplen).y();
    recoinfo.trkfrstHitPosZ=recotraj.position(firstplen).z();
    recoinfo.atSimHit=recpoca1.status().success();
    if(_debugLvl) {
            std::cout<<"first cell len "<<fltlen<<" first gen cell len "<<firstplen<<" estimate "<<simhitfltlen
                            <<" lastfltlen "<<fltlenEnd<<" "<<recpoca1.status().success()<<" "<<fnturns<<" "<<dotin<<" "<<dotout<<" "<<recoinfo.fitnturn<<std::endl;
    }
    fltlen=firstplen;
    recopar = kalrep->helix(fltlen);

    //set to first crossing of beam below first hit 
    double rechelixdz=fabs(2*TMath::Pi()/recopar.omega()*recopar.tanDip());
    recopar.setZ0(recopar.z0()+rechelixdz*floor((recoinfo.trkfrstHitPosZ-recopar.z0())/rechelixdz));

    

    HepVector recoparams = recopar.params();
    HepSymMatrix recocovar = recopar.covariance();


    CLHEP::Hep3Vector fitmom = kalrep->momentum(fltlen);
    BbrVectorErr momerr = kalrep->momentumErr(fltlen);
    recoinfo.fitmom = fitmom.mag();
    Hep3Vector momdir = fitmom.unit();
    HepVector momvec(3);
    for(int icor=0;icor<3;icor++)
      momvec[icor] = momdir[icor];
    recoinfo.fitmomerr = sqrt(momerr.covMatrix().similarity(momvec));
    CLHEP::Hep3Vector seedmom = TrkMomCalculator::vecMom(*(kalrep->seed()),*_bfield/*kalrep->bField()*/,0.0);
    recoinfo.seedmom = seedmom.mag();

    HepVector hvpos(3),hvmom(3);
    HepSymMatrix ehpos(3),ehmom(3);
    HepMatrix ehposmom(3,3);
    kalrep->getAllWeights(fltlen,hvpos,hvmom,ehpos,ehmom,ehposmom);

    TrkLineTraj zaxis(HepPoint(0, 0, -10000), Hep3Vector(0, 0, 1),10000-itr->zHalfLength());

    double zbeam_est=recopar.z0()+rechelixdz*floor((-itr->zHalfLength()-recopar.z0())/rechelixdz);
    double fltlenbeam_est=fltlen+(zbeam_est-recotraj.position(fltlen).z())*sqrt(1+1./(recopar.tanDip()*recopar.tanDip()))*(recopar.tanDip()>0?1.:-1.);

    TrkPoca recpoca(recotraj,fltlenbeam_est,zaxis,zbeam_est+10000, 1e-12);
    double fltlenbeam=fltlenbeam_est;
    if(recpoca.status().success())
      fltlenbeam = recpoca.flt1();
    HelixParams recobeam = kalrep->helix(fltlenbeam);
    HepVector recoparams_beam = recobeam.params();
    HepSymMatrix recocovar_beam = recobeam.covariance();
    double momem_beam=Constants::c*1.0e-3*Bval/recoparams_beam[2]*sqrt(1.+recoparams_beam[4]*recoparams_beam[4]);

//    recoinfo.fitmombeam=momem_beam;
    if (recoinfo.lowFitRange<recoinfo.trk1Hitflght) {
            recoinfo.fitmombeam=kalrep->momentum(recoinfo.lowFitRange).mag();
    } else {
            recoinfo.fitmombeam=-1;
    }

    double helixdz=2*TMath::Pi()/recobeam.omega()*recobeam.tanDip();
    kalrep->getAllWeights(fltlenbeam,hvpos,hvmom,ehpos,ehmom,ehposmom);
    if(_debugLvl) {
      cout<<"gen pos "<<gen_pos<<" gen mom "<<gen_mom
	  <<" reco ("<<hvpos[0]<<","<<hvpos[1]<<","<<hvpos[2]<<") ("<<hvmom[0]<<","<<hvmom[1]<<","<<hvmom[2]<<")"
	  <<" helix dz "<<helixdz<<" fltbeam "<<fltlenbeam<<" fltbeam_est "<<fltlenbeam_est
          <<" genpart "<<GenId(recoinfo.genID).name()<<endl;
    }

    HepVector recoparamsEnd =  kalrep->helix(fltlenEnd).params();
    HepVector startPar=startpar.params();

    double momem=Constants::c*1.0e-3*Bval/recoparams[2]*sqrt(1.+recoparams[4]*recoparams[4]);
    double momemEnd=Constants::c*1.0e-3*Bval/recoparamsEnd[2]*sqrt(1.+recoparamsEnd[4]*recoparamsEnd[4]);
    double momem0=Constants::c*1.0e-3*Bval/startPar[2]*sqrt(1.+startPar[4]*startPar[4]);
    hmom->Fill(momem-momem0);
    if(_debugLvl) {
      cout<<"Fit mom "<<momem<<" min "<<momemEnd<<" beam "<<momem_beam
	  <<" mom0 "<<momem0<<" min "<<recoinfo.momout<<" beam "<<gen_mom.mag()<<" trackerin "<<recoinfo.momtrackerin<<endl;
      recopar.printAll(std::cout);
    }
    recoinfo.fitmomout=momemEnd;

    //recoinfo.fit=1;
    recoinfo.fit=kalrep->fitStatus().success();
    recoinfo.chi2=chi2;
    recoinfo.chi2in=chi2in;
    recoinfo.chi2out=chi2out;
    recoinfo.fitcon=kalrep->chisqConsistency().significanceLevel();
    recoinfo.radlen=kalrep->radiationFraction();

    treefit->Fill();

    hchi2->Fill(chi2/(recoinfo.nhits-5));
    hprob->Fill(TMath::Prob(chi2,(recoinfo.nhits-5)));
    hchi2_0->Fill(chi2out/(recoinfo.nhits-5));

    double maxpull=0;
    for(int i=0;i<5;i++){
      hpars[i]->Fill(recoparams[i]-startPar[i]);
      hdev->Fill(i,(recoparams[i]-startPar[i]));
      double pull=(recoparams[i]-startPar[i])/TMath::Sqrt(recocovar[i][i]);
      hpull->Fill(i,pull);
      if(fabs(pull)>maxpull) maxpull=fabs(pull);
    }
    
    if(maxpull>100){
      cout<<"problem, big pull value "<<maxpull<<endl;
      nbadfit++;
    }
    
    //delete kalrep;
    //delete hotlist;
    //delete retval;

  }
  
  void kalFitOutUtils::finalizeHistos()
  {
    art::ServiceHandle<art::TFileService> tfs;
    cout<<"nbad "<<nbadfit<<std::endl;

    const char *titlenames[]={"pull d0;#delta parameter/#sigma",
			      "pull #phi0;#delta parameter/#sigma",
			      "pull #omega;#delta parameter/#sigma",
			      "pull z0;#delta parameter/#sigma",
			      "pull tan;#delta parameter/#sigma",
			      "chi2/ndf;#chi^{2}/ndf"};

    TCanvas *c1=tfs->make<TCanvas>(/*"c1pull"*/);
    c1->Divide(2,3,0.001,0.001);
    for(int i=0;i<5;i++){
      c1->cd(i+1);
      TH1D *hp=hpull->ProjectionY(Form("_py%i",i),i+1,i+1);
      hp->SetTitle(titlenames[i]);
      hp->SetAxisRange(-10,10);
      hp->Fit("gaus");
    }
    c1->cd(6);
    hchi2->SetAxisRange(0,10);
    hchi2->Draw();
    c1->Write();
    TCanvas *c2=tfs->make<TCanvas>(/*"c1par"*/);
    c2->Divide(2,3,0.001,0.001);

    for(int i=0;i<5;i++){
      c2->cd(i+1);
      hpars[i]->Fit("gaus");
    }
    c2->Write();
    hmom->SetAxisRange(-4,4);
    hmom->Fit("gaus","","",-0.5,0.75);
  }
}  // end namespace mu2e
