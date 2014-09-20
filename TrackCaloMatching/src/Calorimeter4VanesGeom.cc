//
// $Id: Calorimeter4VanesGeom.cc,v 1.15 2014/09/20 18:04:22 murat Exp $
// $Author: murat $
// $Date: 2014/09/20 18:04:22 $
//
// Original author G. Pezzullo & G. Tassielli
//

#include "TrackCaloMatching/inc/Calorimeter4VanesGeom.hh"
#include "BaBar/Constants.hh"
#include <cmath>




namespace mu2e {


//-----------------------------------------------------------------------------
// Construct from a transform.
//-----------------------------------------------------------------------------
  Calorimeter4VanesGeom::Calorimeter4VanesGeom(){
    art::ServiceHandle<GeometryService> geom;
    GeomHandle<Calorimeter> cg;
    _norm.setX(0.0);
    _norm.setY(1.0);
    _norm.setZ(0.0);

    _dU = (cg->caloGeomInfo().crystalHalfLength() + cg->caloGeomInfo().wrapperThickness() ) ;
    _dAPD = cg->caloGeomInfo().roHalfThickness();
    _vaneHalfThickness = _dU + _dAPD;
    
    if(geom->hasElement<VaneCalorimeter>()){
      GeomHandle<VaneCalorimeter> cgVanes;
      _dR = (cgVanes->caloGeomInfo().crystalHalfTrans()+cgVanes->caloGeomInfo().wrapperThickness() + cgVanes->caloGeomInfo().shellThickness() ) * cgVanes->nCrystalR();//*0.5;
      _dZ = (cgVanes->caloGeomInfo().crystalHalfTrans()+cgVanes->caloGeomInfo().wrapperThickness() + cgVanes->caloGeomInfo().shellThickness() ) * cgVanes->nCrystalZ();//*0.5;
    }else if(geom->hasElement<DiskCalorimeter>()){
      _dR = 0.0;
      _dZ = 0.0;
    }
    
               
    _solenoidOffSetX = geom->config().getDouble("mu2e.solenoidOffset");//3904.;//[mm]
    _solenoidOffSetZ = -geom->config().getDouble("mu2e.detectorSystemZ0");//-10200.;
    
    _ZfrontFaceCalo = cg->origin().z() + _solenoidOffSetZ - _dZ;//(cg->nCrystalZ() + 1.0)*cg->crystalHalfSize();
    _ZbackFaceCalo = cg->origin().z() + _solenoidOffSetZ + _dZ;//(cg->nCrystalZ() + 1.0)*cg->crystalHalfSize(); 
    
    if(geom->hasElement<DiskCalorimeter>()){
      GeomHandle<DiskCalorimeter> cgDisks;
//       _ZfrontFaceCalo -= _vaneHalfThickness*2.0 + cgDisks->diskSeparation(1);
//       _ZbackFaceCalo += cgDisks->diskSeparation(1) + _vaneHalfThickness*2.0 + 10.0;
      _ZfrontFaceCalo += cgDisks->diskSeparation(0) - 10.;
      _ZbackFaceCalo  += cgDisks->diskSeparation(1) + 2*cg->caloGeomInfo().crystalHalfLength() + 10.0;
    }
    
    if(geom->hasElement<VaneCalorimeter>()){ 
      GeomHandle<VaneCalorimeter> cg;
      _innerRadius = cg->innerRadius();
      _outherRadius = cg->outherRadius();
      _nVanes = cg->nVane();
    }else if(geom->hasElement<DiskCalorimeter>()){
      GeomHandle<DiskCalorimeter> cg;
      _innerRadius = cg->disk(0).innerRadius();
      _outherRadius = cg->disk(0).outerRadius();
      _nVanes = cg->nDisk();
    }
  }



  CaloVolumeElem* Calorimeter4VanesGeom::vane(int& i){

    GeomHandle<VaneCalorimeter> cg;
    Vane vane = cg->vane(i);

    CaloVolumeType *caloVane = new CaloVolumeType("caloVane", 0);

    //face 0
    const Hep3Vector vec_0(0.0,_vaneHalfThickness , 0.0);
    const Hep3Vector norm_0(0., 0., 1.);

    const HepRotation rot_0(norm_0,0.0);
    const HepTranslation tras_0(vec_0);
    const HepTransformation  transf_0(rot_0, tras_0);

    CaloSurface *face_0 = new CaloSurface(transf_0, _dR, _dZ, _norm);
    caloVane->AddSide(face_0);

    unsigned int idFace_0 = 0;

    caloVane->AddSideCorner( new SurfacePoint(-_dR, _dZ) , idFace_0);
    caloVane->AddSideCorner( new SurfacePoint(_dR, _dZ)  , idFace_0);
    caloVane->AddSideCorner( new SurfacePoint(_dR, -_dZ) , idFace_0);
    caloVane->AddSideCorner( new SurfacePoint(-_dR, -_dZ), idFace_0);

    //face 1
    const Hep3Vector vec_1(_dR, 0.0, 0.0 );
    const Hep3Vector norm_1(0., 0., 1.);

    const HepRotation rot_1(norm_1, -0.5*Constants::pi);
    const HepTranslation tras_1(vec_1);
    const HepTransformation  transf_1(rot_1, tras_1);

    CaloSurface *face_1= new CaloSurface(transf_1, _vaneHalfThickness, _dZ, _norm);
    caloVane->AddSide(face_1);

    unsigned int idFace_1 = 1;
    caloVane->AddSideCorner( new SurfacePoint(-_vaneHalfThickness, _dZ) , idFace_1);
    caloVane->AddSideCorner( new SurfacePoint(_vaneHalfThickness, _dZ)  , idFace_1);//new SurfacePoint(0.0, 0.0), idFace_0 );
    caloVane->AddSideCorner( new SurfacePoint(_vaneHalfThickness, -_dZ) , idFace_1);
    caloVane->AddSideCorner( new SurfacePoint(-_vaneHalfThickness, -_dZ), idFace_1);

    //face 2
    const Hep3Vector vec_2(0.0, -_vaneHalfThickness , 0.0);//(0.0, -_vaneThickOffSet, 0.0);// = tmp_vec_2;
    const Hep3Vector norm_2(0., 0., 1.);//(0.0, -1.0, 0.0);//-y axis

    const HepRotation rot_2(norm_2,  Constants::pi );
    const HepTranslation tras_2(vec_2);
    const HepTransformation  transf_2(rot_2, tras_2);

    // vector and normal ;
    CaloSurface *face_2= new CaloSurface(transf_2, _dR, _dZ, _norm);
    caloVane->AddSide(face_2);

    unsigned int idFace_2 = 2;
    caloVane->AddSideCorner( new SurfacePoint(-_dR, _dZ) , idFace_2);
    caloVane->AddSideCorner( new SurfacePoint(_dR, _dZ)  , idFace_2);
    caloVane->AddSideCorner( new SurfacePoint(_dR, -_dZ) , idFace_2);
    caloVane->AddSideCorner( new SurfacePoint(-_dR, -_dZ), idFace_2);

    //face 3
    const Hep3Vector vec_3(-_dR, 0.0, 0.0);// = tmp_vec_3;
    const Hep3Vector norm_3(0.0, 0.0, 1.0);//(-1.0, 0.0, 0.0);//-x axis

    const HepRotation rot_3(norm_3, 0.5*Constants::pi);
    const HepTranslation tras_3(vec_3);
    const HepTransformation  transf_3(rot_3, tras_3);

    // vector and normal ;
    CaloSurface *face_3= new CaloSurface(transf_3, _vaneHalfThickness, _dZ, _norm);
    caloVane->AddSide(face_3);

    unsigned int idFace_3 = 3;
    caloVane->AddSideCorner( new SurfacePoint(-_vaneHalfThickness, _dZ) , idFace_3);
    caloVane->AddSideCorner( new SurfacePoint(_vaneHalfThickness, _dZ)  , idFace_3);
    caloVane->AddSideCorner( new SurfacePoint(_vaneHalfThickness, -_dZ) , idFace_3);
    caloVane->AddSideCorner( new SurfacePoint(-_vaneHalfThickness, -_dZ), idFace_3);

    //face 4
    // vector and normal ;
    const Hep3Vector vec_4(0.0, 0.0,-_dZ);// = tmp_vec_4;
    const Hep3Vector norm_4(1.0, 0.0, 0.0);//-z axis

    const HepRotation rot_4(norm_4, -0.5*Constants::pi);
    const HepTranslation tras_4(vec_4);
    const HepTransformation  transf_4(rot_4, tras_4);
    CaloSurface *face_4= new CaloSurface(transf_4, _dR, _vaneHalfThickness, _norm);
    caloVane->AddSide(face_4);

    unsigned int idFace_4 = 4;
    caloVane->AddSideCorner( new SurfacePoint(-_dR, _vaneHalfThickness) , idFace_4);
    caloVane->AddSideCorner( new SurfacePoint(_dR, _vaneHalfThickness)  , idFace_4);
    caloVane->AddSideCorner( new SurfacePoint(_dR, -_vaneHalfThickness) , idFace_4);
    caloVane->AddSideCorner( new SurfacePoint(-_dR, -_vaneHalfThickness), idFace_4);

    //face 5
    const Hep3Vector vec_5(0.0, 0.0, _dZ);// = tmp_vec_5;
    const Hep3Vector norm_5(1.0, 0.0, 0.0);//z axis

    const HepRotation rot_5(norm_5, 0.5*Constants::pi);
    const HepTranslation tras_5(vec_5);
    const HepTransformation  transf_5(rot_5, tras_5);

    CaloSurface *face_5= new CaloSurface(transf_5, _dR, _vaneHalfThickness, _norm);
    caloVane->AddSide(face_5);

    unsigned int idFace_5 = 5;
    caloVane->AddSideCorner( new SurfacePoint(-_dR, _vaneHalfThickness) , idFace_5);
    caloVane->AddSideCorner( new SurfacePoint(_dR, _vaneHalfThickness)  , idFace_5);
    caloVane->AddSideCorner( new SurfacePoint(_dR, -_vaneHalfThickness) , idFace_5);
    caloVane->AddSideCorner( new SurfacePoint(-_dR, -_vaneHalfThickness), idFace_5);

    Hep3Vector tmp_vec = vane.origin();
    double b0 = tmp_vec.getX() + _solenoidOffSetX;
    double c0 = tmp_vec.getZ() + _solenoidOffSetZ;
    tmp_vec.setZ(c0);
    tmp_vec.setX(b0);
    //tmp_vec.setY(0.0);

    const Hep3Vector vector = tmp_vec;
    const Hep3Vector normal(0.0, 0.0, 1.0);


    const HepRotation rot(normal, 0.0);
    const HepTranslation tras(tmp_vec);

    const Hep3Vector axes(vane.rotation().getAxis());

    //double delta = vane.getRotation()->getDelta() - Constants::pi*0.5;
    double delta = Constants::pi*0.5*i;
    const HepRotation rot2(normal, delta );

    const HepTransformation  theAlignment(/*rot*/rot2, tras);


    int caloId = i;
    std::string name= "calo";

    CaloVolumeElem *calo = new CaloVolumeElem(caloVane, name.c_str(), caloId, theAlignment);
    return calo;
  }
  
  CLHEP::Hep3Vector Calorimeter4VanesGeom::fromTrkToMu2eFrame(CLHEP::Hep3Vector  &vec){
    art::ServiceHandle<GeometryService> geom;
    double solenoidOffSetX = geom->config().getDouble("mu2e.solenoidOffset");
    double solenoidOffSetZ = -geom->config().getDouble("mu2e.detectorSystemZ0");
    CLHEP::Hep3Vector res;
    
    res.setX(vec.x() - solenoidOffSetX);
    res.setZ(vec.z() - solenoidOffSetZ);
    res.setY(vec.y());
    return res;
  }
  
  bool Calorimeter4VanesGeom::behindVane(double& posX, double& posY, int& vane){
    if(vane == 3){
      if( fabs( posX ) <= _vaneHalfThickness && posY < _outherRadius && posY > _innerRadius){
	return true;
      }else return false;
    }else if(vane == 2){

      if( fabs( posY ) <= _vaneHalfThickness && posX < _outherRadius && posX > _innerRadius){
	return true;
      }else return false;
    }else if(vane == 1){
      if( fabs( posX ) <= _vaneHalfThickness && posY < -_innerRadius && posY > -_outherRadius){
	return true;
      }else return false;
    }else if(vane == 0){
      if( fabs( posY ) <= _vaneHalfThickness && posX < -_innerRadius && posX > -_outherRadius){
	return true;
      }else return false;
    }else return false;
  }
  Calorimeter4VanesGeom::~Calorimeter4VanesGeom (){}
  //-----------------------------------------------------------------------------
  bool Calorimeter4VanesGeom::behindVane(HepPoint pos, int& vane){

    bool rc(false);

    if (pos.z() >= _ZfrontFaceCalo && (pos.z() <= _ZbackFaceCalo)) {

      if (vane == 3) {
	if (fabs(pos.x()) <= _vaneHalfThickness && pos.y() <= _outherRadius && pos.y() >= _innerRadius) rc = true;
      }
      else if (vane == 2) {
	if (fabs(pos.y()) <= _vaneHalfThickness && pos.x() <= _outherRadius && pos.x() >= _innerRadius) rc = true;
      } 
      else if(vane == 1){
	if (fabs( pos.x() ) <= _vaneHalfThickness && pos.y() <= -_innerRadius && pos.y() >= -_outherRadius) rc = true;
      }
      else if (vane == 0) {
	if (fabs(pos.y() ) <= _vaneHalfThickness && pos.x() <= -_innerRadius && pos.x() >= -_outherRadius) rc = true;
      } 
    }
    return rc;
  }

  //-----------------------------------------------------------------------------
  void Calorimeter4VanesGeom::caloExtrapol(int&             diagLevel,
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
    
    static const char* oname = "Calorimeter4VanesGeom::caloExtrapol";

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

    int nVanes = _nVanes;

    double *entr   = new double[nVanes];
    double *ex     = new double[nVanes];
    bool *isInside = new bool[nVanes];

    for(int jVane=0; jVane<nVanes; ++jVane){
      isInside[jVane] = false;
      entr[jVane] = 0.0;
      ex[jVane] = 0.0;
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
 //    if(geom->hasElement<VaneCalorimeter>()){
//       GeomHandle<VaneCalorimeter> cg;// 05 - 21 - 2013 gianipez
//       for(int iStep = 0; iStep< nAngleSteps; ++iStep){
// 	for(int jVane=0; jVane<nVanes; ++jVane){
// 	  // if(diagLevel>4){
// // 	    cout<<" tmpRange = "<< tmpRange<<
// // 	      ", trj.position(tmpRange) = "<<traj.position(tmpRange)<<endl;
// // 	  }
// 	  trjPoint = traj.position(tmpRange);
// 	  trjVec.setX(trjPoint.x());// 05 - 21 - 2013 gianipez
// 	  trjVec.setY(trjPoint.y());// 05 - 21 - 2013 gianipez
// 	  trjVec.setZ(trjPoint.z());// 05 - 21 - 2013 gianipez

// 	  trjVec = fromTrkToMu2eFrame(trjVec);// 05 - 21 - 2013 gianipez
// 	  if(diagLevel > 4){
// 	    printf("\nrange = %10.3f\n trj.position(%10.3f, %10.3f, %10.3f)\n"
// 		    , tmpRange, trjVec.x(), trjVec.y(), trjVec.z() );
// 	  }
// 	  if( cg->isInsideVane(jVane,trjVec ) ){// if(behindVane(traj.position(tmpRange), jVane) ){
// 	    if(!isInside[jVane]){
// 	      if(diagLevel>4){
// 		cout<<"Event Number : "<< evtNumber<< endl;
// 		cout<<" vane "<<jVane<<
// 		  "isInside : true"<<
// 		  "pathLength entrance = "<<tmpRange<<endl;
// 	      }
// 	      isInside[jVane] = true;
// 	      if(fdir.dzdt() == 1.0){
// 		entr[jVane] = tmpRange - pathStepSize;
// 	      }else if(fdir.dzdt() == -1.0){
// 		entr[jVane] = tmpRange + pathStepSize;
// 	      }
// 	    }
// 	  }else if(isInside[jVane]){
// 	    if(fdir.dzdt() == 1.0){
// 	      ex[jVane] = tmpRange + pathStepSize;
// 	    }else if(fdir.dzdt() == -1.0){
// 	      ex[jVane] = tmpRange - pathStepSize;
// 	    }
	    
// 	    if(diagLevel>4){
// 	      cout<<"Event Number : "<< evtNumber<< endl;
// 	      cout<<" vane "<<jVane<<
// 		"isInside : true"<<
// 		"hasExit : true"<<
// 		"pathLength entrance = "<<entr[jVane]<<
// 		"pathLength exit = "<<tmpRange<<endl;
// 	    }
// 	    isInside[jVane] = false;

// 	    if (NIntersections < 100) {
// 	      Intersection[NIntersections].fVane  = jVane;
// 	      Intersection[NIntersections].fRC    = 0;
// 	      Intersection[NIntersections].fSEntr = entr[jVane];
// 	      Intersection[NIntersections].fSExit = ex  [jVane];
// 	      NIntersections++;
// 	    }
// 	    else {
// 	      printf("%s ERROR: NIntersections > 100, TRUNCATE LIST\n",oname);
// 	    }
// 	  }
// 	}
// 	if(fdir.dzdt() == 1.0){
// 	  tmpRange += pathStepSize;
// 	}else if(fdir.dzdt() == -1.0){
// 	  tmpRange -= pathStepSize;
// 	}
//       }
//     }else 

//    if(geom->hasElement<DiskCalorimeter>()){

    GeomHandle<DiskCalorimeter> cg;
    for(int iStep = 0; iStep< nAngleSteps; ++iStep){
      for(int jVane=0; jVane<nVanes; ++jVane){
	trjPoint = traj.position(tmpRange);
	if(diagLevel>4){
	  cout<<" tmpRange = "<< tmpRange<<
	    ", trj.position(tmpRange) = "<<trjPoint<<endl;
	}
	  
	trjVec.setX(trjPoint.x());
	trjVec.setY(trjPoint.y());
	trjVec.setZ(trjPoint.z());
	
	trjVec = fromTrkToMu2eFrame(trjVec);
	
	if( cg->isInsideDisk(jVane,trjVec ) ){
	  if(!isInside[jVane]){
	    if(diagLevel>4){
	      cout<<"Event Number : "<< evtNumber<< endl;
	      cout<<" vane "<<jVane<<
		"isInside : true"<<
		"pathLength entrance = "<<tmpRange<<endl;
	    }
	    isInside[jVane] = true;
	    if(fdir.dzdt() == 1.0){
	      entr[jVane] = tmpRange - pathStepSize;
	    }else if(fdir.dzdt() == -1.0){
	      entr[jVane] = tmpRange + pathStepSize;
	    }
	  }
	}else if(isInside[jVane]){
	  ex[jVane] = tmpRange + pathStepSize;
	  if(diagLevel>4){
	    cout<<"Event Number : "<< evtNumber<< endl;
	    cout<<" vane "<<jVane<<
	      "isInside : true"<<
	      "hasExit : true"<<
	      "pathLength entrance = "<<entr[jVane]<<
	      "pathLength exit = "<<tmpRange<<endl;
	  }
	  isInside[jVane] = false;
	  
	  if (NIntersections < 100) {
	    Intersection[NIntersections].fVane  = jVane;
	    Intersection[NIntersections].fRC    = 0;
	    Intersection[NIntersections].fSEntr = entr[jVane];
	    Intersection[NIntersections].fSExit = ex  [jVane];
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
      cout<<"end search behindVane(), position is : "<<traj.position(tmpRange)<<endl;
    }

    //    int        resT = -1;
    double     lrange, hrange;
    TrkErrCode trk_rc;

    for (int i=0; i<NIntersections; i++) {
      lrange = Intersection[i].fSEntr;
      //      hrange = Intersection[i].fSEntr;
      trk_rc = Krep->extendThrough(lrange);
      if (trk_rc.success() != 1) { 
	//-----------------------------------------------------------------------------
	// failed to extend 
	//-----------------------------------------------------------------------------
	Intersection[i].fRC = -1;
	if (diagLevel>2) {
	  printf("%s ERROR vane = %2i FAILED to EXTEND TRAJECTORY, rc = %i\n",
		 oname,Intersection[i].fVane,trk_rc.success());
	}
      }
      //       else {
      // 	TrkDifTraj const& traj =  Krep->traj();
	
      // 	//	DetIntersection intersec0( vane(Intersection[i].fVane), &traj, lrange, lrange, hrange );
      // 	// turn off the BaBar extrapolation, it seems to be getting confused
      // 	//	resT = 1; // resT =  vane(jVane)->intersect(&traj, intersec0);
      // 	//	if (resT == 1){
      // 	//                          // update path lengths
      // 	// 	  Intersection[i].fSEntr = intersec0.pathrange[0];
      // 	// 	  Intersection[i].fSExit = intersec0.pathrange[1];
      // 	//      }
      //       }
    }

    delete [] isInside ;
    delete [] entr ;
    delete [] ex;

  }//end proce_dUre


  //the following procedure is different from the previous one for two reasons:
  // -> doesn't have the _diagLevel integer control;
  // -> doesn't have the event index;
  // -> it returns the DetIntersection element relative to the first intersection of the trajectory with the volume,
  //    the previous one returns the DetIntersection element relative to the last intersection.
  void Calorimeter4VanesGeom::caloExtrapol(KalRep* Krep, 
					   double& lowrange, 
					   double& highrange,
					   HelixTraj &trkHel, 
					   int &res0, 
					   DetIntersection &intersec0, 
					   Length *pathLengths)
  {
    //if(Krep->extendThrough(lowrange).success() != 1) return;

    TrkDifTraj const &traj = Krep->traj();

    double circleRadius;
    double startLowrange = lowrange;//, startHighrange = highrange;
    circleRadius = fabs(1.0/trkHel.omega());

    const int nVanes = _nVanes;

    double *entr = new double[nVanes];
    double *ex = new double[nVanes];
    bool *isInside = new bool[nVanes];
    for(int jVane=0; jVane<nVanes; ++jVane){
      isInside[jVane] = false;
      entr[jVane] = 0.0;
      ex[jVane] = 0.0;
    }
    int nAngleSteps = 500;

    double pathStepSize = Constants::twoPi / (double) nAngleSteps;
    nAngleSteps *= 2.0;

    pathStepSize *= circleRadius/fabs(trkHel.cosDip());

    double tmpRange = startLowrange;
    int resT = -1;
    Length tmpPathLengths[nVanes];

    for(int iStep = 0; iStep< nAngleSteps; ++iStep){
      for(int jVane=0; jVane<nVanes; ++jVane){

	if(behindVane(traj.position(tmpRange), jVane) ){
	  if(!isInside[jVane]){

	    isInside[jVane] = true;
	    entr[jVane] = tmpRange - pathStepSize;

	  }
	}else if(isInside[jVane]){
	  ex[jVane] = tmpRange + pathStepSize;
	  isInside[jVane] = false;
	  tmpPathLengths[jVane].push_back(std::pair<double, double>(entr[jVane], ex[jVane]) );
	}

      }
      tmpRange += pathStepSize;
    }

    bool isFirst = false;

    for(int jVane=0; jVane<nVanes; ++jVane){
      isFirst = false;
      for(Length::iterator it = tmpPathLengths[jVane].begin(); it != tmpPathLengths[jVane].end(); ++it){
	lowrange = it->first;
	highrange = it->second;


	//if(Krep->extendThrough(lowrange).success() != 1) continue;

	TrkDifTraj const& traj =  Krep->traj();

	DetIntersection tmp3( vane(jVane), &traj, lowrange, lowrange, highrange );
	//intersec0 = tmp3;
	resT =  vane(jVane)->intersect(&traj, tmp3/*intersec0*/);
	if(resT == 1){

	  pathLengths[jVane].push_back(std::pair<double, double>(intersec0.pathrange[0], intersec0.pathrange[1]) );
	}

	if(resT==1 && !isFirst){
	  isFirst = true;
	  res0 = resT;
	  intersec0 = tmp3;
	}

      }
    }


    delete [] isInside ;
    delete [] entr ;
    delete [] ex;
  }//end proce_dUre


  void Calorimeter4VanesGeom::minimumPathLength(Length *length, int& vane, double& lowrange, double& highrange){
  
    const int nVanes = _nVanes;
    Length::iterator it = length[0].begin();
    lowrange = it->first;
    for(int jVane=0; jVane<nVanes; ++jVane){
      for( it = length[jVane].begin(); it != length[jVane].end(); ++it){
	if(it->first < lowrange){
	  lowrange = it->first;
	  vane = jVane;
	  highrange = it->second;
	}

      }
    }
  }




};


