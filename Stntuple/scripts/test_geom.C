int test_geom_1() {

  TTrajectoryPoint point[2];

  double x[8];
  int    side(-1),wedge(-1);
  double u(1.e6),v(1.e6);


  TSimpleExtrapolator* ex = new TSimpleExtrapolator();

  double px = 0.2;
  double py = 0;
  double pz = 0.4;

  double pmom = sqrt(px*px+py*py+pz*pz);

  
  x[0] = 0;
  x[1] = 0;
  x[2] = 0;
  x[3] = px/pmom;
  x[4] = py/pmom;
  x[5] = pz/pmom;
  x[6] = 0;
  x[7] = pmom;

  point[0].SetPoint(x);

  point[0].Print("banner");

  ex->SwimToPes(point,1,side,wedge,u,v);
  point[1].Print("data");

  ex->SwimToPes(point,-1,side,wedge,u,v);
  point[1].Print("data");


  x[5] = -x[5];
  point[0].SetPoint(x);

  ex->SwimToPes(point,1,side,wedge,u,v);
  point[1].Print("data");

  ex->SwimToPes(point,-1,side,wedge,u,v);
  point[1].Print("data");

}

//_____________________________________________________________________________
int test_bob() {

  SimpleExtrapolatedTrack bob;


  TVector   par(5);
  TVector3 x1;
  TVector3 n1;

  double px = 0.2;
  double py = 0;
  double pz = 0.4;

  double pt = sqrt(px*px+py*py);

  double pmom = sqrt(px*px+py*py+pz*pz);

  //  Double_t pt =  0.5/C0*0.0029979*BField;

  par[0] = pz/pt;
  par[1] = 0.5*0.0029979*1.4116/pt ; // c0
  par[2] = 0 ; // z0
  par[3] = 0 ; // d0
  par[3] = asin(py/pt) ; // phi0 

  bob.loadTrack(par);
  bob.extrapolateZ(185.4);

  x1 = bob.currentSpacePoint();
  n1 = bob.currentMomentum().Unit();

  printf(" %11.4f %11.4f %11.4f %11.4f %11.4f %11.4f\n",
	 x1.X(),x1.Y(),x1.Z(),n1.X(),n1.Y(),n1.Z());

  par[1] = -par[1];
  bob.loadTrack(par);
  bob.extrapolateZ(185.4);


  x1 = bob.currentSpacePoint();
  n1 = bob.currentMomentum().Unit();

  printf(" %11.4f %11.4f %11.4f %11.4f %11.4f %11.4f\n",
	 x1.X(),x1.Y(),x1.Z(),n1.X(),n1.Y(),n1.Z());



  par[1] = -par[1];

  bob.loadTrack(par);
  bob.extrapolateZ(-185.4);


  x1 = bob.currentSpacePoint();
  n1 = bob.currentMomentum().Unit();

  printf(" %11.4f %11.4f %11.4f %11.4f %11.4f %11.4f\n",
	 x1.X(),x1.Y(),x1.Z(),n1.X(),n1.Y(),n1.Z());

  
  par[1] = -par[1];

  bob.loadTrack(par);
  bob.extrapolateZ(-185.4);


  x1 = bob.currentSpacePoint();
  n1 = bob.currentMomentum().Unit();

  printf(" %11.4f  %11.4f  %11.4f  %11.4f  %11.4f  %11.4f \n",
	 x1.X(),x1.Y(),x1.Z(),n1.X(),n1.Y(),n1.Z());

}
