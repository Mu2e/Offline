///////////////////////////////////////////////////////////////////////////////
//  $Id: 
//  $Author: 
//  $Date: 
//
// cloned from Vadim's code
//
// 2015-07-10 P.Murat: default condiguration is stored in ParticleID/fcl/prolog.fcl
//
// assume that both electron and muon track reconstruction have been attempted,
// so expect two track collections on input
///////////////////////////////////////////////////////////////////////////////

// using namespace std;

#include "ParticleID/inc/AvikPID_module.hh"
#include "TGraphErrors.h"
#include "KalmanTests/inc/KalFitResult.hh"

namespace mu2e {

  TGraphErrors *error;

  float AvikPID::_pathbounds[AvikPID::nbounds] = {0.5,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.};

//-----------------------------------------------------------------------------
// Vadim's fitting: compute sum of squares of residuals
//-----------------------------------------------------------------------------
  void AvikPID::myfcn(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t) {
    //minimisation function computing the sum of squares of residuals
    Int_t np = error->GetN();
    f = 0;
    Double_t *x = error->GetX();
    Double_t *y = error->GetY();
    Double_t *ey = error->GetEY();
    //   Double_t *ey = error->GetEY();
    
    for (Int_t i=0;i<np;i++) {
      Double_t dr = y[i] - (par[0]+par[1]*x[i]);
      f += (dr*dr)/(ey[i]*ey[i]);
    }
  }

//-----------------------------------------------------------------------------
  int AvikPID::findlowhist(float d){

    if (d<=AvikPID::_pathbounds[0]) return 0;

    else if (d>_pathbounds[0] && d<=_pathbounds[1]) return 0;
    else if (d>_pathbounds[1] && d<=_pathbounds[2]) return 1;
    else if (d>_pathbounds[2] && d<=_pathbounds[3]) return 2;
    else if (d>_pathbounds[3] && d<=_pathbounds[4]) return 3;
    else if (d>_pathbounds[4] && d<=_pathbounds[5]) return 4;
    else if (d>_pathbounds[5] && d<=_pathbounds[6]) return 5;
    else if (d>_pathbounds[6] && d<=_pathbounds[7]) return 6;
    else if (d>_pathbounds[7] && d<=_pathbounds[8]) return 7;
    else if (d>_pathbounds[8] && d<=_pathbounds[9]) return 8;
    else if (d>_pathbounds[9] && d<=_pathbounds[10]) return 9;
    else if (d>_pathbounds[10]) return 10;
  
    else { 
      cout<<"Out of  bounds. Should never end up here\n"; 
    }

    return -9999.;
  }


//-----------------------------------------------------------------------------
  AvikPID::AvikPID(fhicl::ParameterSet const& pset):
    _debugLevel(pset.get<int>("debugLevel")),
    _verbosity (pset.get<int>("verbosity" )),
    _diagLevel (pset.get<int>("diagLevel" )),

    _trkPatRecDemModuleLabel(pset.get<string>("trkPatRecDemModuleLabel")),
    _trkPatRecDmmModuleLabel(pset.get<string>("trkPatRecDmmModuleLabel")),

    _eleDedxTemplateFile(pset.get<std::string>("EleDedxTemplateFile")),
    _muoDedxTemplateFile(pset.get<std::string>("MuoDedxTemplateFile")),
    _pidtree(0)
  {
    _processed_events = -1;

    _iname = _fdir.name() + _tpart.name();
    produces<AvikPIDProductCollection>();

    // location-independent files
    ConfigFileLookupPolicy configFile;

    _eleTemplates = configFile(_eleDedxTemplateFile);
    _muoTemplates = configFile(_muoDedxTemplateFile);

    char name[50];

    TFile* eleDedxTemplateFile = TFile::Open(_eleTemplates.c_str());
    for (int i = 0; i < nbounds; i++){
      sprintf(name,"htempe%d",i);
      eleDedxTemplateFile->GetObject(name,_heletemp[i]);
    }

    TFile* muoDedxTemplateFile = TFile::Open(_muoTemplates.c_str());
    for (int i = 0; i < nbounds; i++){
      sprintf(name,"htempm%d",i);
      muoDedxTemplateFile->GetObject(name,_hmuotemp[i]);
    }
//-----------------------------------------------------------------------------
// all electron and muon De/Dx template histograms are supposed to have the 
// same limits and number of bins
//-----------------------------------------------------------------------------
    _templatesnbins   = _heletemp[0]->GetNbinsX();
    _templateslastbin = _heletemp[0]->GetBinLowEdge(_templatesnbins)+_heletemp[0]->GetBinWidth(1);
    _templatesbinsize = _heletemp[0]->GetBinWidth(1);
  }


//-----------------------------------------------------------------------------
  void AvikPID::beginJob(){

    // histograms

    art::ServiceHandle<art::TFileService> tfs;

    if (_diagLevel) {
      _pidtree = tfs->make<TTree>("PID", "PID info");

      _pidtree->Branch("trkid"          , &_trkid            , "trkid/I");
      _pidtree->Branch("p"              , &_trkmom           , "trkmom/D");
      _pidtree->Branch("drdsVadimEle"   , &_drdsVadimEle     , "drdsVadimEle/D");
      _pidtree->Branch("drdsVadimEleErr", &_drdsVadimEleErr  , "drdsVadimEleErr/D");
      _pidtree->Branch("drdsVadimMuo"   , &_drdsVadimMuo     , "drdsVadimMuo/D");
      _pidtree->Branch("drdsVadimMuoErr", &_drdsVadimMuoErr  , "drdsVadimMuoErr/D");
      _pidtree->Branch("logDedxProbEle" , &_logDedxProbEle   , "logDedxProbEle/D");
      _pidtree->Branch("logDedxProbMuo" , &_logDedxProbMuo   , "logDedxProbMuo/D");
    }
  }

//-----------------------------------------------------------------------------
  void AvikPID::beginRun(art::Run & run){
    if (_verbosity>=2) cout << "AvikPID: From beginRun: " << run.id().run() << endl;


  }

//-----------------------------------------------------------------------------
  void AvikPID::beginSubRun(art::SubRun & lblock ) {
    if (_verbosity>=2) cout << "AvikPID: From beginSubRun. " << endl;
  }

//-----------------------------------------------------------------------------
  void AvikPID::endJob(){
    if (_verbosity>=2) cout << "AvikPID: From endJob. " << endl;
  }

//-----------------------------------------------------------------------------
// Avik's weighted residual
//-----------------------------------------------------------------------------
  double AvikPID::weightedResidual(double R) {
    double wr(1.e10);

    double ar = fabs(R);
    if (ar < 0.2) wr = pow(ar, 1.7);
    else          wr = pow(0.2,1.3)*pow(ar,0.4);

    return wr;
  }

//-----------------------------------------------------------------------------
  void AvikPID::doubletMaker(const KalRep* ele_Trk, const KalRep* muo_Trk) {
            
    art::Handle<mu2e::KalRepPtrCollection> eleKrepsHandle;
    art::Handle<mu2e::KalRepPtrCollection> muoKrepsHandle;
    
    const TrkHotList* ele_hot_list = ele_Trk->hotList();
    const TrkHotList* muo_hot_list = muo_Trk->hotList();
//-----------------------------------------------
// ELECTRONS:
//-----------------------------------------------
    int ele_nhits =0;  // total number of hits

    int ele_ncount=0;  // number of special events (i.e. events in doublets (or triplets etc))
    
    for(TrkHotList::hot_iterator it=ele_hot_list->begin(); it<ele_hot_list->end(); it++) ele_nhits+=1;
    
    int   ele_iamb0 = 0;
    float ele_resall[ele_nhits];
    int   ele_resgood[ele_nhits];
    float ele_res[ele_nhits];
    int   ele_secall[ele_nhits];
    int   ele_devall[ele_nhits];
    int   ele_layall[ele_nhits];
    int   ele_Nall[ele_nhits];
    int   ele_strawall[ele_nhits];
    int   ele_straw[ele_nhits];
    int   ele_iamball[ele_nhits];
    int   ele_nmatchall[ele_nhits] = {0};   //number of matches for each hit
    int   ele_nnlet[6] = {ele_nhits};    //array giving number of singlets, doublets, triplets, ...

    float h1_fltlen = ele_Trk->firstHit()->kalHit()->hitOnTrack()->fltLen() - 10;
    float hn_fltlen = ele_Trk->lastHit ()->kalHit()->hitOnTrack()->fltLen() - 10;
    float entlen = std::min(h1_fltlen,hn_fltlen);
    BbrVectorErr momerr = ele_Trk->momentumErr(entlen);
    CLHEP::Hep3Vector fitmom = ele_Trk->momentum(entlen);
    CLHEP::Hep3Vector momdir = fitmom.unit();
    HepVector momvec(3);
    for (int i=0; i<3; i++) momvec[i] = momdir[i];

    int k = 0;
    Hep3Vector pos;
    double     hitres, hiterr;
    
    for (TrkHotList::hot_iterator it=ele_hot_list->begin(); it<ele_hot_list->end(); it++) {
      
      const mu2e::TrkStrawHit* hit = (const mu2e::TrkStrawHit*) &(*it);
      
      mu2e::Straw*   straw = (mu2e::Straw*) &hit->straw();
      
      ele_secall[k]=straw->id().getSector();
      ele_devall[k]=straw->id().getDevice();
      ele_layall[k]=straw->id().getLayer();
      ele_Nall[k]=straw->id().getStraw();
      ele_strawall[k]= straw->index().asInt();
      ele_iamball[k]= hit->ambig();
      hit->hitPosition(pos);
      ele_resgood[k] = hit->resid(hitres,hiterr,1);
      ele_resall[k]=hitres;
      
      k += 1;
    }
//-----------------------------------------------------------------------------
// loop over the electron track hits
//-----------------------------------------------------------------------------
    for(int i=0; i<ele_nhits; i++) {
      
      if(ele_iamball[i]==0) ele_iamb0+=1;
      
      for(int j=0; j<i; j++) {     //this loop checks current hit against all previous hits in list
	
	if((ele_secall[i]==ele_secall[i-(j+1)]) && 
	   (ele_devall[i]==ele_devall[i-(j+1)]) && 
	   (ele_layall[i]!=ele_layall[i-(j+1)]) && 
	   (abs(ele_Nall[i]-ele_Nall[i-(j+1)])<=2) && 
	   (ele_iamball[i]!=ele_iamball[i-(j+1)]) &&
	   (ele_resgood[i]==1) && (ele_resgood[i-(j+1)]==1)) {   //doublet conditions
	  
	  ele_nmatchall[i]+=1;    //number of matches to other hits
	  if((ele_nmatchall[i-(j+1)]==0) && (ele_nmatchall[i]==1)) {  //i.e. has found a doublet
	    ele_res[ele_ncount+1]=ele_resall[i];   // i is current hit
	    ele_straw[ele_ncount+1]=ele_strawall[i];
	    ele_res[ele_ncount]=ele_resall[i-(j+1)];   // i-(j+1) is successful matched (previous) hit
	    ele_straw[ele_ncount]=ele_strawall[i-(j+1)];
	    ele_ncount+=2;
	    ele_nmatchall[i-(j+1)]+=1;    //number of matches to other hits
	    ele_nnlet[1]+=1;
	    ele_nnlet[0]-=2;
	  }

	  else if (ele_nmatchall[i-(j+1)]==0) {
	    ele_res[ele_ncount]=ele_resall[i-(j+1)];   // i-(j+1) is successful matched (previous) hit
	    ele_straw[ele_ncount]=ele_strawall[i-(j+1)];
	    ele_ncount+=1;
	    ele_nnlet[ele_nmatchall[i]]+=1;
	    ele_nnlet[ele_nmatchall[i]-1]-=1;
	    ele_nmatchall[i-(j+1)]+=1;    //number of matches to other hits
	    ele_nnlet[0]-=1;
	  }
	  
	  else {
	    ele_res[ele_ncount]=ele_resall[i];   // i is current hit
	    ele_straw[ele_ncount]=ele_strawall[i];
	    ele_ncount+=1;
	    ele_nnlet[ele_nmatchall[i-(j+1)]+ele_nmatchall[i]]+=1;
	    ele_nnlet[ele_nmatchall[i-(j+1)]]-=1;
	    ele_nnlet[ele_nmatchall[i]-1]-=1;
	    ele_nmatchall[i-(j+1)]+=1;    //number of matches to other hits
	  }	
	}      
      }    
    }


    /*    printf("ELECTRONS \n");
	  for(int i=0; i<ele_ncount; i++) {
	  printf( "res is:  %4.3f %s %i %s %i \n", ele_res[i], "  straw is: ", ele_straw[i], " and match num is: ", ele_nmatch[i]);
	  }
	  
	  for(int i=0; i<4; i++) {
	  printf( "number of %i%s %i \n", i+1, "-lets is: ", ele_nnlet[i]);
	  }
	  
	  printf("ncount is: %i \n", ele_ncount);
	  printf("nhits is: %i \n", ele_nhits);
	  printf("iamb0 is: %i \n \n", ele_iamb0); */
    
    //---------------------------------------------------------
    
    //MUONS:
    
    int muo_nhits=0;  //total number of hits
    int muo_ncount=0;  //number of special events (i.e. events in doublets (or triplets etc))
    
    for(TrkHotList::hot_iterator it=muo_hot_list->begin(); it<muo_hot_list->end(); it++) muo_nhits+=1;
    
    int   muo_iamb0=0;
    float muo_resall   [muo_nhits];
    float muo_res      [muo_nhits];
    int   muo_resgood  [muo_nhits];
    int   muo_secall   [muo_nhits];
    int   muo_devall   [muo_nhits];
    int   muo_layall   [muo_nhits];
    int   muo_Nall     [muo_nhits];
    int   muo_strawall [muo_nhits];
    int   muo_straw    [muo_nhits];
    int   muo_iamball  [muo_nhits];
    int   muo_nmatchall[muo_nhits] = {0};   //number of matches for each hit
    int   muo_nnlet    [6] = {muo_nhits};    //array giving number of singlets, doublets, triplets, ...
    //    float muo_chi2 = muo_Trk->chisq();
    
    k = 0;
    for(TrkHotList::hot_iterator it=muo_hot_list->begin(); it<muo_hot_list->end(); it++) {
      
      const mu2e::TrkStrawHit* hit = (const mu2e::TrkStrawHit*) &(*it);
      
      mu2e::Straw*   straw = (mu2e::Straw*) &hit->straw();
      
      muo_secall[k]=straw->id().getSector();
      muo_devall[k]=straw->id().getDevice();
      muo_layall[k]=straw->id().getLayer();
      muo_Nall[k]=straw->id().getStraw();
      muo_strawall[k]= straw->index().asInt();
      muo_iamball[k]= hit->ambig();
      hit->hitPosition(pos);
      muo_resgood[k] = hit->resid(hitres,hiterr,1);
      muo_resall[k]=hitres;
      
      k += 1;
    }
    
    for(int i=0; i<muo_nhits; i++) {  //loop goes through each hit in the track
      
      if(muo_iamball[i]==0) muo_iamb0+=1;
      
      for(int j=0; j<i; j++){     //this loop checks current hit against all previous hits in list
	
	if((muo_secall[i]==muo_secall[i-(j+1)]) && 
	   (muo_devall[i]==muo_devall[i-(j+1)]) && 
	   (muo_layall[i]!=muo_layall[i-(j+1)]) && 
	   (abs(muo_Nall[i]-muo_Nall[i-(j+1)])<=2) && 
	   (muo_iamball[i]!=muo_iamball[i-(j+1)]) &&
	   (muo_resgood[i]==1) && (muo_resgood[i-(j+1)]==1)) {   //doublet conditions
	
	  muo_nmatchall[i]+=1;    //number of matches to other hits

	  if((muo_nmatchall[i-(j+1)]==0) && (muo_nmatchall[i]==1)) {  //i.e. has found a doublet
	    muo_res[muo_ncount+1]=muo_resall[i];   // i is current hit
	    muo_straw[muo_ncount+1]=muo_strawall[i];
	    muo_res[muo_ncount]=muo_resall[i-(j+1)];   // i-(j+1) is successful matched (previous) hit
	    muo_straw[muo_ncount]=muo_strawall[i-(j+1)];
	    muo_ncount+=2;
	    muo_nmatchall[i-(j+1)]+=1;    //number of matches to other hits
	    muo_nnlet[1]+=1;
	    muo_nnlet[0]-=2;
	  }
	  else if (muo_nmatchall[i-(j+1)]==0) {
	    muo_res[muo_ncount]=muo_resall[i-(j+1)];   // i-(j+1) is successful matched (previous) hit
	    muo_straw[muo_ncount]=muo_strawall[i-(j+1)];
	    muo_ncount+=1;
	    muo_nnlet[muo_nmatchall[i]]+=1;
	    muo_nnlet[muo_nmatchall[i]-1]-=1;
	    muo_nmatchall[i-(j+1)]+=1;    //number of matches to other hits
	    muo_nnlet[0]-=1;
	  }

	  else {
	    muo_res[muo_ncount]=muo_resall[i];   // i is current hit
	    muo_straw[muo_ncount]=muo_strawall[i];
	    muo_ncount+=1;
	    muo_nnlet[muo_nmatchall[i-(j+1)]+muo_nmatchall[i]]+=1;
	    muo_nnlet[muo_nmatchall[i-(j+1)]]-=1;
	    muo_nnlet[muo_nmatchall[i]-1]-=1;
	    muo_nmatchall[i-(j+1)]+=1;    //number of matches to other hits
	  }      
	}    
      }
    }

    //-----------------------------------------------------
    //ANALYSIS:
    
    float res_ele   [ele_ncount]; 
    float res_muo   [muo_ncount];
    float res_ele_sq[ele_ncount]; 
    float res_muo_sq[muo_ncount];
    float resall_ele[ele_nhits] = {0}; 
    float resall_muo[muo_nhits] = {0};
    float res_ele_sq2[ele_nhits] = {0}; 
    float res_muo_sq2[muo_nhits] = {0};
    float res_ele_sum = 0;
    float res_muo_sum = 0;
    float res_ele_sum2 = 0;
    float res_muo_sum2 = 0;
    float ratio = 1.0;
    float ratio2 = 1.0;
    float logratio;
    float logratio2;    
    float logratio3;
    float matchhits = 0.0;
    float matchhits_all = 0.0;

    for(int i=0; i<ele_ncount; i++) {
      res_ele   [i]=0.0;
      res_ele_sq[i]=0.0;
    }

    for(int i=0; i<muo_ncount; i++) {
      res_muo   [i]=0.0;
      res_muo_sq[i]=0.0;
    }

    for (int i=0; i<ele_nhits; i++) {
      resall_ele [i] = 0.;
      res_ele_sq2[i] = 0.;
    }
    
    for (int i=0; i<muo_nhits; i++) {
      resall_muo [i] = 0.;
      res_muo_sq2[i] = 0.;
    }
    
    for(int i=0; i<ele_ncount; i++) { 
      for(int j=0; j<muo_ncount; j++) {
	if(ele_straw[i]==muo_straw[j]) {   // if a hit is part of doublets in both hypotheses...
	  res_ele[i]=ele_res[i];           // enter residuals for both respectively
	  res_muo[j]=muo_res[j];
	}
      }
    }
    
    for(int i=0; i<ele_nhits; i++) { 
      for(int j=0; j<muo_nhits; j++) {
	if(ele_strawall[i]==muo_strawall[j]) {   // if hit exists in both hypotheses...
	  resall_ele[i]=ele_resall[i];           // enter residuals for both respectively
	  resall_muo[j]=muo_resall[j];
	}
      }
    }
    
    for (int i=0; i<ele_ncount; i++) {
      if (abs(res_ele[i])>1.5) {   //if for either hypothesis residuals are too high
	res_ele[i]=0;
	for (int j=0; j<muo_ncount; j++) {
	  if (ele_straw[i]==muo_straw[j]) res_muo[j]=0;   //set them both to zero
	}
      }
    }
    
    for (int i=0; i<muo_ncount; i++) {
      if (abs(res_muo[i])>1.5) {   //if for either hypothesis residuals are too high
	res_muo[i]=0;
	for (int j=0; j<ele_ncount; j++) {
	  if (muo_straw[i]==ele_straw[j]) res_ele[j]=0;   //set them both to zero
	}
      }
    }
    
    for (int i=0; i<ele_nhits; i++) {
      if (abs(resall_ele[i])>1.5) {   //if for either hypothesis residuals are too high
	resall_ele[i]=0;
	for (int j=0; j<muo_nhits; j++) {
	  if (ele_strawall[i]==muo_strawall[j]) resall_muo[j]=0;   //set them both to zero
	}
      }
    }
    
    for (int i=0; i<muo_nhits; i++) {
      if (abs(resall_muo[i])>1.5) {   // if for either hypothesis residuals are too high
	resall_muo[i]=0;
	for (int j=0; j<ele_nhits; j++) {
	  if (muo_strawall[i]==ele_strawall[j]) resall_ele[j]=0;   //set them both to zero
	}
      }
    }
    
    for(int i=0; i<ele_ncount; i++) {
      if (res_ele[i]!=0) matchhits+=1.0;

      res_ele_sq[i] = weightedResidual(res_ele[i]);
      res_ele_sum  += res_ele_sq[i]; 
    }
    
    
    for(int i=0; i<muo_ncount; i++) {
      res_muo_sq[i] = weightedResidual(res_muo[i]);
      res_muo_sum  += res_muo_sq[i];
    }
    
    for(int i=0; i<ele_nhits; i++) {
      if (resall_ele[i]!=0) matchhits_all+=1.0;

      res_ele_sq2[i] = weightedResidual(resall_ele[i]);
      res_ele_sum2  += res_ele_sq2[i];
    }
    
    for(int i=0; i<muo_nhits; i++) {
      res_muo_sq2[i] = weightedResidual(resall_muo[i]);
      res_muo_sum2  += res_muo_sq2[i];
    }
    
    _sumAvikEle  = res_ele_sum;
    _sumAvikMuo  = res_muo_sum;

    _sq2AvikEle  = res_ele_sum2;
    _sq2AvikMuo  = res_muo_sum2;

    _nMatched    = matchhits;
    _nMatchedAll = matchhits_all;

    if ((res_ele_sum!=0) && (res_muo_sum!=0)) {
      ratio = (res_muo_sum)/(res_ele_sum);
    }
    if ((res_ele_sum2!=0) && (res_muo_sum2!=0)) {
      ratio2 = (res_muo_sum2)/(res_ele_sum2);
    }
    
    logratio  = log(ratio);
    logratio2 = log(ratio2);   
    logratio3 = (matchhits/matchhits_all)*log(ratio) + log(ratio2);


    /* printf("ANALYSIS \n");
       
    for(int i=0; i<ele_ncount; i++) {
    printf( "res_ele is:  %4.4f %s %4.4f \n", res_ele[i], "  res_ele_sq is: ", res_ele_sq[i]);
    }
    
    for(int i=0; i<muo_ncount; i++) {
    printf( "res_muo is:  %4.4f %s %4.4f \n", res_muo[i], "  res_muo_sq is: ", res_muo_sq[i]);
    }
    */
    printf( "res_ele_sum is:  %4.4f %s %4.4f \n", res_ele_sum, "  res_muo_sum is: ", res_muo_sum);
    printf( "res_ele_sum2 is:  %4.4f %s %4.4f \n", res_ele_sum2, "  res_muo_sum2 is: ", res_muo_sum2);
    printf("logratio is: %8.4f \n", logratio);
    printf("logratio2 is: %8.4f \n", logratio2);
    printf("logratio3 is: %8.4f \n", logratio3);
    /*    
    printf("t0true is: %8.4f \n", t0true);
    printf("ele_t0 is: %8.4f \n", ele_t0);
    printf("t0est is %8.4f \n", t0est);
    printf("muo_t0 is: %8.4f \n", muo_t0);*/
  }

//-----------------------------------------------------------------------------
  double AvikPID::calculateDedxProb(std::vector<double> gaspaths , 
				    std::vector<double> edeps    , 
				    TH1D**              templates) {

    static const double _minpath = 0.5;
    static const double _maxpath = 10.;

    double thisprob = 1;

    for (unsigned int ipath = 0; ipath < gaspaths.size(); ipath++){
      double thispath = gaspaths.at(ipath);
      double thisedep = edeps.at(ipath);      

      double tmpprob = 0;
      if (thispath > _minpath && thispath<=_maxpath){
	  int lowhist = findlowhist(thispath);

	  PIDUtilities util;
	  TH1D* hinterp = util.th1dmorph(templates[lowhist],
					 templates[lowhist+1],
					 _pathbounds[lowhist],
					 _pathbounds[lowhist+1],
					 thispath,1,0);
//-----------------------------------------------------------------------------
// probability for this edep
//-----------------------------------------------------------------------------
	  int thisedepbin = -999;
	  if (thisedep > _templateslastbin) thisedepbin = _templatesnbins;
	  else                              thisedepbin = int(thisedep/_templatesbinsize)+1;
	  
	  tmpprob=hinterp->GetBinContent(thisedepbin);

	  hinterp->Delete();
	}

      if (tmpprob> 0) thisprob = thisprob * tmpprob;
    }

    return thisprob;
  }
  

//-----------------------------------------------------------------------------
  bool AvikPID::calculateVadimSlope(std::vector<double>  vresd , 
				    std::vector<double>  vflt  , 
				    std::vector<double>  evresd, 
				    std::vector<double>  evflt ,  
				    double               *slope, 
				    double               *eslope) {

    error = new TGraphErrors(vresd.size(),vflt.data(),vresd.data(),evflt.data(),evresd.data());

    TMinuit *gmMinuit = new TMinuit(2); 
    gmMinuit->SetPrintLevel(-1);
    gmMinuit->SetFCN(myfcn);
    const int dim(2);
    const char par_name[dim][20]={"offset","slope"};
    static Double_t step[dim] = {0.001,0.001};
    Double_t sfpar[dim]={0.0,0.005};
    Double_t errsfpar[dim]={0.0,0.0};
    int ierflg = 0;
    for (int ii = 0; ii<dim; ii++) {    
      gmMinuit->mnparm(ii,par_name[ii],sfpar[ii], step[ii], 0,0,ierflg);
    }
    gmMinuit->FixParameter(0);
    gmMinuit->Migrad();
    bool converged = gmMinuit->fCstatu.Contains("CONVERGED");
    if (!converged) 
      {
        cout <<"-----------TOF Linear fit did not converge---------------------------" <<endl;
        return converged;
      }
    for (int i = 0;i<dim;i++) {
      gmMinuit->GetParameter(i,sfpar[i],errsfpar[i]);
    } 

    *slope = sfpar[1];
    *eslope = errsfpar[1];

    delete error;
  
    return converged;
  }

//-----------------------------------------------------------------------------
// calculate parameters of the straight line fit
//-----------------------------------------------------------------------------
  int AvikPID::CalculateSlope(vector<double>& Fltlen  , 
			      vector<double>& Resid   , 
			      double&         Slope   , 
			      double&         SlopeErr) {
    
    double dl, dr, fltSum(0), resSum(0), fltMean, resMean, fltVar(0), resVar(0), fltRes(0), fltDev, resDev, rCoeff;

    int n = Fltlen.size();

    for (int i=0; i<n; i++) {
      fltSum += Fltlen[i];
      resSum += Resid[i];
    }
    
    fltMean = fltSum/n;
    resMean = resSum/n;
    
    for (int i=0; i<n; i++) {
      dl      = Fltlen[i] - fltMean;
      dr      = Resid [i] - resMean;
      fltVar += dl*dl;
      resVar += dr*dr;
      fltRes += dl*dr;
    }
    
    fltDev = sqrt(fltVar/n);              // sigxx
    resDev = sqrt(resVar/n);              // sigyy
    rCoeff = fltRes/sqrt(fltVar*resVar);  // 

    Slope    = rCoeff*resDev/fltDev;
    SlopeErr = 0.1/sqrt(fltVar);

    return 0;
  }


//-----------------------------------------------------------------------------
  int AvikPID::AddHits(const Doublet* Multiplet, vector<double>& Fltlen, vector<double>& Resid) {
    //    int nhits = Multiplet->fNstrawHits;
    int ihit;
    double res, flt, reserr;

    for (int i=0; i<2; i++) {
      ihit = Multiplet->fHitIndex[i];
      flt  = Multiplet->fHit[ihit]->fltLen();
      Multiplet->fHit[ihit]->resid(res, reserr, true);
      res  = (Multiplet->fHit[ihit]->poca()->doca()>0?res:-res);
      Fltlen.push_back(flt);
      Resid.push_back(res);
    }

    return 0;
  }

//-----------------------------------------------------------------------------
  int AvikPID::AddSsMultiplets(const vector<Doublet>* ListOfDoublets,
			       vector<double>&        Fltlen        , 
			       vector<double>&        Resid         ) {

    const mu2e::Doublet  *multiplet, *mj;
    int ndblts    = ListOfDoublets->size();
    int best, bestj;
    double   trkdxdzj, bestdxdzj, dxdzresidj;
      
    for (int i=0; i<ndblts; i++) {
      multiplet       = &ListOfDoublets->at(i);
      int sid    = multiplet->fStationId/2;
      int nhits       = multiplet->fNstrawHits;
      if (nhits > 1) {
	best         = multiplet->fIBest;
	double trkdxdz   = multiplet->fTrkDxDz;
	double bestdxdz  = multiplet->fDxDz[best];
	double dxdzresid = fabs(trkdxdz - bestdxdz);
				// always require the local doublet slope to be close to that of the track
	if (dxdzresid < .1) {
	  if ((best == 1) || (best == 3)) {
	    // OS doublet
	    AddHits(multiplet,Fltlen,Resid);
	  }
	  else {
//-----------------------------------------------------------------------------
// SS doublet: check if there is a OS doublet in the same station and add this 
// double only if the OS one exists
//-----------------------------------------------------------------------------
	    for (int j=0; j<ndblts; j++) {
	      mj            = &ListOfDoublets->at(j);
	      int sidj      = mj->fStationId/2;
	      if (sid == sidj) {
		int nhj     = mj->fNstrawHits;
		if (nhj > 1) {
//-----------------------------------------------------------------------------
// to begin with, don't use single hit multiplets
// however, in principle they could be used
//-----------------------------------------------------------------------------
		  bestj      = mj->fIBest;
		  trkdxdzj   = mj->fTrkDxDz;
		  bestdxdzj  = mj->fDxDz[bestj];
		  dxdzresidj = fabs(trkdxdzj - bestdxdzj);
		
		// always require the local doublet slope to be close to that of the track
	      
		  if ((dxdzresidj < .1) && ((bestj == 1) || (bestj == 3))) {
		    AddHits(multiplet,Fltlen,Resid);
		    break;
		  }
		}
	      }
	    }
	  }
	}
      }
    }

    return 0;
  }


//--------------------------------------------------------------------------------
  void AvikPID::calculateSsSums(const vector<Doublet>* ele_LOD, 
				const vector<Doublet>* muo_LOD) {
    
    //    const mu2e::Doublet  *multiplet, *mj;
    int ele_used      = 0;
    int muo_used      = 0;
    vector<double> ele_fltLen;
    vector<double> muo_fltLen;
    vector<double> ele_resid;
    vector<double> muo_resid;

    _logRatioSs = 0.;

    AddSsMultiplets(ele_LOD,ele_fltLen,ele_resid);
    AddSsMultiplets(muo_LOD,muo_fltLen,muo_resid);

    ele_used = ele_fltLen.size();
    muo_used = muo_fltLen.size();

    if ((ele_used > 1) && (muo_used > 1)) {

      CalculateSlope(ele_fltLen,ele_resid,_drdsSsEle,_drdsSsEleErr);
      CalculateSlope(muo_fltLen,muo_resid,_drdsSsMuo,_drdsSsMuoErr);

      _logRatioSs = log(fabs(_drdsSsMuo/_drdsSsEle));
    }	 

    // printf("ele_resSlope: %10.4e  muo_resSlope: %10.4e \n", ele_slope, muo_slope);
    // printf("ele_used: %3i  muo_used: %3i \n", ele_used, muo_used);
    // printf("logratio: %6.4f \n", ssLogratio);

    _ele_nusedSs = ele_used;
    _muo_nusedSs = muo_used;
  }

//------------------------------------------------------------------------------
// calculate sums over the local doublet residuals
// residuals are weighted, as Avik is trying to de-weight the tails
//-----------------------------------------------------------------------------
  void AvikPID::calculateOsSums(const vector<Doublet>* ele_LOD, 
				const vector<Doublet>* muo_LOD) 
  {
    const Doublet* multiplet;

    _ele_nusedOs  = 0;
    _muo_nusedOs  = 0;
    _ele_resSumOs = 0.;
    _muo_resSumOs = 0.;

    int ele_nnlets    = ele_LOD->size();
    int muo_nnlets    = muo_LOD->size();
    int ele_ndblts    = 0.;
    int muo_ndblts    = 0.;

    if ((ele_nnlets != muo_nnlets) || (ele_nnlets == 0)) {
      printf("We have problems! BE CAREFUL...ele_ndblts: %i %s %i \n", ele_ndblts, "  muo_ndblts: ", muo_ndblts);
    }

    double pow1  = 1.5;
    double bound = 0.0625;
    double pow2  = 0.25;

    for (int i=0; i<ele_nnlets; i++) {
      multiplet       = &ele_LOD->at(i);
      int nhits       = multiplet->fNstrawHits;
      //      ele_dbltsize[i] = nhits;
      if (nhits >= 2) {
	int best    = multiplet->fIBest;
	//	ele_best[i] = best;
	ele_ndblts += 1;
	if ((best == 1) || (best == 3)) {
	  double trkdxdz   = multiplet->fTrkDxDz;
	  double bestdxdz  = multiplet->fDxDz[best];
	  double dxdzresid = trkdxdz - bestdxdz;
	  if (dxdzresid < 0.5) {
	    //	    fHist4.fele_resids->Fill(dxdzresid);
	    double residsq = 0;
	    if (abs(dxdzresid) < bound) residsq = (1/pow(bound, pow1))*pow(abs(dxdzresid), pow1);
	    else residsq = (1/pow(bound, pow1))*pow(bound, (pow1 - pow2))*pow(abs(dxdzresid), pow2);
	    //	  printf("ele_resid: %6.4f %s %6.4f \n", dxdzresid, "   ele_residsq: ", residsq);	  
	    _ele_resSumOs += residsq;
	    _ele_nusedOs  += 1.0;
	  }
	}
      }
    }

    for (int i=0; i<muo_nnlets; i++) {
      multiplet       = &muo_LOD->at(i);
      int nhits       = multiplet->fNstrawHits;
      //      muo_dbltsize[i] = nhits;
      if (nhits >= 2) {
	int best    = multiplet->fIBest;
	//	muo_best[i] = best;
	muo_ndblts += 1;
	if ((best == 1) || (best == 3)) {
	  double trkdxdz   = multiplet->fTrkDxDz;
	  double bestdxdz  = multiplet->fDxDz[best];
	  double dxdzresid = trkdxdz - bestdxdz;
	  if (dxdzresid < 0.5) {
	    // fHist4.fmuo_resids->Fill(dxdzresid);
	    double residsq = 0;
	    if (abs(dxdzresid) < bound) residsq = (1/pow(bound, pow1))*pow(abs(dxdzresid), pow1);
	    else residsq = (1/pow(bound, pow1))*pow(bound, (pow1 - pow2))*pow(abs(dxdzresid), pow2);
	    //	  printf("muo_resid: %6.4f %s %6.4f \n", dxdzresid, "   muo_residsq: ", residsq);	  
	    _muo_resSumOs += residsq;
	    _muo_nusedOs   += 1.0;
	  }
	}
      }
    }

    _logRatioOs = 1.e16;

    if (_muo_resSumOs != 0.) {
      _logRatioOs = log(_muo_resSumOs/_ele_resSumOs*_ele_nusedOs/_muo_nusedOs);
    }

    // printf("ele_resSum: %7.4f %s %7.4f \n", ele_resSum, "   muo_resSum: ", muo_resSum);
    // printf("ele_used: %6.0f %s %6.0f \n", ele_used, "   muo_used: ", muo_used);
    // printf("logratio: %6.4f \n", osLogratio);

  }



//-----------------------------------------------------------------------------
  void AvikPID::produce(art::Event& event) {

    art::Handle<mu2e::KalRepPtrCollection> eleHandle, muoHandle;
    
    double         resid, residerr, aresd, normflt, normresd, firsthitfltlen, lasthitfltlen, entlen;

    int n_ele_trk, n_muo_trk;

    vector<double> vresd , vresd_muo;
    vector<double> vflt  , vflt_muo;
    vector<double> evflt , mvflt;
    vector<double> evresd, mvresd;

    vector<double> gaspaths;
    vector<double> edeps;
    
    const KalRep           *ele_Trk, *muo_Trk;
    const KalFitResult     *ele_kfres, *muo_kfres;
    const vector<Doublet>  *ele_listOfDoublets;
    const vector<Doublet>  *muo_listOfDoublets;
    
    const TrkHotList       *hots;
    mu2e::TrkStrawHit      *hit;
      
    art::Handle<mu2e::KalFitResultCollection> eleKfresHandle;
    art::Handle<mu2e::KalFitResultCollection> muoKfresHandle;
    
    double eprob, muprob;

    _evtid = event.id().event();
    ++_processed_events;
    
    if (_processed_events%100 == 0) {
      if (_verbosity>=1) cout << "AvikPID: processing " << _processed_events << "-th event at evtid=" << _evtid << endl;
    }
    
    if (_verbosity>=2) cout << "AvikPID: processing " << _processed_events << "-th event at evtid=" << _evtid << endl;
    
    unique_ptr<AvikPIDProductCollection> pids(new AvikPIDProductCollection );
    
    art::Selector  ele_selector(art::ProcessNameSelector("")         && 
				art::ModuleLabelSelector(_trkPatRecDemModuleLabel));
    
    art::Selector  muo_selector(art::ProcessNameSelector("")         && 
				art::ModuleLabelSelector(_trkPatRecDmmModuleLabel));
    
    event.get(ele_selector,eleHandle);
    event.get(muo_selector,muoHandle);
    
    if (! eleHandle.isValid()) {
      printf("TAnaDump::printKalRepPtrCollection: no ELE KalRepPtrCollection for module, BAIL OUT\n");
      goto END;
    }

    if (! muoHandle.isValid()) {
      printf("TAnaDump::printKalRepPtrCollection: no MUO KalRepPtrCollection for module, BAIL OUT\n");
      goto END;
    }

    _listOfEleTracks = eleHandle.product();
    _listOfMuoTracks = muoHandle.product();

    n_ele_trk = _listOfEleTracks->size();
    n_muo_trk = _listOfMuoTracks->size();
//-----------------------------------------------------------------------------
// get KalFitResult's with lists of doublets
//-----------------------------------------------------------------------------
    event.get(ele_selector,eleKfresHandle);
    event.get(muo_selector,muoKfresHandle);

    if (! eleKfresHandle.isValid()) {
      printf("TAnaDump::printKalFitResultCollection: no KalFitResultCollection for module, BAIL OUT\n");
      goto END;
    }
    
    if (! muoKfresHandle.isValid()) {
      printf("TAnaDump::printKalFitResultCollection: no KalFitResultCollection for module, BAIL OUT\n");
      goto END;
    }
//-----------------------------------------------------------------------------
// proceed further
//-----------------------------------------------------------------------------
    if ((n_ele_trk != n_muo_trk) || (n_ele_trk == 0)) {
      printf("We are in trouble! BAIL OUT...n_ele_trk: %i %s %i \n", n_ele_trk, "  n_muo_trk: ", n_muo_trk);
      goto END;
    }

    for (int i=0; i<n_ele_trk; i++) {
      _trkid             = i;

      ele_Trk            = _listOfEleTracks->at(i).get();
      muo_Trk            = _listOfMuoTracks->at(i).get();
      ele_kfres          = &eleKfresHandle->at(i);
      muo_kfres          = &muoKfresHandle->at(i);
      ele_listOfDoublets = &ele_kfres->_listOfDoublets;
      muo_listOfDoublets = &muo_kfres->_listOfDoublets;

      firsthitfltlen = ele_Trk->firstHit()->kalHit()->hitOnTrack()->fltLen() - 10;
      lasthitfltlen  = ele_Trk->lastHit()->kalHit()->hitOnTrack()->fltLen() - 10;
      entlen         = std::min(firsthitfltlen,lasthitfltlen);
      _trkmom        = ele_Trk->momentum(entlen).mag();
//-----------------------------------------------------------------------------
// calculate De/Dx for the electron track
//-----------------------------------------------------------------------------
      hots = ele_Trk->hotList();
      for (TrkHotList::hot_iterator ihot=hots->begin(); ihot != hots->end(); ++ihot) {
	hit = (mu2e::TrkStrawHit*) &(*ihot);
	if (hit->isActive()) {
//-----------------------------------------------------------------------------
// 'unbiased' residual, signed with the radius - if the drift radius is greater than 
// the track-to-wire distance, the residual is positive (need to double check the sign!)
// use active hits only
//-----------------------------------------------------------------------------
	  hit->resid(resid,residerr,true);

	  aresd    = (hit->poca()->doca()>0?resid:-resid);
	  normflt  = hit->fltLen() -  ele_Trk->flt0();
	  normresd = aresd/residerr;
	  
	  vresd.push_back(normresd);
	  vflt.push_back(normflt);
	  evresd.push_back(1.);
	  evflt.push_back(0.1);
	  
	  // 2. * here because KalmanFit reports half the path through gas.
	  
	  gaspaths.push_back(2. * hit->gasPath(hit->driftRadius(),hit->trkTraj()->direction( hit->fltLen() )));
	  
	  edeps.push_back(hit->strawHit().energyDep());
	}
      }
      
      calculateVadimSlope(vresd,vflt,evresd,evflt,&_drdsVadimEle,&_drdsVadimEleErr);
//-----------------------------------------------------------------------------
// calculate De/Dx for the electron track
//-----------------------------------------------------------------------------
      hots = muo_Trk->hotList();
      for (TrkHotList::hot_iterator ihot=hots->begin(); ihot != hots->end(); ++ihot) {
	mu2e::TrkStrawHit* hit = (mu2e::TrkStrawHit*) &(*ihot);
	if (hit->isActive()) {
//-----------------------------------------------------------------------------
// 'unbiased' residual, signed with the radius - if the drift radius is greater than 
// the track-to-wire distance, the residual is positive (need to double check the sign!)
// use active hits only
//-----------------------------------------------------------------------------
	  hit->resid(resid,residerr,true);

	  aresd    = (hit->poca()->doca()>0?resid:-resid);
	  normflt  = hit->fltLen() -  ele_Trk->flt0();
	  normresd = aresd/residerr;
	  
	  vresd_muo.push_back(normresd);
	  vflt_muo.push_back(normflt);
	  mvresd.push_back(1.);
	  mvflt.push_back(0.1);
	}
      }
      
      calculateVadimSlope(vresd_muo,vflt_muo,mvresd,mvflt,&_drdsVadimMuo,&_drdsVadimMuoErr);

      eprob  = calculateDedxProb(gaspaths, edeps, _heletemp);
      muprob = calculateDedxProb(gaspaths, edeps, _hmuotemp);
      
      _logDedxProbEle = log(eprob);
      _logDedxProbMuo = log(muprob);
//-----------------------------------------------------------------------------
// calculate Avik's sums
//-----------------------------------------------------------------------------
      doubletMaker(ele_Trk,muo_Trk);
      //      calculateAvikSums();
//-----------------------------------------------------------------------------
// calculate OS slopes
//-----------------------------------------------------------------------------
      calculateOsSums(ele_listOfDoublets,muo_listOfDoublets);
//-----------------------------------------------------------------------------
// calculate SS slopes
//-----------------------------------------------------------------------------
      calculateSsSums(ele_listOfDoublets,muo_listOfDoublets);
//-----------------------------------------------------------------------------
// form PID object : FIXME
//-----------------------------------------------------------------------------
      _pid.init(_trkid,
		_logDedxProbEle, _logDedxProbMuo  ,
		_drdsVadimEle  , _drdsVadimEleErr ,
		_drdsVadimMuo  , _drdsVadimMuoErr ,
		_nMatched      , _nMatchedAll     ,
		_sumAvikEle    , _sumAvikMuo      ,
		_sq2AvikEle    , _sq2AvikMuo      ,
		_drdsOsEle     , _drdsOsEleErr    ,
		_drdsOsMuo     , _drdsOsMuoErr    ,
		_ele_nusedSs   , _muo_nusedSs     ,
		_drdsSsEle     , _drdsSsEleErr    ,
		_drdsSsMuo     , _drdsSsMuoErr    ,
		_ele_nusedOs   , _muo_nusedOs     ,
		_ele_resSumOs  , _muo_resSumOs
		);

      pids->push_back(_pid);
//-----------------------------------------------------------------------------
// fill ntuple
//-----------------------------------------------------------------------------
      if (_diagLevel) {
	_pidtree->Fill();
      }
    }

  END: ;
    event.put(std::move(pids));
  }

} // end namespace mu2e

using mu2e::AvikPID;
DEFINE_ART_MODULE(AvikPID);
