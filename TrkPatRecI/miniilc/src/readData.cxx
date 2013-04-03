#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"
#include "readData.h"
#include "StrawHit.hh"
#include "IlcDCHcluster.h"
#include "StrawHitMCTruth.hh"
#include <iostream>
#include <vector>

// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
//geometry
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "ITrackerGeom/inc/ITracker.hh"

using namespace std; 

TObjArray *mu2eHits2ilc(const vector<mu2e::StrawHit> *shits,const vector<mu2e::StrawHitMCTruth> *shitmcs,const vector<int>* select,double t0, bool useZCoordinate){
  TObjArray* array=new TObjArray;
  int nhits=select?select->size():shits->size();
  for(int ii=0;ii<nhits;ii++){
    int i=select?(*select)[ii]:ii;
    //std::cout<<(*shits)[i].time()<<" "<<(*shitmcs)[i].driftDistance()<<endl;
    int lab[3]={0,0,0};
    unsigned int id=(*shits)[i].strawIndex().asUint();
    //    cout<<id<<" "<<(id/10000)<<" "<<id%10000<<endl;
    float hit[5]={0,(float)((*shitmcs)[i].distanceToMid()*0.1+gRandom->Gaus()*1),
		  0.005*0.005,1.*1.,(float)(((*shitmcs)[i].driftDistance()*0.1)+gRandom->Gaus()*0.02)};
    double sigimp2=0.02*0.02;
    if(t0>0){
      hit[4]=((*shits)[i].time()-t0)*0.035/10.;
      //      std::cout<<i<<" hits "<<((*shitmcs)[i].driftDistance()*0.1)<<" "<<hit[4]<<endl;
      sigimp2+=0.023*0.023;
    }
    //if(hit[1]>-50) continue;
    IlcDCHcluster *dhit=new IlcDCHcluster(lab,hit);
    dhit->SetId(int(id/10000)*100000+(id%10000));
    dhit->SetSigmaImP2(sigimp2);

    if (useZCoordinate) {
            mu2e::ConditionsHandle<mu2e::TrackerCalibrations> tcal("ignored");
            const mu2e::Tracker& tracker = mu2e::getTrackerOrThrow();
	    mu2e::Straw const& straw= tracker.getStraw((*shits)[i].strawIndex());
	    mu2e::SHInfo shinfo;
            tcal->StrawHitInfo(straw,(*shits)[i],shinfo);
            dhit->SetZ(shinfo._pos.z()*0.1); //in cm
            const mu2e::ITracker &itracker = static_cast<const mu2e::ITracker&>( tracker );
            mu2e::CellGeometryHandle *_itwp = itracker.getCellGeometryHandle();
            _itwp->SelectCellDet(id);
            double z_err = shinfo._tdres*cos(_itwp->GetWireEpsilon())*0.1; //in cm
            dhit->SetSigmaZ2(z_err*z_err);
    }

    dhit->SetQ(1.);
    dhit->SetIndex(i);
    array->AddLast(dhit);
  }
  return array;
}

void readData(){
  vector<mu2e::StrawHit> *shits=new vector<mu2e::StrawHit>;
  vector<mu2e::StrawHitMCTruth> *shitmcs=new vector<mu2e::StrawHitMCTruth>;

  TFile *_file0 = TFile::Open("/data/mu2e/fnal/Sep17/data_it_v42_lght.root");
  TTree *tr=(TTree*)_file0->Get("Events");
  /*tr->SetBranchStatus("*",0);
  tr->SetBranchStatus("mu2e::StrawHits_makeDcH__KALMANFit.*",1);
  tr->SetBranchStatus("mu2e::StrawHitMCTruths_makeDcH__KALMANFit.*",1);*/
  TBranch *brnch=new TBranch;
  TBranch *brnch2=new TBranch;
  tr->SetBranchAddress("mu2e::StrawHits_makeDcH__KALMANFit.obj",shits,&brnch);
  tr->SetBranchAddress("mu2e::StrawHitMCTruths_makeDcH__KALMANFit.obj",shitmcs,&brnch2);
  for(int i=0;i<10;i++){
    int rsize=brnch->GetEvent(i);
    rsize+=brnch2->GetEvent(i);
    std::cout<<shits<<" "<<shitmcs
	     <<" "<<shits->size()<<" "<<shitmcs->size()<<" "<<rsize<<std::endl;
    for(size_t j=0;j<shits->size();j++){
      std::cout<<(*shits)[j].time()<<" "<<(*shitmcs)[j].driftDistance()<<endl;
    }
  }

}
