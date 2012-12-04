#include <TGeoMatrix.h> 
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "GeometryService/inc/GeomHandle.hh"
#include "ITrackerGeom/inc/ITracker.hh"
#include "TrkPatRecI/inc/fillWires.hh"

using namespace std; 
using namespace mu2e;

IlcDCHwireposition * fillWires() { 

  art::ServiceHandle<GeometryService> geom;
  GeomHandle<ITracker> itr;
  CellGeometryHandle *itwp = itr->getCellGeometryHandle();
  int nlayers=0;
  vector<int> num_wires;
  vector<double> epsilon1s,rwires,phi0s,dphis,zs;

  int _nLayers=itr->nSuperLayers();
  for(int i=0;i<_nLayers;i++){
    boost::shared_ptr<ITLayer> ily;
    for (int il=0; il<itr->getSuperLayer(i)->nLayers(); il++) {
      ily = itr->getSuperLayer(i)->getLayer(il);
      if (ily->nCells()!=0){
	itwp->SelectCell(i,ily->getCell(0)->Id().getLayer(),0);
	Float_t xup[3],xdown[3];
	itwp->WirePosAtEndcap(xup,xdown);
	double rwire=itwp->GetWireCenter().rho();
	double phi0=itwp->GetWireCenter().phi();
	double dphi=TMath::TwoPi()/ily->nCells();
	double sangle=-itwp->GetWireEpsilon();
	double z=fabs(xdown[2]);
	double num_wire=ily->nCells();
	
	rwires.push_back(rwire/10.);
	phi0s.push_back(phi0);
	dphis.push_back(dphi);
	epsilon1s.push_back(sangle);
	zs.push_back(z/10.);
	num_wires.push_back(num_wire);
	nlayers++;
	break;
      }
    }
  }

  std::cout<<"Loaded nlayers "<<nlayers<<std::endl;

  Double_t gra2deg = TMath::Pi()/180.; 
  TGeoRotation *pMatrix_SenseII = new TGeoRotation("rot_SenseII");
  TGeoRotation *pMatrix_epsilon1 = new TGeoRotation("epsilon1");
  IlcDCHwireposition *storeWireData = new IlcDCHwireposition(-1);

  int nsuperlayer=nlayers;

  for (Int_t superlayer=0;superlayer<nsuperlayer;superlayer++) {
    Int_t ncel=1;
    int num_wire=num_wires[superlayer];
    storeWireData->AddSuperLayer( ncel, num_wire);
    for(Int_t iring=0; iring<ncel ; iring++){
      double epsilon1=epsilon1s[superlayer],rwire=rwires[superlayer],
	phi0=phi0s[superlayer],dphi=dphis[superlayer],z=zs[superlayer];
      double alfa1=TMath::ATan(fabs(z*TMath::Tan(epsilon1)/rwire));
      pMatrix_epsilon1->SetAngles(0.0,-epsilon1/gra2deg,0.0); 
      TGeoCombiTrans *pMat;
      TGeoHMatrix *comb;
      pMat = new TGeoCombiTrans("id_comb_Sense",rwire,0,0,0);
      comb= new TGeoHMatrix((*pMat));
      comb->SetName("rot");
      comb->RegisterYourself();
      storeWireData->InsertAlfaEpsilon(iring, alfa1, epsilon1);
      storeWireData->InsertRadius(iring,rwire);
      for (Int_t i = 0; i< num_wire; i++){ //num_wire
	pMatrix_SenseII->SetAngles((phi0+i*dphi)/gra2deg,0.0,0.0);
	storeWireData->InsertWireMatrix(iring,i,new TGeoHMatrix((*pMatrix_SenseII)*(*comb)*(*pMatrix_epsilon1)));
      }
    }
  
    storeWireData->FillData(superlayer);
  }
  return storeWireData;
}
