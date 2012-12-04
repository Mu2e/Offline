#include <TGeoMatrix.h> 
#include "IlcDCHwireposition.h" 
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
using namespace std; 

void createWires() { 


  ifstream in("wires.txt");
  std::string tmp;
  int nlayers;
  vector<int> num_wires;
  vector<double> epsilon1s,rwires,phi0s,dphis,zs;

  double rwire,phi0,dphi,sangle,z;
  int num_wire;
  while(in>>nlayers>>tmp>>rwire>>tmp>>phi0>>tmp>>dphi>>tmp>>sangle>>tmp>>z>>
	tmp>>num_wire){
    getline(in,tmp);
    num_wires.push_back(num_wire);
    epsilon1s.push_back(sangle);
    rwires.push_back(rwire/10.);
    phi0s.push_back(phi0);
    dphis.push_back(dphi);
    zs.push_back(z/10.);
  }
  cout<<"Number of layers "<<nlayers+1<<" = "<<num_wires.size()<<endl;
       
  Double_t gra2deg = TMath::Pi()/180.; 
  TGeoRotation *pMatrix_SenseII = new TGeoRotation("rot_SenseII");
  TGeoRotation *pMatrix_epsilon1 = new TGeoRotation("epsilon1");
  IlcDCHwireposition *storeWireData = new IlcDCHwireposition(kTRUE);

  int nsuperlayer=nlayers+1;

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
  
    storeWireData->FillData();
  }
  
  storeWireData->WriteData();
  delete storeWireData;

}
