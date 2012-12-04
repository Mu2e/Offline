/**************************************************************************
 * Copyright(c) 2005-2006, ILC Project Experiment, All rights reserved.   *
 *                                                                        *
// Author: The ILC Off-line Project. 
 // Part of the code has been developed by Alice Off-line Project. 
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
 
/* $Id: IlcDCHgeometry.cxx,v 1.1 2012/12/04 00:51:27 tassiell Exp $ */ 
 
/////////////////////////////////////////////////////////////////////////////// 
//                                                                           // 
//  DCH geometry class                                                       // 
//                                                                           // 
/////////////////////////////////////////////////////////////////////////////// 

 
 
#include <TError.h> 

#include <TGeoPhysicalNode.h> 
#include <TGeoMatrix.h> 
#include <TGeoManager.h>  
#include <TGeoShape.h>
#include <TGeoTube.h>
#include <TGeoHype.h>
#include <TGeoSphere.h>
#include <TGeoArb8.h>
#include <TGeoCompositeShape.h>

#include <iostream>
using namespace std; 
// #include "IlcRunLoader.h" 
#include "IlcDCHgeometry.h" 
#include "IlcDCHwireposition.h" 
#include "IlcLog.h"
//#include "IlcDCHpadPlane.h" 
 
//#include "IlcAlignObj.h" 
//#include "IlcAlignObjAngles.h" 
 
// #include "IlcRun.h" 
// #include "IlcDCH.h" 
//#include "IlcDCHcalibDB.h" 
//#include "IlcDCHCommonParam.h" 
//#include "IlcMAG.h"

 
//ClassImp(IlcDCHgeometry) 
 
//_____________________________________________________________________________ 
 
  // 
  // The geometry constants 
  // 
/*   const Int_t   IlcDCHgeometry::fgkNsect   = kNsect; 
  const Int_t   IlcDCHgeometry::fgkNplan   = kNplan; 
  const Int_t   IlcDCHgeometry::fgkNcham   = kNcham; 
  const Int_t   IlcDCHgeometry::fgkNdet    = kNdet; 
 
  // 
  // Dimensions of the detector 
  // 
 
  // Inner and outer radius of the mother volumes  
//const Float_t IlcDCHgeometry::fgkRmin    = 350.0; 
const Float_t IlcDCHgeometry::fgkRmin    = 400.0; 
 
const Float_t IlcDCHgeometry::fgkRmax    = 550.0; 
const Float_t IlcDCHgeometry::fgkCrack    = 0.09; // percentage of cracks 
 
  // Upper and lower length of the mother volumes  
const Float_t IlcDCHgeometry::fgkZmax1   = 600.0;  
const Float_t IlcDCHgeometry::fgkZmax2   = 600.1;  
 
const Float_t IlcDCHgeometry::fgkTradius = 2.3;  
const Float_t  IlcDCHgeometry::fgkstep=4.6; 
 
 // Parameter of the BTR2 and Box1  mother volumes  
  const Float_t IlcDCHgeometry::fgkSheight[2] = { (fgkRmax - fgkRmin)* TMath::Cos(30.*TMath::Pi()/180.), fgkRmax - fgkRmin};   //Delta R 
 
const Float_t IlcDCHgeometry::fgkSwidth1 = 2*fgkRmin/TMath::Cos(GetAlpha2())*TMath::Sin(GetAlpha1());  //Deltax min 
 const Float_t IlcDCHgeometry::fgkSwidth2 = fgkSwidth1+2*fgkSheight[0]*TMath::Tan(30./180.*TMath::Pi());   //Deltax max 
const Float_t IlcDCHgeometry::fgkSlenTR1 = fgkZmax1*2;//Deltaz for DCH  
const Float_t IlcDCHgeometry::fgkSlenTR2 = 2*186.0; //Dz for each wheel 
 // Trapezius Parameters 
 const Float_t IlcDCHgeometry::fgkTRwidth1 = fgkSwidth1;  //Deltax inner Trapezius 
 const Float_t IlcDCHgeometry::fgkTRwidth2 = fgkSwidth2;   //Deltax outer 
 const Float_t IlcDCHgeometry::fgkTRheight = fgkSheight[0];   //Deltar  
  
// Box Parameters 
const Float_t IlcDCHgeometry::fgkBwidth1 =2*fgkRmin*TMath::Tan(GetAlpha2());     //Deltax Box1 
const Float_t IlcDCHgeometry::fgkBSheight = fgkSheight[1];//height Super Box 
const Float_t IlcDCHgeometry::fgkBheight = fgkBSheight/4;//height each Box1 
 
 // The TUBE (Al)  // Amplification region 
  const Float_t IlcDCHgeometry::fgkCamRmin=0.95*fgkTradius; 
  const Float_t IlcDCHgeometry::fgkCamRmax=fgkTradius; 
  
  const Float_t IlcDCHgeometry::fgkCamH    =  2*fgkCamRmin; 
  // Total height for each chamber 
  const Float_t IlcDCHgeometry::fgkCH[2]      = {fgkTradius*2*kNdets,fgkTradius*2*kNdets+2*fgkTradius}; //Tube box height for each plane 
 
  // Vertical spacing of the chambers 
  const Float_t IlcDCHgeometry::fgkVspace[2] = {(fgkSheight[0]-fgkCH[0]*kNplan)/2., (fgkSheight[1]-fgkCH[1]*kNplan)/2.}; 
  // Horizontal spacing of the chambers 
  const Float_t IlcDCHgeometry::fgkHspace  =   2.0; 
 
const Int_t  IlcDCHgeometry::fgkNtub[2] ={int(fgkTRwidth2/fgkTradius/2),int(fgkBwidth1/fgkTradius/2)}; // //Tubes number 4 plane (tr,box) 

    // Thicknesses of different parts of the chamber frame 
  // Lower aluminum frame 
  const Float_t IlcDCHgeometry::fgkCalT    =   0.3; 
 
  // Additional width of the readout wires 
  const Float_t IlcDCHgeometry::fgkCroW    =   0.005; 
 
  // 
  // Thickness of the the material layers 
  // 
  const Float_t IlcDCHgeometry::fgkAmThick = IlcDCHgeometry::fgkCamH; 
 
  // 
  // Position of the material layers 
  // 
  const Float_t IlcDCHgeometry::fgkAmZpos  =  0.0; 
   
  const Double_t IlcDCHgeometry::fgkTime0Base = Rmin() + CamHght()/2.; 
  const Float_t  IlcDCHgeometry::fgkTime0[kNplan]  = { fgkTime0Base + 0 * (Cheight() + Cspace()),  
                                                  fgkTime0Base + 1 * (Cheight() + Cspace()),  
                                                  fgkTime0Base + 2 * Cheight(),  
                                                  fgkTime0Base + 3 * (Cheight() + Cspace())}; 
 */ 
//_____________________________________________________________________________ 
IlcDCHgeometry::IlcDCHgeometry()
{ 
  // 
  // IlcDCHgeometry default constructor 
  // 
 
  fMatrixArray           = 0; 
  fMatrixCorrectionArray = 0;
  fDCHParam              = 0;
  
//  Init(); 
 
} 
 
//_____________________________________________________________________________ 
IlcDCHgeometry::~IlcDCHgeometry() 
{ 
  // 
  // IlcDCHgeometry destructor 
  // 
 
  delete fMatrixArray; 
  delete fMatrixCorrectionArray; 
 
} 
 


//_____________________________________________________________________________ 
 
 void IlcDCHgeometry::Init() 

{ 
 Double_t gra2deg = TMath::Pi()/180.;
 if ( !fDCHParam ) {
   IlcError("No Parameters found!!!!");
   IlcError("Geometry construction aborted!!!!");
   return;
 }
 
  Double_t inner_radius=fDCHParam->GetInnerRadius(); //142.5549;
//  Double_t outer_radius=fDCHParam->GetOuterRadius(); //148.0343;
  Double_t fieldwire_diameter=fDCHParam->GetFWireDiameter();//.01;
  Double_t envelop_Inner_thickness=fDCHParam->GetInnerWallThickness();//0.02;
//  Double_t envelop_Outer_thickness=fDCHParam->GetOuterWallThickness(); //.;
//  Double_t envelop_EndCap_thickness=fDCHParam->GetEndCapWallThickness();//2.;
  Int_t num_wire_sense=fDCHParam->GetSWireNum(); //80;
  Int_t delta_num_wire_sense=fDCHParam->GetSDeltaWireNum(); //20;
  Int_t nsuperlayer=fDCHParam->GetSuperLayerNum(); //20;
  Int_t nring=fDCHParam->GetRingNum(); //11;
  Double_t drop=fDCHParam->GetDrop(); //2.;
  Double_t length=fDCHParam->GetLength(); //0.5*300.;
//  Double_t extra_EndCap_dist = fDCHParam->GetExtraEndCapDist();
//  Int_t EndCap_type=fDCHParam->GetEndCapType();
//  Double_t max_EndCap_dim=fDCHParam->GetMaxEndCapDim(); //212. used if EndCap_type!=0
//  if(EndCap_type==0) {
//    length = length-envelop_EndCap_thickness;
//    extra_EndCap_dist=0.;
//  }
//  else if(EndCap_type==1){
//    length=(max_EndCap_dim-envelop_EndCap_thickness)*length/max_EndCap_dim; ///TMath::Sqrt(2.);
//  /*  extra_EndCap_dist=TMath::Sqrt(pow(max_EndCap_dim-envelop_EndCap_thickness,2)-pow(inner_radius,2))-length;
//    EndCap_Wall_theta_outer=TMath::ACos(length/max_EndCap_dim)/gra2deg;
//    EndCap_Wall_theta_inner=TMath::ACos((length+extra_EndCap_dist)/(max_EndCap_dim-envelop_EndCap_thickness))/gra2deg;
//*/
//  }
  Double_t radius_ring, radius_ring_0, radius_ring2, radius_ring2_0, theta_ring, delta_radius_ring, ringangle;
  Double_t alfa1, alfa2;
  Double_t epsilon1, epsilon2;
//  Double_t thz_Sense_0,thz_Sense_2_0;
  Double_t thz_Sense,z= 0.000000;
  Int_t num_wire;
  Double_t kbase_exag=0;
  radius_ring_0 = inner_radius + envelop_Inner_thickness + 0.5*fieldwire_diameter;
  radius_ring = radius_ring_0 + drop;
  radius_ring2 = radius_ring;
  radius_ring2_0 = radius_ring_0;
  TGeoRotation *pMatrix_SenseII = new TGeoRotation("rot_SenseII");
  TGeoRotation *pMatrix_epsilon1 = new TGeoRotation("epsilon1");
  TGeoRotation *pMatrix_epsilon2 = new TGeoRotation("epsilon2");
  Int_t sign_epsilon=-1;
//  Double_t wre[24];
  Double_t radiuss=0;
//  Double_t Drd=0.000001;
  //Double_t lenwire_endcap=0.;
//  Double_t wri=0.;
//  Double_t gri=0.;
//  Double_t gre=0.;
  Double_t rwire=0.;
  Double_t phi=0.;
  IlcDCHwireposition *storeWireData = new IlcDCHwireposition(kTRUE);
  Int_t ncel;
  for (Int_t superlayer=0;superlayer<nsuperlayer;superlayer++) {
    // count_node = 0;
//     count_subsector = 0;
    num_wire=num_wire_sense+superlayer*delta_num_wire_sense;

    if(superlayer>=10 && superlayer<18)  {num_wire=num_wire-4*delta_num_wire_sense;
    }
    if(superlayer>=18 && superlayer<24)  num_wire=num_wire-6*delta_num_wire_sense;

    theta_ring=2.*(Double_t) TMath::Pi()/(3*num_wire);
    phi=2.*(Double_t) TMath::Pi()/(num_wire);
    sign_epsilon*=-1;
    ncel = nring-1;
    storeWireData->AddSuperLayer( ncel, num_wire);
    for(Int_t iring=0; iring< nring+1 ; iring++){
      if ((iring%2)==0){
  	ringangle = 0.;
  	
      }
      else{
  	ringangle = -(1.5*theta_ring);
  	
      }
      kbase_exag = 2.*radius_ring_0*TMath::Sin(theta_ring*0.5);
      fBaseExag[superlayer] = kbase_exag;
      delta_radius_ring = kbase_exag * TMath::Cos(30. * gra2deg);
      if(iring!=nring){
	radius_ring2 += delta_radius_ring;
	radius_ring2_0 += delta_radius_ring;
      }
      alfa1 = TMath::ACos(1.-(drop/radius_ring));
      alfa2 = TMath::ACos(1.-(drop/radius_ring2));
      epsilon1 = TMath::ATan(radius_ring/length*TMath::Sin(alfa1));
      epsilon2 = TMath::ATan(radius_ring2/length*TMath::Sin(alfa2));
      pMatrix_epsilon1->SetAngles(0.0,-sign_epsilon*epsilon1/gra2deg,0.0); 
      pMatrix_epsilon2->SetAngles(0.0,-sign_epsilon*epsilon2/gra2deg,0.0); 
//      thz_Sense_0 = ringangle;
//      thz_Sense_2_0 = ((theta_ring*1.5) + ringangle);
      if(iring==0){
	radiuss=radius_ring_0;
//	wri=radiuss;
//	wre[iring]=radiuss+2*fieldwire_diameter;
	rwire=radiuss+fieldwire_diameter;
      }else{ 
	radiuss=radius_ring_0;
	if(iring!=nring){
//	  wri=radiuss-fieldwire_diameter;
//	  wre[iring]=radiuss+fieldwire_diameter;
	  rwire=radiuss;
	}else{
//	  wri=radiuss-2*fieldwire_diameter;
//	  wre[iring]=radiuss;
	  rwire=radiuss-fieldwire_diameter;
	}
//	gri=wre[iring-1]+Drd;
//	gre=wri-Drd;
      }
//	if(iring==nring && superlayer!=nsuperlayer-1){
//	  gri=wre[iring]+Drd;
//	  gre=gri+fieldwire_diameter;
//	}
//	if(superlayer==nsuperlayer-1 && iring==nring ) {
//	  gri=wre[iring]+Drd;//+2*delta_radius_ring;
//	  gre=outer_radius-envelop_Outer_thickness;
//	}
//	if(superlayer==0 && iring==0) {
//	  gri=inner_radius + envelop_Inner_thickness;
//	  gre=radius_ring_0;
//	}
	TGeoCombiTrans *pMat;
	TGeoHMatrix *comb;
	phi=ringangle;
	Double_t x1_Sense_2 = rwire;
	Double_t y1_Sense_2 = 0;
	pMat = new TGeoCombiTrans("id_comb_Sense",x1_Sense_2,y1_Sense_2,z,0);
	comb= new TGeoHMatrix((*pMat));
	comb->SetName("rot");
	comb->RegisterYourself();
	if (iring != 0 && iring!=nring) {
	  //	  storeWireData->InsertWireMatrix(iring-1,0,comb);
	  storeWireData->InsertAlfaEpsilon(iring-1, alfa1, sign_epsilon*epsilon1);
	  storeWireData->InsertRadius(iring-1, radius_ring_0);
	}
	for (Int_t i = 0; i< num_wire; i++){ //num_wire
	  thz_Sense = phi+i*3*theta_ring;
	  pMatrix_SenseII->SetAngles(thz_Sense/gra2deg,0.0,0.0);
	  if (iring != 0 && iring!=nring) 
	    storeWireData->InsertWireMatrix(iring-1,i,new TGeoHMatrix((*pMatrix_SenseII)*(*comb)*(*pMatrix_epsilon1)));
	}
	if(iring!=nring){
	  radius_ring = radius_ring2;
	  radius_ring_0 = radius_ring2_0;  
	}
    }
  
    radius_ring+=fieldwire_diameter;
    radius_ring_0+=fieldwire_diameter;
  
    storeWireData->FillData();
  }
  
  storeWireData->WriteData();
  delete storeWireData;

}

//_____________________________________________________________________________ 
void IlcDCHgeometry::CreateGeometryCluCou(Int_t *idtmed) 
{ 

 // TGeoMedium *gas_mix=gGeoManager->GetMedium("DCH_Gas-mix"); 
 // TGeoMedium *Al=gGeoManager->GetMedium("DCH_Al"); 
 // TGeoMedium *Kapton=gGeoManager->GetMedium("DCH_Kapton"); 
 // TGeoMedium *Air=gGeoManager->GetMedium("DCH_Air"); 
 // TGeoMedium *CarbonFiber=gGeoManager->GetMedium("DCH_CarbonFiber"); 
 // TGeoMedium *Tungsten=gGeoManager->GetMedium("DCH_Tungsten"); 
 // TGeoMedium *Polypropylene=gGeoManager->GetMedium("DCH_Polypropylene"); 
 // TGeoMedium *Vacuum=gGeoManager->GetMedium("DCH_Vacuum"); 


 TGeoMedium *FWmed=0;

  if ( !fDCHParam ) {
    IlcError("No Parameters found!!!!");
    IlcError("Geometry construction aborted!!!!");
    return;
  }
//  TGeoVolume *topILC=0x0;
  TGeoVolume *top=0x0;

  if (fDCHParam->GetExperiment()==0) {
    top = gGeoManager->GetVolume("ILCC");
    if ( !top ) {
      top = gGeoManager->MakeBox("TOP", 0, 1000., 1000.,1000.);
      gGeoManager->SetTopVolume(top);
    }
  }
  else if (fDCHParam->GetExperiment()==1) {

    // topILC = gGeoManager->GetVolume("ILCC");
    // if ( !topILC ) {
    //   topILC = gGeoManager->MakeBox("TOP", 0, 1000., 1000.,1000.);
    //   gGeoManager->SetTopVolume(topILC);
    // }

    // top = new TGeoVolume("DCH_top",topILC->GetShape(),Vacuum);
    // topILC->AddNode(top,1,0);
  }

  Double_t gra2deg = TMath::Pi()/180.;
  
  // if (fDCHParam->GetFWMaterialType()==0) FWmed = Al;
  // else if (fDCHParam->GetFWMaterialType()==1) FWmed = Kapton;

  
  //Int_t Flag_geoDCHg4=fDCHParam->GetFlag_geoDCHg4(); //If the flag is 0  you run Axial DCH with Geant4
  Double_t inner_radius=fDCHParam->GetInnerRadius(); //142.5549;
//  Double_t endcap_inner_radius=inner_radius;
  Double_t outer_radius=fDCHParam->GetOuterRadius(); //148.0343;
  Double_t fieldwire_diameter=fDCHParam->GetFWireDiameter();//.01;
//  Double_t sensewire_diameter=fDCHParam->GetSWireDiameter(); //.0025;
  Double_t envelop_Inner_thickness=fDCHParam->GetInnerWallThickness();//0.02;
//  Double_t envelop_Outer_thickness=fDCHParam->GetOuterWallThickness(); //.;
  Double_t envelop_EndCap_thickness=fDCHParam->GetEndCapWallThickness();//2.;
  Double_t extra_EndCap_dist = fDCHParam->GetExtraEndCapDist();
  Double_t length=fDCHParam->GetLength(); //0.5*300.;
  
  Int_t num_wire_sense=fDCHParam->GetSWireNum(); //80;
  
  Int_t delta_num_wire_sense=fDCHParam->GetSDeltaWireNum(); //20;
  Int_t nsuperlayer=fDCHParam->GetSuperLayerNum(); //20;
  Int_t nring=fDCHParam->GetRingNum(); //11;
  Double_t drop=fDCHParam->GetDrop(); //2.;
  
  //Double_t extra_EndCap_dist=fDCHParam->GetExtraEndCapDist();  
  
  
  Double_t extra_dim=0.0;//fDCHParam->GetExtraDim(); //100.;
  Double_t EndCap_Wall_theta_inner=fDCHParam->GetEndCapWallThetaIn(); //0.;
  Double_t EndCap_Wall_theta_outer=fDCHParam->GetEndCapWallThetaOut();//0.;

  
  Int_t EndCap_type=fDCHParam->GetEndCapType();

  Double_t max_EndCap_dim=fDCHParam->GetMaxEndCapDim(); //212. used if EndCap_type!=0

  //if(Flag_geoDCHg4==0) drop=0;
  
  if(EndCap_type==0) {
    length = length-envelop_EndCap_thickness;
    extra_EndCap_dist=0.;
  }
  else if(EndCap_type==1){
    length=(max_EndCap_dim-envelop_EndCap_thickness)*length/max_EndCap_dim; ///TMath::Sqrt(2.);

/*    extra_EndCap_dist=TMath::Sqrt(pow(max_EndCap_dim-envelop_EndCap_thickness,2)-pow(inner_radius,2))-length;
    // extra_EndCap_dist=(max_EndCap_dim-envelop_EndCap_thickness)-length;
    EndCap_Wall_theta_outer=TMath::ACos(length/(max_EndCap_dim-envelop_EndCap_thickness))/gra2deg;
    //EndCap_Wall_theta_inner=TMath::ACos((length+extra_EndCap_dist)/(max_EndCap_dim-envelop_EndCap_thickness))/gra2deg;
    EndCap_Wall_theta_inner=TMath::ASin((inner_radius+ envelop_Inner_thickness)/(max_EndCap_dim-envelop_EndCap_thickness))/gra2deg;    
*/
  }

  Double_t z=0.;
  Double_t radius_ring, radius_ring_0, radius_ring2, radius_ring2_0, theta_ring, delta_radius_ring, ringangle;
  Double_t alfa1, alfa2;
  Double_t epsilon1, epsilon2;
  Double_t thz_Sense;
  Int_t num_wire/*,color_ring1,color_ring2,idcopy_cell,idcopy_fw1,idcopy_fw2*/;
//  Int_t count_node, count_subsector;
  Double_t kbase_exag=0;
  Char_t shape_name_FD[50], shape_name_SD[50];
  
  //TGeoVolumeAssembly *sector;

  radius_ring_0 = inner_radius + envelop_Inner_thickness + 0.5*fieldwire_diameter;
  radius_ring = radius_ring_0 + drop;
  radius_ring2 = radius_ring;
  radius_ring2_0 = radius_ring_0;
  TGeoRotation *pMatrix1_SenseII = new TGeoRotation("rot1_SenseII");
  TGeoRotation *pMatrix2_SenseII = new TGeoRotation("rot2_SenseII");
  TGeoRotation *pMatrix3_SenseII = new TGeoRotation("rot3_SenseII");
  TGeoRotation *pMatrix_epsilon1 = new TGeoRotation("epsilon1");
  TGeoRotation *pMatrix_epsilon2 = new TGeoRotation("epsilon2");
  
  Char_t name_tr_rot_SenseII[50];


  TGeoShape *shape_Field;
//  TGeoShape *shape_Sense;

  TGeoShape *shape_Cut=0;
  TGeoTranslation *posCut_EndCap=0;

  if(EndCap_type==0){
    shape_Cut = new TGeoTube("EndCap_cut",inner_radius,outer_radius, extra_dim);
    posCut_EndCap = new TGeoTranslation("TR_EndCap",0.,0., /*extra_dim+*/length);
    posCut_EndCap->RegisterYourself();
  }

  else if (EndCap_type==1){
    shape_Cut = new TGeoSphere("EndCap_cut",max_EndCap_dim-envelop_EndCap_thickness,max_EndCap_dim/*+extra_dim*/,EndCap_Wall_theta_inner,EndCap_Wall_theta_outer,0.,360.);
    posCut_EndCap = new TGeoTranslation("TR_EndCap",0.,0.,0.);
    posCut_EndCap->RegisterYourself();
  }
  

  TGeoVolume *vol_tube_FD=0;
//  TGeoVolume *vol_tube_SD=0;

  TGeoRotation *leftEndCap = new TGeoRotation("rot_EndCap",0.,180.,0.);
  TGeoCombiTrans *comb_tr_rot_EndCap = new TGeoCombiTrans("tr_r_EndCap");
  (*comb_tr_rot_EndCap) = (*leftEndCap) * (*posCut_EndCap);
  TGeoHMatrix *tr_rot_EndCap = new TGeoHMatrix((*comb_tr_rot_EndCap));
  tr_rot_EndCap->SetName("TR_ROT_EndCap");
  tr_rot_EndCap->RegisterYourself();

  Int_t sign_epsilon=-1;



//******************************new geometry***********
//  TGeoShape*wring;
//  TGeoShape*gring;
//  TGeoVolume*vol_gring;
//  TGeoVolume*vol_wring;
  Char_t wshape[50];
  Char_t gshape[50];
  Char_t wvol[50];
  Char_t gvol[50];
    
//  Double_t wre[20];
//  Double_t epsilonRing[20];
  Double_t radiuss=0;
  Double_t Drd=0.000001;
  Double_t lenwire_endcap=0.,zlenwire_endcap=0.;
//  Double_t wri=0.;
//  Double_t gri=0.;
//  Double_t gre=0.;
  Double_t rwire=0.;
  Double_t phi=0.;
      //******************************************end********

  


  
  IlcDCHwireposition *storeWireData = new IlcDCHwireposition(kTRUE);
  
  Int_t ncel;
  for (Int_t superlayer=0;superlayer<nsuperlayer;superlayer++) {
    cout <<"Building super layer: "<<superlayer+1<<endl;

//    count_node = 0;
//    count_subsector = 0;
    //sector = new TGeoVolumeAssembly(Form("SuperLayer_%02d_%02d",superlayer,count_subsector));
     num_wire=num_wire_sense+superlayer*delta_num_wire_sense;
    
     if(superlayer>=10 && superlayer<18) num_wire=num_wire-4*delta_num_wire_sense;
     if(superlayer>=18 && superlayer<24) num_wire=num_wire-6*delta_num_wire_sense;
    
     //     cout<<" SL "<<superlayer<<" nwire "<<num_wire<<endl;
    theta_ring=2.*(Double_t) TMath::Pi()/(3*num_wire);


    phi=2.*(Double_t) TMath::Pi()/(num_wire);
    
    sign_epsilon*=-1;

    
    ncel = nring-1;
    storeWireData->AddSuperLayer( ncel, num_wire);

    for(Int_t iring=0; iring< nring+1 ; iring++){

      if ((iring%2)==0){
  	ringangle = 0.;
//  	color_ring1=1 + (((superlayer%2)==0) ? 0 : 3) ;
//  	color_ring2=2 + (((superlayer%2)==0) ? 0 : 3) ;
      }
      else{
  	ringangle = -(1.5*theta_ring);
//  	color_ring1=2 + (((superlayer%2)==0) ? 0 : 3) ;
//  	color_ring2=1 + (((superlayer%2)==0) ? 0 : 3) ;
      }

      kbase_exag = 2.*radius_ring_0*TMath::Sin(theta_ring*0.5);
      fBaseExag[superlayer] = kbase_exag;
      delta_radius_ring = kbase_exag * TMath::Cos(30. * gra2deg);
      if(iring!=nring){
	radius_ring2 += delta_radius_ring;
	radius_ring2_0 += delta_radius_ring;
      }
      alfa1 = TMath::ACos(1.-(drop/radius_ring));
      alfa2 = TMath::ACos(1.-(drop/radius_ring2));
      epsilon1 = TMath::ATan(radius_ring/length*TMath::Sin(alfa1));
      epsilon2 = TMath::ATan(radius_ring2/length*TMath::Sin(alfa2));
//      epsilonRing[iring]=epsilon1;
//      epsilonRing[iring+1]=epsilon2;
      pMatrix_epsilon1->SetAngles(0.0,-sign_epsilon*epsilon1/gra2deg,0.0); 
      pMatrix_epsilon2->SetAngles(0.0,-sign_epsilon*epsilon2/gra2deg,0.0); 
      

      Double_t cord[3]= {radius_ring_0+2*fieldwire_diameter,0.,0.};
      Double_t dir[3]={0,TMath::Sin(epsilon1),TMath::Cos(epsilon1)};
      Int_t iact=1;
      Double_t step=TGeoShape::Big();
      Double_t* safe = 0;
      if(EndCap_type==0) { 
        zlenwire_endcap=length;
	lenwire_endcap=zlenwire_endcap/TMath::Cos(epsilon1);
	}
      else if (EndCap_type==1) {
        lenwire_endcap=shape_Cut->DistFromOutside(cord,dir,iact,step, safe);
        zlenwire_endcap=TMath::Cos(epsilon1)*lenwire_endcap;
      }
      lenwire_endcap-=fieldwire_diameter*TMath::Sin(epsilon1);


      sprintf(wshape,"wS%dR%d",superlayer,iring);
      sprintf(gshape,"gS%dR%d",superlayer,iring);
      sprintf(wvol,"wvolS%02dR%02d",superlayer,iring);
      sprintf(gvol,"gvolS%02dR%02d",superlayer,iring);

      if(iring==0){
	radiuss=radius_ring_0;
//	wri=radiuss+Drd;
//	wre[iring]=radiuss+2*fieldwire_diameter;
	rwire=radiuss+fieldwire_diameter;
	//	wring=new TGeoTube(wshape,wri,wre[iring],lenwire_endcap);
	// wring=new TGeoHype(wshape,wri,epsilon1/gra2deg,wre[iring],epsilon1/gra2deg,zlenwire_endcap);
	// vol_wring=new TGeoVolume(wvol,wring,gas_mix);
	// vol_wring->SetLineColor(5);
	// vol_wring->SetVisibility(1);
	// top->AddNode(vol_wring,1,0);
      }else{ 
	radiuss=radius_ring_0;
	if(iring!=nring){
//	  wri=radiuss-fieldwire_diameter;
//	  wre[iring]=radiuss+fieldwire_diameter;
	  rwire=radiuss;
	}else{
//	  wri=radiuss-2*fieldwire_diameter;
//	  wre[iring]=radiuss;
	  rwire=radiuss-fieldwire_diameter;
	}
//	gri=wre[iring-1]+Drd;
//	gre=wri-Drd;

	//	gring=new TGeoTube(gshape,gri,gre,lenwire_endcap);
	// gring=new TGeoHype(gshape,gri,epsilonRing[iring-1]/gra2deg,gre,epsilon1/gra2deg,zlenwire_endcap);
	// vol_gring=new TGeoVolume(gvol,gring,gas_mix);
	// vol_gring->SetLineColor(3);
	// vol_gring->SetVisibility(1);
	// top->AddNode(vol_gring,1,0);

	//	wring=new TGeoTube(wshape,wri,wre[iring],lenwire_endcap);
	// wring=new TGeoHype(wshape,wri,epsilon1/gra2deg,wre[iring],epsilon1/gra2deg,zlenwire_endcap);
	// vol_wring=new TGeoVolume(wvol,wring,gas_mix);
	// vol_wring->SetLineColor(5);
	// vol_wring->SetVisibility(1);
	// top->AddNode(vol_wring,1,0);
     
      }
      if(iring==0 && superlayer!=0){
//	gri=radius_ring_0-fieldwire_diameter;
//	gre=radius_ring_0;
	sprintf(gvol,"gvolE%02dR%02d",superlayer,iring+1);
	// cout<<" ************* "<<gvol<<end;
	//	gring=new TGeoTube(gshape,gri,gre,lenwire_endcap-Drd2);
	// gring=new TGeoHype(gshape,gri,epsilonRing[nring]/gra2deg,gre,epsilon1/gra2deg,zlenwire_endcap);
	// vol_gring=new TGeoVolume(gvol,gring,gas_mix);
	// vol_gring->SetLineColor(4);
	// vol_gring->SetVisibility(1);
	// top->AddNode(vol_gring,1,0);

      }
     

//	if(superlayer==nsuperlayer-1 && iring==nring ) {
//	  gri=wre[iring]+Drd;//+2*delta_radius_ring;
//	  gre=outer_radius-envelop_Outer_thickness;
//
//	  //	  gring=new TGeoTube(gshape,gri,gre,length);
//	  // gring=new TGeoHype(gshape,gri,epsilon1/gra2deg,gre,0.,length);
//	  // vol_gring=new TGeoVolume("gvoES23R11",gring,gas_mix);
//	  // vol_gring->SetLineColor(2);
//	  // vol_gring->SetVisibility(1);
//	  // top->AddNode(vol_gring,1,0);
//	}
//
//	if(superlayer==0 && iring==0) {
//	  gri=inner_radius + envelop_Inner_thickness;
//	  gre=radius_ring_0;
//
//	  //	  gring=new TGeoTube(gshape,gri,gre,length);
//	  // gring=new TGeoHype(gshape,gri,0,gre,epsilon1/gra2deg,zlenwire_endcap);
//	  // vol_gring=new TGeoVolume("gvoES00R00",gring,gas_mix);
//	  // vol_gring->SetLineColor(2);
//	  // vol_gring->SetVisibility(1);
//	  // top->AddNode(vol_gring,1,0);
//	}
     
	//cout<<gre<<endl;
	//matrix cell shape for geant4------------------------------------

	TGeoCombiTrans *pMat;
	TGeoHMatrix *comb;

	phi=ringangle;

	pMat = new TGeoCombiTrans("id_comb_Sense",rwire,0.,z,0);
	comb= new TGeoHMatrix((*pMat));
	comb->SetName("rot");
	comb->RegisterYourself();

     
	if (iring != 0 && iring!=nring) {
	  storeWireData->InsertAlfaEpsilon(iring-1, alfa1,sign_epsilon*epsilon1);
	  storeWireData->InsertRadius(iring-1, radius_ring_0);
	}
	if(iring==nring)
	  cout<<"Max R on the EndCap for Layer "<<rwire<<endl;




	//---------------End--------------------------------------------
      shape_Field=new TGeoTube("tube_Field", 1.0E-10, fieldwire_diameter*0.5, lenwire_endcap);
//      shape_Sense=new TGeoTube("tube_Sense", 1.0E-10, sensewire_diameter*0.5, lenwire_endcap);

      if(iring!=nring){
	sprintf(shape_name_FD,"tubeFD_%d_%d",superlayer,iring);
	vol_tube_FD=new TGeoVolume(shape_name_FD,shape_Field,FWmed);
	vol_tube_FD->SetLineColor(3);
	vol_tube_FD->SetVisibility(1);
      }

      if(iring!=0&&iring!=nring){
	sprintf(shape_name_SD,"tubeSD_%d_%d",superlayer,iring);
	// vol_tube_SD=new TGeoVolume(shape_name_SD,shape_Sense,Tungsten);
	// vol_tube_SD->SetLineColor(2);
      }

      for (Int_t i = 0; i< num_wire; i++){ //num_wire

  	thz_Sense = phi+i*3*theta_ring;
      
  	sprintf(name_tr_rot_SenseII,"combTR_SenseII_%d_%d_%d",superlayer,iring,i);
  	
  	pMatrix1_SenseII->SetAngles(thz_Sense/gra2deg,0.0,0.0);
  	pMatrix2_SenseII->SetAngles((thz_Sense+theta_ring)/gra2deg,0.0,0.0);
  	pMatrix3_SenseII->SetAngles((thz_Sense+2*theta_ring)/gra2deg,0.0,0.0);

	//matrix of  cell volume for geant4---------------------- 

	TGeoHMatrix *vcomb_tr_rot_SenseII;
	TGeoHMatrix *v2comb_tr_rot_SenseII;
	TGeoHMatrix *v3comb_tr_rot_SenseII;

	vcomb_tr_rot_SenseII = new TGeoHMatrix((*pMatrix1_SenseII)*(*comb)*(*pMatrix_epsilon1));
  	vcomb_tr_rot_SenseII->SetName(name_tr_rot_SenseII);
  	vcomb_tr_rot_SenseII->RegisterYourself();

	v2comb_tr_rot_SenseII = new TGeoHMatrix((*pMatrix2_SenseII)*(*comb)*(*pMatrix_epsilon1));
  	v2comb_tr_rot_SenseII->SetName(name_tr_rot_SenseII);
  	v2comb_tr_rot_SenseII->RegisterYourself();

	v3comb_tr_rot_SenseII = new TGeoHMatrix((*pMatrix3_SenseII)*(*comb)*(*pMatrix_epsilon1));
  	v3comb_tr_rot_SenseII->SetName(name_tr_rot_SenseII);
  	v3comb_tr_rot_SenseII->RegisterYourself();
  
	//end----------------------------------------------------------


	if (iring != 0 && iring!=nring) 
	  storeWireData->InsertWireMatrix(iring-1,i,new TGeoHMatrix(*vcomb_tr_rot_SenseII));
	
//        idcopy_cell = i+1;
//
//        idcopy_fw1 = 2*i+1;
//
//        idcopy_fw2 = idcopy_fw1 + 1;
//	if (iring != 0 && iring!=nring)
//	  //vol_wring->AddNode(vol_tube_SD,idcopy_cell,vcomb_tr_rot_SenseII);
//	if (iring!=nring){
//	  //vol_wring->AddNode(vol_tube_FD,idcopy_fw1,v2comb_tr_rot_SenseII);
//	  //vol_wring->AddNode(vol_tube_FD,idcopy_fw2,v3comb_tr_rot_SenseII);
//	}

      }
      if(iring!=nring){
	radius_ring = radius_ring2;
	radius_ring_0 = radius_ring2_0;  
      }
    }// end cycle over rings
  
    radius_ring+=fieldwire_diameter;
    radius_ring_0+=fieldwire_diameter;
    storeWireData->FillData();
  }
  
  storeWireData->WriteData();
  delete storeWireData;

  double lengthInnerWall=length+extra_EndCap_dist;
  if(EndCap_type==1){
    double zl=TMath::Sqrt(TMath::Power(max_EndCap_dim-envelop_EndCap_thickness,2)
			  -TMath::Power(inner_radius+envelop_Inner_thickness,2));
    if(zl<lengthInnerWall) lengthInnerWall=zl-Drd;
  }

  if (fDCHParam->GetExperiment()==0) {
    // TGeoShape *innerWall=new TGeoTube("InnerWall",inner_radius,inner_radius+envelop_Inner_thickness, lengthInnerWall);

    // TGeoVolume *vol_innerWall=new TGeoVolume("InnerWall",innerWall,CarbonFiber);
    // vol_innerWall->SetLineColor(10);
    // vol_innerWall->SetVisibility(1);
    // top->AddNode(vol_innerWall,1,0);
    // if (gIlc->GetModule("VXD")==0x0) {
    //   TGeoVolume *innerHole=new TGeoVolume("InnerHole",new TGeoTube("InnerHole",0.0,inner_radius, lengthInnerWall),Air);
    //   top->AddNode(innerHole,1,0);
    // }
  }
  else if (fDCHParam->GetExperiment()==1) {
//    if (fDCHParam->GetExperimentSubVer()==0){
//      // TGeoVolume *innerHole=new TGeoVolume("InnerHole",new TGeoTube("InnerHole",0.0,inner_radius+envelop_Inner_thickness, lengthInnerWall),gas_mix);
//      // top->AddNode(innerHole,1,0);
//      endcap_inner_radius=0.0;
//    }
//    else if (fDCHParam->GetExperimentSubVer()==1){
//      // TGeoShape *innerWall=new TGeoTube("InnerWall",inner_radius,inner_radius+envelop_Inner_thickness, lengthInnerWall);
//      // TGeoShape *innerWall_foam=new TGeoTube("InnerWall_foam",inner_radius-0.48,inner_radius, lengthInnerWall);
//      // TGeoShape *innerWall_1=new TGeoTube("InnerWall_1",inner_radius-0.48-envelop_Inner_thickness,inner_radius-0.48, lengthInnerWall);
//
//      // TGeoVolume *vol_innerWall=new TGeoVolume("InnerWall",innerWall,CarbonFiber);/*Kapton*/
//      // TGeoVolume *vol_innerWall_foam=new TGeoVolume("InnerWall_foam",innerWall_foam,Polypropylene);
//      // TGeoVolume *vol_innerWall_1=new TGeoVolume("InnerWall_1",innerWall_1,CarbonFiber);
//      // TGeoVolume *innerHole=new TGeoVolume("InnerHole",new TGeoTube("InnerHole",0.0,inner_radius-0.48-envelop_Inner_thickness, lengthInnerWall),Vacuum);
//
//      // vol_innerWall->SetLineColor(10);
//      // vol_innerWall->SetVisibility(1);
//      // top->AddNode(vol_innerWall,1,0);
//      // //  vol_innerWall_foam->SetLineColor(10);
//      // vol_innerWall_foam->SetVisibility(1);
//      // top->AddNode(vol_innerWall_foam,1,0);
//      //  vol_innerWall_1->SetLineColor(10);
//      // vol_innerWall_1->SetVisibility(1);
//      // top->AddNode(vol_innerWall_1,1,0);
//
//      // top->AddNode(innerHole,1,0);
//    }

    //mu2e absorber

    //TGeoMedium *Polyethilene=gGeoManager->GetMedium("DCH_Radiator"); 

    //TGeoShape *absorber=new TGeoTube("absorber",outer_radius+2.,outer_radius+22., length+extra_EndCap_dist);
    //TGeoVolume *vol_absorber=new TGeoVolume("vol_absorber",absorber,Polyethilene);

    //vol_absorber->SetLineColor(3);
    //vol_absorber->SetVisibility(1);
    //top->AddNode(vol_absorber,1,0);
  }
  /*
  TGeoShape *outerWall=new TGeoTube("OuterWall",outer_radius-envelop_Outer_thickness,outer_radius, length);

  TGeoVolume *vol_outerWall=new TGeoVolume("OuterWall",outerWall,CarbonFiber);
  vol_outerWall->SetLineColor(10);
  vol_outerWall->SetVisibility(1);
  top->AddNode(vol_outerWall,1,0);
  TGeoShape *EndCapWall;
  TGeoVolume *vol_EndCapWall_Right;
  TGeoVolume *vol_EndCapWall_Left;
  if(EndCap_type==0){
    EndCapWall=new TGeoTube("EndCapWall",endcap_inner_radius,outer_radius, 0.5*envelop_EndCap_thickness);
    vol_EndCapWall_Right=new TGeoVolume("RightWall",EndCapWall,CarbonFiber);

    TGeoTranslation *pos_EndCapWall = new TGeoTranslation("TR_EndCapWall",0.,0., length + 0.5*envelop_EndCap_thickness);
    TGeoTranslation *pos_EndCapWallL = new TGeoTranslation("TR_EndCapWallL",0.,0., -length -0.5*envelop_EndCap_thickness);
    vol_EndCapWall_Right->SetLineColor(10);
    vol_EndCapWall_Right->SetVisibility(1);
    top->AddNode(vol_EndCapWall_Right,1,pos_EndCapWall);

    vol_EndCapWall_Left=new TGeoVolume("LeftWall",EndCapWall,CarbonFiber);
    vol_EndCapWall_Left->SetLineColor(10);
    vol_EndCapWall_Left->SetVisibility(1);
    top->AddNode(vol_EndCapWall_Left,1,pos_EndCapWallL);
  }

  else if(EndCap_type==1){
    EndCapWall=new TGeoSphere("EndCapWall",max_EndCap_dim-envelop_EndCap_thickness,max_EndCap_dim,EndCap_Wall_theta_inner,EndCap_Wall_theta_outer,0.,360.);
    vol_EndCapWall_Right=new TGeoVolume("RightWall",EndCapWall,CarbonFiber);
    vol_EndCapWall_Right->SetLineColor(3);
    vol_EndCapWall_Right->SetVisibility(1);
    top->AddNode(vol_EndCapWall_Right,1,0);

    vol_EndCapWall_Left=new TGeoVolume("LeftWall",EndCapWall,CarbonFiber);
    vol_EndCapWall_Left->SetLineColor(3);
    vol_EndCapWall_Left->SetVisibility(1);
    top->AddNode(vol_EndCapWall_Left,1,leftEndCap);
  }
  */
//  gGeoManager->CloseGeometry();
//  gGeoManager->CheckOverlaps();
//  top->Raytrace();

/*      if ((s%2)==0)
       vol_trap->SetLineColor(2);
      else
       vol_trap->SetLineColor(1);*/


}


//_____________________________________________________________________________ 
void IlcDCHgeometry::CreateGeometry(Int_t *idtmed) 
{ 
/*   // 
  // Create the DCH geometry without hole 
  // 
  // 
  // Names of the DCH volumina (xx = detector number): 
  // 
  //      Volume (Air) wrapping the readout chamber components 
  //        UTxx    includes: UAxx, UDxx, UFxx, UUxx 
  //      Obs: 
  //        UUxx    the services volume has been reduced by 7.42 mm 
  //                in order to allow shifts in radial direction 
  // 
  //      Lower part of the readout chambers (gas volume + radiator) 
  // 
  //        UAxx    Aluminum frames             (Al) 
  //        UBxx    G10 frames                  (C) 
  //        UCxx    Inner volumes               (Air) 
  // 
  //      Upper part of the readout chambers (readout plane + fee) 
  // 
  //        UDxx    G10 frames                  (C) 
  //        UExx    Inner volumes of the G10    (Air) 
  //        UFxx    Aluminum frames             (Al) 
  //        UGxx    Inner volumes of the Al     (Air) 
  // 
  //      Inner material layers 
  //0 
  //        UHxx    Radiator                    (Rohacell) 
  //        UIxx    Entrance window             (Mylar) 
  //        UJxx    Drift volume                (Xe/CO2) 
  //        UKxx    Amplification volume        (Xe/CO2) 
  //        ULxx    Pad plane                   (Cu) 
  //        UMxx    Support structure           (Rohacell) 
 
 
  //BTR1 mother volumes 
  Char_t  cTagV[5]; 
  Char_t  cTagVB[5]; 

  Float_t parTrdVM[3]; 
  parTrdVM[0] = fgkRmin; 
  parTrdVM[1] = 600;//kMAGr2; 
  parTrdVM[2] = fgkZmax1; 

  gMC->Gsvolu("BTR1","TUBE",idtmed[1302-1],parTrdVM,3); 
  gMC->Gspos("BTR1",1, "ILCC",0.,0.,0.,0, "ONLY")  ; 
  
  const Float_t krad2deg = 180./TMath::Pi(); 
  const Float_t kdeg2rad = 1./krad2deg; 
  
  Float_t pmud1[4]; 

///Big Trapezius for one sector
  pmud1[0] =fgkTRwidth1/2.; 
  pmud1[1] =fgkTRwidth2/2.; 
  pmud1[2] =fgkSlenTR2/2;  
  pmud1[3] =fgkSheight[0]/2.; 

 
  gMC->Gsvolu("BTR2", "TRD1", idtmed[1302-1], pmud1, 4); 
 
 
///Big Box for one sector
  pmud1[0] =fgkBwidth1/2; 
  pmud1[1] =fgkSheight[1]/2; 
  pmud1[2] =fgkSlenTR2/2;  
  pmud1[3] =0.; 
 
  gMC->Gsvolu("Box1", "BOX", idtmed[1302-1], pmud1, 3); 

  Float_t r      = (fgkRmax - fgkRmin)/2.; 
  Float_t dx, dy,dx1,dy1; 
  Int_t idrotm[18]; 
  Int_t idrotm_1[18]; 
  Int_t isec= Nsect();  
  Int_t idet= Ncham(); 
  for (Int_t j=0; j<idet; j++){ 
     for (Int_t i=0; i<isec; i+=2){ 
       Float_t phi = i*360/isec; 
       Float_t phi2 =270+phi; 
       if (phi2 >= 360.) phi2-=360.;  
       dx =  TMath::Sin(phi*kdeg2rad)*(fgkRmin*TMath::Cos(GetAlpha1())/TMath::Cos(GetAlpha2())+r*TMath::Cos(30.*kdeg2rad)); 
       dy = -TMath::Cos(phi*kdeg2rad)*(fgkRmin*TMath::Cos(GetAlpha1())/TMath::Cos(GetAlpha2())+r*TMath::Cos(30.*kdeg2rad)); 
       gMC->Matrix(idrotm[i],  90.0, phi, 0., 0., 90., phi2);
       gMC->Gspos("BTR2", j*isec+i+1, "BTR1", dx, dy,(j-1)*fgkSlenTR1/3, idrotm[i], "ONLY")  ; 
       Float_t phi1=phi+30; 
 
       dx1 =  TMath::Sin(phi1*kdeg2rad)*(r+fgkRmin);
       dy1 = -TMath::Cos(phi1*kdeg2rad)*(r+fgkRmin);  
       gMC->Matrix(idrotm_1[i],  90.,phi1, 90.,90.+phi1,0.,0.); 
       gMC->Gspos("Box1",j*isec+i+2, "BTR1", dx1, dy1,(j-1)*fgkSlenTR1/3 , idrotm_1[i], "ONLY"); 
  } 
  } 
 
  Float_t xpos, ypos, zpos; 
  Float_t parTrd[3]; 
  Float_t parCha1[4]; 
  Float_t parBCha[3]; 

  for (Int_t ip=0; ip<kNplan; ip++) { 
     parCha1[0] =fCWidth[ip][0]/2; 
     parCha1[1] =(fCWidth[ip][0]+Cwidcha())/2; 
     parCha1[3] =Cheight()/2; 
     parCha1[2] =GetChamberLength(ip,0);  
 
     sprintf(cTagV,"CTR%01d",ip); 
      zpos=+ip*Cheight()+Cheight()/2 -fgkTRheight/2; 
     // zpos=-ip*fgkTRheight/4+fgkTRheight/2-Cheight()/2; 
     gMC->Gsvolu(cTagV,"TRD1",idtmed[1302-1],parCha1,4); 
     gMC->Gspos(cTagV, 1, "BTR2", 0., 0.,zpos, 0, "ONLY")  ; 
   
     parBCha[0] =GetChamberWidth(ip,1)/2; 
     parBCha[1] =Cheight(1)/2; 
     parBCha[2] =GetChamberLength(ip,0);  
 
     sprintf(cTagV,"CBx%01d",ip); 
     zpos=+ip*fgkSheight[1]/4-fgkSheight[1]/2+Cheight(1)/2;
     gMC->Gsvolu(cTagV,"BOX",idtmed[1302-1],parBCha,3); 
     gMC->Gspos(cTagV, 1, "Box1", 0.,zpos,0., 0, "ONLY"); 
 
} 
 
  // The TUBE (Al) 
parTrd[0] =0; 
parTrd[1] = fgkCamRmax; 
parTrd[2] = fgkSlenTR2/2;  
gMC->Gsvolu("TTR1","TUBE",idtmed[1301-1],parTrd,3); 
 // The gas-mixure tube (He(90%)C4H10(10) 
parTrd[0] =0; 
parTrd[1] =fgkCamRmin; 
parTrd[2] =fgkSlenTR2/2;  
gMC->Gsvolu("DTR1","TUBE",idtmed[1309-1],parTrd,3); 
gMC->Gspos("DTR1",1, "TTR1", 0, 0, 0, 0, "ONLY"); 
 
 // The sense wire(W) 
parTrd[0] =0; 
parTrd[1] =fgkCroW ; 
parTrd[2] =fgkSlenTR2/2; 
gMC->Gsvolu("TTR2","TUBE",idtmed[1316-1],parTrd,3); 
gMC->Gspos("TTR2",1, "DTR1", 0, 0, 0, 0, "ONLY"); 
   
     
 Float_t dxx=Step(); 
 Float_t dz=Step(); 
Double_t d =fgkTradius+fgkTradius/TMath::Cos(30.*TMath::Pi()/180.); 

Int_t imax=Ntub();
Int_t jmax=Ntub(1);
Int_t ntp=int(2*(fgkVspace[0]/(dxx*(1+TMath::Cos(30*TMath::Pi()/180.)))));

Int_t f; 
gMC->Matrix(f,  90., 0., 0., 0., 90., 90.); 
for (Int_t ip=kNplan-1;ip>-1;ip--){
     //CTR dispari
   if(ip%2==0) imax=(imax-1)-ntp;

   sprintf(cTagV,"CTR%01d",ip);
   if(ip==1)   imax=imax+1;
   for (Int_t ipt=0;ipt<kNdets;ipt+=2){ 
     for (Int_t i=0;i<imax;i++){
   
       xpos =-(fCWidth[ip][0]+Cwidcha())/2+d+i*dxx+ipt*dxx/2;
       ypos = 0; 
       zpos =-ipt*dz+Cheight()/2-fgkTradius;
       if(ip%2==0){
         xpos =-(imax)*dxx/2-dxx/2+d+i*dxx;
         zpos=zpos-fgkVspace[0]/2;
       }
       gMC->Gspos("TTR1",i+ipt*100+kNdets*100*ip,cTagV,xpos,ypos,zpos,f,"ONLY");
       if(i<imax-1 && ipt<kNdets-1){
       xpos =xpos+dxx/2; 
       zpos=zpos-dz;
       gMC->Gspos("TTR1",i+(ipt+1)*100+kNdets*100*ip,cTagV,xpos,ypos,zpos,f,"ONLY");
       }
     } 
     imax-=2;
   }
   if(ip==1)imax=imax-1; 
 
  //BOXes
   sprintf(cTagVB,"CBx%01d",ip);
   for (Int_t ipt=0;ipt<kNdets;ipt++) { 
     for (Int_t i=0;i<jmax;i++){
       xpos=-GetChamberWidth(ip,1)/2+dxx/2+i*dxx;
       ypos=-ipt*dz+Cheight(1)/2-fgkTradius-fgkTradius/2;
       if(i%2==0)ypos=ypos+fgkTradius;
       zpos=0;
       if(ip%2==0) {
        xpos=-GetChamberWidth(ip,1)/2+dxx/2+i*dxx;
        ypos=-ipt*dz+Cheight(1)/2-dxx-fgkVspace[1]/2-fgkTradius-fgkTradius/2;
        if(i%2==0) ypos=ypos+fgkTradius;
       }
       gMC->Gspos("TTR1",i+ipt*jmax+kNdets*jmax*ip+20000,cTagVB,xpos,ypos,zpos,0,"ONLY"); 
     }}

}
 */
}
/* 
double IlcDCHgeometry::GetXLayer(int nlayer,int plane,int nsec){
  int layer=kNlayers-nlayer-1;
  int sec=nsec+1;
  Float_t dxx=Step(); 
  Float_t dz=Step(); 
  if(sec%2==0){
    double zpos =-layer*dz+Cheight()/2-fgkTradius;
    if(plane%2==0){
      zpos=zpos-fgkVspace[0]/2;
    }
    zpos+=plane*Cheight()+Cheight()/2 -fgkTRheight/2;
    //    zpos+=Rmin()+Sheight(0);
    return zpos;
  }else{
    double ypos=-int(layer/2)*dz+Cheight(1)/2-fgkTradius-fgkTradius/2;
    if(layer%2==0)ypos=ypos+fgkTradius;
    if(plane%2==0) {
      ypos=-int(layer/2)*dz+Cheight(1)/2-dxx-fgkVspace[1]/2-fgkTradius-fgkTradius/2;
      if(layer%2==0) ypos=ypos+fgkTradius;
    }
    ypos+=plane*fgkSheight[1]/4-fgkSheight[1]/2+Cheight(1)/2;
    //    ypos+=Rmin()+Sheight(1);
    return ypos;
  }

}
double IlcDCHgeometry::GetXSector(int sec){
  return Rmin()*(((sec+1)%2==0)?(TMath::Cos(GetAlpha1())/TMath::Cos(GetAlpha2())):1.)+Sheight((sec+1)%2)*0.5;
}
 */
//_____________________________________________________________________________ 
/* 
void IlcDCHgeometry::CreateFrame(Int_t *idtmed) 
{ 
  // 
  // Create the geometry of the frame of the supermodule 
  // 
  // Names of the DCH services volumina 
  // 
  //        USRL    Support rails for the chambers (Al) 
  //        USxx    Support cross bars between the chambers (Al) 
  // 
 
  Int_t   iplan = 0; 
 
  Float_t xpos  = 0.0; 
  Float_t ypos  = 0.0; 
  Float_t zpos  = 0.0; 
 
  Char_t  cTagV[5]; 
 
  // 
  // The chamber support rails 
  // 
 
  const Float_t kSRLwid  = 2.0; 
  const Float_t kSRLhgt  = 2.3; 
  const Float_t kSRLdst  = 0.6; 
  const Int_t   kNparSRL = 3; 
  Float_t parSRL[kNparSRL]; 
  parSRL[0] = kSRLwid/2.; 
  parSRL[1] = fgkSlenTR1/2.; 
  parSRL[2] = kSRLhgt/2.; 
  gMC->Gsvolu("USRL","BOX ",idtmed[1301-1],parSRL,kNparSRL); 
 
  xpos  = 0.0; 
  ypos  = 0.0; 
  zpos  = 0.0; 
  for (iplan = 0; iplan < kNplan; iplan++) { 
     
    xpos  = fCWidth[iplan]/2. + kSRLwid/2. + kSRLdst; 
    ypos  = 0.0; 
    zpos  = fgkCraH + fgkCdrH - fgkSheight[1]/2. - kSRLhgt/2.  
          + iplan * (fgkCH + fgkVspace); 
    //  gMC->Gspos("USRL",iplan+1         ,"UTI1", xpos,ypos,zpos,0,"ONLY"); 
    //gMC->Gspos("USRL",iplan+1+  kNplan,"UTI1",-xpos,ypos,zpos,0,"ONLY"); 
 
  } 
 
  // 
  // The cross bars between the chambers 
  // 
 
  const Float_t kSCBwid  = 1.0; 
  const Int_t   kNparSCB = 3; 
  Float_t parSCB[kNparSCB]; 
  parSCB[1] = kSCBwid/2.; 
  parSCB[2] = fgkCH/2.; 
 
  xpos  = 0.0; 
  ypos  = 0.0; 
  zpos  = 0.0; 
  for (iplan = 0; iplan < kNplan; iplan++) { 
 
    parSCB[0] = fCWidth[iplan]/2. + kSRLdst/2.; 
 
    sprintf(cTagV,"US0%01d",iplan); 
    gMC->Gsvolu(cTagV,"BOX ",idtmed[1301-1],parSCB,kNparSCB); 
    xpos  = 0.0; 
    ypos  =   fgkSlenTR1/2. - kSCBwid/2.; 
    zpos  = fgkCH/2. - fgkSheight/2. + iplan * (fgkCH + fgkVspace); 
    // gMC->Gspos(cTagV,1,"UTI1", xpos,ypos,zpos,0,"ONLY"); 
 
    sprintf(cTagV,"US1%01d",iplan); 
    gMC->Gsvolu(cTagV,"BOX ",idtmed[1301-1],parSCB,kNparSCB); 
    xpos  = 0.0; 
    ypos  = fClength[iplan][2]/2. + fClength[iplan][1]; 
    zpos  = fgkCH/2. - fgkSheight/2. + iplan * (fgkCH + fgkVspace); 
    // gMC->Gspos(cTagV,1,"UTI1", xpos,ypos,zpos,0,"ONLY"); 
 
    sprintf(cTagV,"US2%01d",iplan); 
    gMC->Gsvolu(cTagV,"BOX ",idtmed[1301-1],parSCB,kNparSCB); 
    xpos  = 0.0; 
    ypos  = fClength[iplan][2]/2.; 
    zpos  = fgkCH/2. - fgkSheight/2. + iplan * (fgkCH + fgkVspace); 
    // gMC->Gspos(cTagV,1,"UTI1", xpos,ypos,zpos,0,"ONLY"); 
 
    sprintf(cTagV,"US3%01d",iplan); 
    gMC->Gsvolu(cTagV,"BOX ",idtmed[1301-1],parSCB,kNparSCB); 
    xpos  = 0.0; 
    ypos  = - fClength[iplan][2]/2.; 
    zpos  = fgkCH/2. - fgkSheight/2. + iplan * (fgkCH + fgkVspace); 
    // gMC->Gspos(cTagV,1,"UTI1", xpos,ypos,zpos,0,"ONLY"); 
 
    sprintf(cTagV,"US4%01d",iplan); 
    gMC->Gsvolu(cTagV,"BOX ",idtmed[1301-1],parSCB,kNparSCB); 
    xpos  = 0.0; 
    ypos  = - fClength[iplan][2]/2. - fClength[iplan][1]; 
    zpos  = fgkCH/2. - fgkSheight/2. + iplan * (fgkCH + fgkVspace); 
    //gMC->Gspos(cTagV,1,"UTI1", xpos,ypos,zpos,0,"ONLY"); 
 
    sprintf(cTagV,"US5%01d",iplan); 
    gMC->Gsvolu(cTagV,"BOX ",idtmed[1301-1],parSCB,kNparSCB); 
    xpos  = 0.0; 
    ypos  = - fgkSlenTR1/2. + kSCBwid/2.; 
    zpos  = fgkCH/2. - fgkSheight/2. + iplan * (fgkCH + fgkVspace); 
    // gMC->Gspos(cTagV,1,"UTI1", xpos,ypos,zpos,0,"ONLY"); 
 
  } 
} 
*/ 
 
//_____________________________________________________________________________ 
/* 
void IlcDCHgeometry::CreateServices(Int_t *idtmed) 
{ 
  // 
  // Create the geometry of the services 
  // 
  // Names of the DCH services volumina 
  // 
  //        UTCL    Cooling arterias (Al) 
  //        UTCW    Cooling arterias (Water) 
  //        UUxx    Volumes for the services at the chambers (Air) 
  //        UTPW    Power bars       (Cu) 
  //        UTCP    Cooling pipes    (Al) 
  //        UTCH    Cooling pipes    (Water) 
  //        UTPL    Power lines      (Cu) 
  //        UMCM    Readout MCMs     (G10/Cu/Si) 
  // 
  Int_t   iplan = 0; 
  Int_t   icham = 0; 
 
  Float_t xpos  = 0.0; 
  Float_t ypos  = 0.0; 
  Float_t zpos  = 0.0; 
 
  Char_t  cTagV[5]; 
 
  // The rotation matrices 
  const Int_t kNmatrix = 3; 
  Int_t   matrix[kNmatrix]; 
  gMC->Matrix(matrix[0],100.0,  0.0, 90.0, 90.0, 10.0,  0.0); 
  gMC->Matrix(matrix[1], 80.0,  0.0, 90.0, 90.0, 10.0,180.0); 
  gMC->Matrix(matrix[2],  0.0,  0.0, 90.0, 90.0, 90.0,  0.0); 
 
  IlcDCHCommonParam* commonParam = IlcDCHCommonParam::Instance(); 
  if (!commonParam) 
  { 
    IlcError("Could not get common params\n"); 
    return; 
  } 
     
  // 
  // The cooling arterias 
  // 
 
  // Width of the cooling arterias 
  const Float_t kCOLwid  =  0.5;  
  // Height of the cooling arterias 
  const Float_t kCOLhgt  =  5.5; 
  // Positioning of the cooling  
  const Float_t kCOLposx =  1.6; 
  const Float_t kCOLposz = -0.2; 
  // Thickness of the walls of the cooling arterias 
  const Float_t kCOLthk  =  0.1; 
  const Int_t   kNparCOL = 3; 
  Float_t parCOL[kNparCOL]; 
  parCOL[0]  = kCOLwid/2.; 
  parCOL[1]  = fgkSlenTR1/2.; 
  parCOL[2]  = kCOLhgt/2.; 
  gMC->Gsvolu("UTCL","BOX ",idtmed[1324-1],parCOL,kNparCOL); 
  parCOL[0] -= kCOLthk; 
  parCOL[1]  = fgkSlenTR1/2.; 
  parCOL[2] -= kCOLthk; 
  gMC->Gsvolu("UTCW","BOX ",idtmed[1314-1],parCOL,kNparCOL); 
 
  xpos  = 0.0; 
  ypos  = 0.0; 
  zpos  = 0.0; 
  // gMC->Gspos("UTCW",1,"UTCL", xpos,ypos,zpos,0,"ONLY"); 
 
  for (iplan = 0; iplan < kNplan; iplan++) { // CHECK FOR OVERLAPS !!!  
    //for (iplan = 1; iplan < kNplan; iplan++) { 
     
    xpos  = fCWidth[iplan]/2. + kCOLwid/2. + kCOLposx; 
    ypos  = 0.0; 
    zpos  = kCOLhgt/2. - fgkSheight/2. + kCOLposz + iplan * (fgkCH + fgkVspace); 
    if (iplan == 0) zpos += 0.25;  // To avoid overlaps ! 
    // gMC->Gspos("UTCL",iplan+1         ,"UTI1", xpos,ypos,zpos,matrix[0],"ONLY"); 
    //  gMC->Gspos("UTCL",iplan+1+  kNplan,"UTI1",-xpos,ypos,zpos,matrix[1],"ONLY"); 
 
  } 
 
  // 
  // The power bars 
  // 
 
  const Float_t kPWRwid  =  0.6; 
  const Float_t kPWRhgt  =  4.5; 
  const Float_t kPWRposx =  1.05; 
  const Float_t kPWRposz =  0.9; 
  const Int_t   kNparPWR = 3; 
  Float_t parPWR[kNparPWR]; 
  parPWR[0] = kPWRwid/2.; 
  parPWR[1] = fgkSlenTR1/2.; 
  parPWR[2] = kPWRhgt/2.; 
  gMC->Gsvolu("UTPW","BOX ",idtmed[1325-1],parPWR,kNparPWR); 
   
  for (iplan = 0; iplan < kNplan; iplan++) { // CHECK FOR OVERLAPS !!!  
    //for (iplan = 1; iplan < kNplan; iplan++) { 
     
    xpos  = fCWidth[iplan]/2. + kPWRwid/2. + kPWRposx; 
    ypos  = 0.0; 
    zpos  = kPWRhgt/2. - fgkSheight/2. + kPWRposz + iplan * (fgkCH + fgkVspace); 
    // gMC->Gspos("UTPW",iplan+1         ,"UTI1", xpos,ypos,zpos,matrix[0],"ONLY"); 
    // gMC->Gspos("UTPW",iplan+1+  kNplan,"UTI1",-xpos,ypos,zpos,matrix[1],"ONLY"); 
 
  } 
 
  // 
  // The volumes for the services at the chambers 
  // 
 
  const Int_t kNparServ = 3; 
  Float_t parServ[kNparServ]; 
 
  for (icham = 0; icham < kNcham; icham++) { 
    for (iplan = 0; iplan < kNplan; iplan++) { 
    // Take out upper plane until DCH mothervolume is adjusted 
    //for (iplan = 0; iplan < kNplan-1; iplan++) { 
 
      Int_t iDet = GetDetectorSec(iplan,icham); 
 
      sprintf(cTagV,"UU%02d",iDet); 
      parServ[0] = fCWidth[iplan]/2.; 
      parServ[1] = fClength[iplan][icham]/2. - fgkHspace/2.; 
      parServ[2] = fgkVspace/2. - 0.742/2.;       
      fChamberUUboxd[iDet][0] = parServ[0]; 
      fChamberUUboxd[iDet][1] = parServ[1]; 
      fChamberUUboxd[iDet][2] = parServ[2]; 
       
      gMC->Gsvolu(cTagV,"BOX",idtmed[1302-1],parServ,kNparServ); 
      xpos  = 0.; 
      ypos  = - fClength[iplan][0] - fClength[iplan][1] - fClength[iplan][2]/2.; 
      for (Int_t ic = 0; ic < icham; ic++) { 
        ypos += fClength[iplan][ic];         
      } 
      ypos += fClength[iplan][icham]/2.; 
      zpos  = fgkCH + fgkVspace/2. - fgkSheight/2. + iplan * (fgkCH + fgkVspace); 
      zpos -= 0.742/2.; 
      fChamberUUorig[iDet][0] = xpos; 
      fChamberUUorig[iDet][1] = ypos; 
      fChamberUUorig[iDet][2] = zpos; 
 
    } 
  } 
 
  // 
  // The cooling pipes inside the service volumes 
  // 
 
  const Int_t kNparTube = 3; 
  Float_t parTube[kNparTube]; 
  // The aluminum pipe for the cooling 
  parTube[0] = 0.0; 
  parTube[1] = 0.0; 
  parTube[2] = 0.0; 
  gMC->Gsvolu("UTCP","TUBE",idtmed[1324-1],parTube,0); 
  // The cooling water 
  parTube[0] =  0.0; 
  parTube[1] =  0.2/2.; 
  parTube[2] = -1.; 
  gMC->Gsvolu("UTCH","TUBE",idtmed[1314-1],parTube,kNparTube); 
  // Water inside the cooling pipe 
  xpos = 0.0; 
  ypos = 0.0; 
  zpos = 0.0; 
  //gMC->Gspos("UTCH",1,"UTCP",xpos,ypos,zpos,0,"ONLY"); 
 
  // Position the cooling pipes in the mother volume 
  const Int_t kNpar = 3; 
  Float_t par[kNpar]; 
  for (icham = 0; icham < kNcham;   icham++) { 
    for (iplan = 0; iplan < kNplan; iplan++) { 
    // Take out upper plane until DCH mothervolume is adjusted 
    //for (iplan = 0; iplan < kNplan-1; iplan++) {  
      Int_t   iDet    = GetDetectorSec(iplan,icham); 
      //   Int_t   iCopy   = GetDetector(iplan,icham,0) * 100; 
      Int_t   nMCMrow = commonParam->GetRowMax(iplan,icham,0); 
      Float_t ySize   = (GetChamberLength(iplan,icham) - 2.*fgkRpadW)  
                      / ((Float_t) nMCMrow); 
      sprintf(cTagV,"UU%02d",iDet); 
      for (Int_t iMCMrow = 0; iMCMrow < nMCMrow; iMCMrow++) { 
        xpos   = 0.0; 
        ypos   = (0.5 + iMCMrow) * ySize - 1.9  
               - fClength[iplan][icham]/2. + fgkHspace/2.; 
        zpos   = 0.0 + 0.742/2.;                  
        par[0] = 0.0; 
        par[1] = 0.3/2.; // Thickness of the cooling pipes 
        par[2] = fCWidth[iplan]/2.; 
	// gMC->Gsposp("UTCP",iCopy+iMCMrow,cTagV,xpos,ypos,zpos 
	//  ,matrix[2],"ONLY",par,kNpar); 
      } 
    } 
  } 
 
  // 
  // The power lines 
  // 
 
  // The copper power lines 
  parTube[0] = 0.0; 
  parTube[1] = 0.0; 
  parTube[2] = 0.0; 
  gMC->Gsvolu("UTPL","TUBE",idtmed[1305-1],parTube,0); 
 
  // Position the power lines in the mother volume 
  for (icham = 0; icham < kNcham;   icham++) { 
    for (iplan = 0; iplan < kNplan; iplan++) { 
    // Take out upper plane until DCH mothervolume is adjusted 
    //for (iplan = 0; iplan < kNplan-1; iplan++) {  
      Int_t   iDet    = GetDetectorSec(iplan,icham); 
      Int_t   iCopy   = GetDetector(iplan,icham,0) * 100; 
      Int_t   nMCMrow = commonParam->GetRowMax(iplan,icham,0); 
      Float_t ySize   = (GetChamberLength(iplan,icham) - 2.*fgkRpadW)  
                      / ((Float_t) nMCMrow); 
      sprintf(cTagV,"UU%02d",iDet); 
      for (Int_t iMCMrow = 0; iMCMrow < nMCMrow; iMCMrow++) { 
        xpos   = 0.0; 
        ypos   = (0.5 + iMCMrow) * ySize - 1.0  
               - fClength[iplan][icham]/2. + fgkHspace/2.; 
        zpos   = -0.4 + 0.742/2.; 
        par[0] = 0.0; 
        par[1] = 0.2/2.; // Thickness of the power lines 
        par[2] = fCWidth[iplan]/2.; 
        gMC->Gsposp("UTPL",iCopy+iMCMrow,cTagV,xpos,ypos,zpos 
                          ,matrix[2],"ONLY",par,kNpar); 
      } 
    } 
  } 
 
  // 
  // The MCMs 
  // 
 
  // The mother volume for the MCMs (air) 
  const Int_t kNparMCM = 3; 
  Float_t parMCM[kNparMCM]; 
  parMCM[0] = 3.0/2.; 
  parMCM[1] = 3.0/2.; 
  parMCM[2] = 0.14/2.; 
  gMC->Gsvolu("UMCM","BOX",idtmed[1302-1],parMCM,kNparMCM); 
 
  // The MCM carrier G10 layer 
  parMCM[0] = 3.0/2.; 
  parMCM[1] = 3.0/2.; 
  parMCM[2] = 0.1/2.; 
  gMC->Gsvolu("UMC1","BOX",idtmed[1319-1],parMCM,kNparMCM); 
  // The MCM carrier Cu layer 
  parMCM[0] = 3.0/2.; 
  parMCM[1] = 3.0/2.; 
  parMCM[2] = 0.0162/2.; 
  gMC->Gsvolu("UMC2","BOX",idtmed[1318-1],parMCM,kNparMCM); 
  // The silicon of the chips 
  parMCM[0] = 3.0/2.; 
  parMCM[1] = 3.0/2.; 
  parMCM[2] = 0.003/2.; 
  gMC->Gsvolu("UMC3","BOX",idtmed[1320-1],parMCM,kNparMCM); 
 
  // Put the MCM material inside the MCM mother volume 
  xpos  =  0.0; 
  ypos  =  0.0; 
  zpos  = -0.07      + 0.1/2.; 
  gMC->Gspos("UMC1",1,"UMCM",xpos,ypos,zpos,0,"ONLY"); 
  zpos +=  0.1/2.    + 0.0162/2.; 
  gMC->Gspos("UMC2",1,"UMCM",xpos,ypos,zpos,0,"ONLY"); 
  zpos +=  0.00162/2 + 0.003/2.; 
  gMC->Gspos("UMC3",1,"UMCM",xpos,ypos,zpos,0,"ONLY"); 
 
  // Position the MCMs in the mother volume 
  for (icham = 0; icham < kNcham;   icham++) { 
    for (iplan = 0; iplan < kNplan; iplan++) { 
    // Take out upper plane until DCH mothervolume is adjusted 
    //for (iplan = 0; iplan < kNplan-1; iplan++) {  
      Int_t   iDet    = GetDetectorSec(iplan,icham); 
      Int_t   iCopy   = GetDetector(iplan,icham,0) * 1000; 
      Int_t   nMCMrow = commonParam->GetRowMax(iplan,icham,0); 
      Float_t ySize   = (GetChamberLength(iplan,icham) - 2.*fgkRpadW)  
                      / ((Float_t) nMCMrow); 
      Int_t   nMCMcol = 8; 
      Float_t xSize   = (GetChamberWidth(iplan) - 2.* fgkCpadW) 
	              / ((Float_t) nMCMcol); 
      sprintf(cTagV,"UU%02d",iDet); 
      for (Int_t iMCMrow = 0; iMCMrow < nMCMrow; iMCMrow++) { 
        for (Int_t iMCMcol = 0; iMCMcol < nMCMcol; iMCMcol++) { 
          xpos   = (0.5 + iMCMcol) * xSize + 1.0  
                 - fCWidth[iplan]/2.; 
          ypos   = (0.5 + iMCMrow) * ySize + 1.0  
                 - fClength[iplan][icham]/2. + fgkHspace/2.; 
          zpos   = -0.4 + 0.742/2.; 
          par[0] = 0.0; 
          par[1] = 0.2/2.; // Thickness of the power lines 
          par[2] = fCWidth[iplan]/2.; 
          gMC->Gspos("UMCM",iCopy+iMCMrow*10+iMCMcol,cTagV 
                           ,xpos,ypos,zpos,0,"ONLY"); 
	} 
      } 
 
    } 
  } 
} 
*/ 
 
//_____________________________________________________________________________ 
 /*
void IlcDCHgeometry::GroupChamber(Int_t iplan, Int_t icham, Int_t *idtmed) 
{ 
  // 
  // Group volumes UA, UD, UF, UU in a single chamber (Air) 
  // UA, UD, UF, UU are boxes 
  // UT will be a box 
  // 
  // ... for the moment there are no services (UU) for the upper plane ! 
  // 
  const Int_t kNparCha = 3; 
 
  Int_t iDet = GetDetectorSec(iplan,icham); 
 
  Float_t xyzMin[3]; 
  Float_t xyzMax[3]; 
  Float_t xyzOrig[3]; 
  Float_t xyzBoxd[3]; 
 
  Char_t  cTagV[5]; 
  Char_t  cTagM[5]; 
 
  for (Int_t i = 0; i < 3; i++) { 
    xyzMin[i] = +9999; xyzMax[i] = -9999; 
  } 
 
  for (Int_t i = 0; i < 3; i++) { 
 
    xyzMin[i] = TMath::Min(xyzMin[i],fChamberUAorig[iDet][i]-fChamberUAboxd[iDet][i]); 
    xyzMax[i] = TMath::Max(xyzMax[i],fChamberUAorig[iDet][i]+fChamberUAboxd[iDet][i]); 
 
    xyzMin[i] = TMath::Min(xyzMin[i],fChamberUDorig[iDet][i]-fChamberUDboxd[iDet][i]); 
    xyzMax[i] = TMath::Max(xyzMax[i],fChamberUDorig[iDet][i]+fChamberUDboxd[iDet][i]); 
 
    xyzMin[i] = TMath::Min(xyzMin[i],fChamberUForig[iDet][i]-fChamberUFboxd[iDet][i]); 
    xyzMax[i] = TMath::Max(xyzMax[i],fChamberUForig[iDet][i]+fChamberUFboxd[iDet][i]); 
 
    // CBL 
    //if (iplan < (kNplan-1)) { 
      xyzMin[i] = TMath::Min(xyzMin[i],fChamberUUorig[iDet][i]-fChamberUUboxd[iDet][i]); 
      xyzMax[i] = TMath::Max(xyzMax[i],fChamberUUorig[iDet][i]+fChamberUUboxd[iDet][i]); 
      //} 
 
    xyzOrig[i] = 0.5*(xyzMax[i]+xyzMin[i]); 
    xyzBoxd[i] = 0.5*(xyzMax[i]-xyzMin[i]); 
 
  } 
   
  sprintf(cTagM,"UT%02d",iDet); 
 
  gMC->Gsvolu(cTagM,"BOX ",idtmed[1302-1],xyzBoxd,kNparCha); 
 
  sprintf(cTagV,"UA%02d",iDet); 
  gMC->Gspos(cTagV,1,cTagM, 
	     fChamberUAorig[iDet][0]-xyzOrig[0], 
	     fChamberUAorig[iDet][1]-xyzOrig[1], 
	     fChamberUAorig[iDet][2]-xyzOrig[2], 
	     0,"ONLY"); 
 
  sprintf(cTagV,"UD%02d",iDet); 
  gMC->Gspos(cTagV,1,cTagM, 
	     fChamberUDorig[iDet][0]-xyzOrig[0], 
	     fChamberUDorig[iDet][1]-xyzOrig[1], 
	     fChamberUDorig[iDet][2]-xyzOrig[2], 
	     0,"ONLY"); 
 
  sprintf(cTagV,"UF%02d",iDet); 
  gMC->Gspos(cTagV,1,cTagM, 
	     fChamberUForig[iDet][0]-xyzOrig[0], 
	     fChamberUForig[iDet][1]-xyzOrig[1], 
	     fChamberUForig[iDet][2]-xyzOrig[2], 
	     0,"ONLY"); 
   
  // CBL 
  //if (iplan < (kNplan-1)) { 
    sprintf(cTagV,"UU%02d",iDet); 
    gMC->Gspos(cTagV,1,cTagM, 
	       fChamberUUorig[iDet][0]-xyzOrig[0], 
	       fChamberUUorig[iDet][1]-xyzOrig[1], 
	       fChamberUUorig[iDet][2]-xyzOrig[2], 
	       0,"ONLY"); 
 
    // } 
 
  sprintf(cTagV,"UT%02d",iDet); 
  gMC->Gspos(cTagV,1,"UTI1",xyzOrig[0],xyzOrig[1],xyzOrig[2],0,"ONLY"); 
} 
 */
 
//_____________________________________________________________________________ 
/* Bool_t IlcDCHgeometry::Local2Global(Int_t idet, Double_t *local 
                                   , Double_t *global) const 
{ 
  // 
  // Converts local pad-coordinates (row,col,time) into  
  // global ILC reference frame coordinates (x,y,z) 
  // 
 
  Int_t icham = GetChamber(idet);    // Chamber info (0-4) 
  Int_t isect = GetSector(idet);     // Sector info  (0-17) 
  Int_t iplan = GetPlane(idet);      // Plane info   (0-5) 
 
  return Local2Global(iplan,icham,isect,local,global); 
 
} 
  
//_____________________________________________________________________________ 
Bool_t IlcDCHgeometry::Local2Global(Int_t iplan, Int_t icham, Int_t isect 
                                  , Double_t *local, Double_t *global) const 
{ 
  // 
  // Converts local pad-coordinates (row,col,time) into  
  // global ILC reference frame coordinates (x,y,z) 
  // 
 
  IlcDCHCommonParam* commonParam = IlcDCHCommonParam::Instance(); 
  if (!commonParam) 
    return kFALSE; 
 
  IlcDCHcalibDB* calibration = IlcDCHcalibDB::Instance(); 
  if (!calibration) 
    return kFALSE;   
   
  IlcDCHpadPlane *padPlane = commonParam->GetPadPlane(iplan,icham); 
 
  // calculate (x,y,z) position in rotated chamber 
  Int_t    row       = ((Int_t) local[0]); 
  Int_t    col       = ((Int_t) local[1]); 
  Float_t  timeSlice = local[2] + 0.5; 
  Float_t  time0     = GetTime0(iplan); 
 
  Int_t idet = GetDetector(iplan, icham, isect); 
 
  Double_t  rot[3]; 
  rot[0] = time0 - (timeSlice - calibration->GetT0(idet, col, row)) 
      * calibration->GetVdrift(idet, col, row)/calibration->GetSamplingFrequency(); 
  rot[1] = padPlane->GetColPos(col) - 0.5 * padPlane->GetColSize(col); 
  rot[2] = padPlane->GetRowPos(row) - 0.5 * padPlane->GetRowSize(row); 
 
  // Rotate back to original position 
  return RotateBack(idet,rot,global); 
 
} 
 */ 
//_____________________________________________________________________________ 
/* Bool_t IlcDCHgeometry::Global2Local(Int_t mode, Double_t *local, Double_t *global 
                                  , Int_t* index) const 
{ 
  // 
  // Converts local pad-coordinates (row,col,time) into  
  // global ILC reference frame coordinates (x,y,z) 
  // 
  // index[0] = plane number 
  // index[1] = chamber number 
  // index[2] = sector number 
  // 
  // mode=0  - local coordinate in y, z,             x - rotated global    
  // mode=2  - local coordinate in pad, and pad row, x - rotated global 
  // 
 
  Int_t idet = GetDetector(index[0],index[1],index[2]); // Detector number 
  RotateBack(idet,global,local); 
  if (mode == 0) return kTRUE; 
 
  return kTRUE; 
 
} 
 */ 
//_____________________________________________________________________________ 
/* Bool_t IlcDCHgeometry::Global2Detector(Double_t global[3], Int_t index[3]) 
{ 
  //   
  //  Find detector for given global point - Ideal geometry  
  //   
  // 
  // input    = global position 
  // output   = index 
  // index[0] = plane number 
  // index[1] = chamber number 
  // index[2] = sector number 
  // 
 
  // 
  // Find sector 
  // 
  Float_t fi = TMath::ATan2(global[1],global[0]); 
  if (fi < 0) { 
    fi += 2*TMath::Pi(); 
  } 
  index[2] = fgkNsect - 1 - TMath::Nint((fi - GetAlpha()/2.)/GetAlpha()); 
 
  // 
  // Find plane 
  // 
  Float_t locx = global[0] * fRotA11[index[2]] + global[1] * fRotA12[index[2]];   
  index[0] = 0; 
  Float_t max = locx - GetTime0(0); 
  for (Int_t iplane=1; iplane<fgkNplan;iplane++){ 
    Float_t dist = TMath::Abs(locx - GetTime0(iplane)); 
    if (dist < max){ 
      index[0] = iplane; 
      max = dist; 
    } 
  } 
 
  // 
  // Find chamber 
  // 
  if (TMath::Abs(global[2]) < 0.5*GetChamberLength(index[0],2)){ 
    index[1]=2; 
  } 
  else{ 
    Double_t localZ = global[2]; 
    if (global[2] > 0){ 
      localZ -= 0.5*(GetChamberLength(index[0],2)+GetChamberLength(index[0],1)); 
      index[1] = (TMath::Abs(localZ) < 0.5*GetChamberLength(index[0],3)) ? 1:0; 
    } 
    else{ 
      localZ += 0.5*(GetChamberLength(index[0],2)+GetChamberLength(index[0],3)); 
      index[1] = (TMath::Abs(localZ) < 0.5*GetChamberLength(index[0],1)) ? 3:4; 
    } 
  }   
 
  return kTRUE; 
 
} 
 */ 
//_____________________________________________________________________________ 
/* Bool_t IlcDCHgeometry::Rotate(Int_t d, Double_t *pos, Double_t *rot) const 
{ 
  // 
  // Rotates all chambers in the position of sector 0 and transforms 
  // the coordinates in the ILC restframe <pos> into the  
  // corresponding local frame <rot>. 
  // 
 
  Int_t sector = GetSector(d); 
  rot[0] =  pos[0] * fRotA11[sector] + pos[1] * fRotA12[sector]; 
  rot[1] =  -pos[0] * fRotA21[sector] + pos[1] * fRotA22[sector]; 
  rot[2] =  pos[2]; 
 
  return kTRUE; 
 
} 
 
//_____________________________________________________________________________ 
Bool_t IlcDCHgeometry::RotateBack(Int_t d, Double_t *rot, Double_t *pos) const 
{ 
  // 
  // Rotates a chambers from the position of sector 0 into its 
  // original position and transforms the corresponding local frame  
  // coordinates <rot> into the coordinates of the ILC restframe <pos>. 
  // 
 
  Int_t sector = GetSector(d); 
  pos[0] =  rot[0] * fRotB11[sector] + rot[1] * fRotB12[sector]; 
  pos[1] =  -rot[0] * fRotB21[sector] + rot[1] * fRotB22[sector]; 
  pos[2] =  rot[2]; 
 
  return kTRUE; 
 
} 
 
//_____________________________________________________________________________ 
Int_t IlcDCHgeometry::GetDetectorSec(Int_t p, Int_t c) 
{ 
  // 
  // Convert plane / chamber into detector number for one single sector 
  // 
 
  return (p + c * fgkNplan); 
 
} 
 
//_____________________________________________________________________________ 
Int_t IlcDCHgeometry::GetDetector(Int_t p, Int_t c, Int_t s) 
{ 
  // 
  // Convert plane / chamber / sector into detector number 
  // 
 
  return (p + c * fgkNplan + s * fgkNplan * fgkNcham); 
 
} 
 
//_____________________________________________________________________________ 
Int_t IlcDCHgeometry::GetPlane(Int_t d) const 
{ 
  // 
  // Reconstruct the plane number from the detector number 
  // 
 
  return ((Int_t) (d % fgkNplan)); 
 
} 
 
//_____________________________________________________________________________ 
Int_t IlcDCHgeometry::GetChamber(Int_t d) const 
{ 
  // 
  // Reconstruct the chamber number from the detector number 
  // 
 
  return ((Int_t) (d % (fgkNplan * fgkNcham)) / fgkNplan); 
 
} 
 
//_____________________________________________________________________________ 
Int_t IlcDCHgeometry::GetSector(Int_t d) const 
{ 
  // 
  // Reconstruct the sector number from the detector number 
  // 
 
  
  return ((Int_t) (d / (fgkNplan * fgkNcham)));


} 
 */ /*
//_____________________________________________________________________________ 
IlcDCHgeometry* IlcDCHgeometry::GetGeometry(IlcRunLoader* runLoader) 
{ 
  // 
  // load the geometry from the gilc file 
  // 
 
  if (!runLoader) runLoader = IlcRunLoader::GetRunLoader(); 
  if (!runLoader) { 
    ::Error("IlcDCHgeometry::GetGeometry", "No run loader"); 
    return NULL; 
  } 
 
  TDirectory* saveDir = gDirectory; 
  runLoader->CdGAFile(); 
  //cout<<gIlc->GetRunLoader->CdGAFile()<<endl;
  // Try from the gilc.root file 
  IlcDCHgeometry* geom =(IlcDCHgeometry*) gDirectory->Get("DCHgeometry"); 
  if (!geom) { 
    // It is not in the file, try to get it from gIlc,  
    // which corresponds to the run loader  
    IlcDCH * dch = (IlcDCH*)runLoader->GetIlcRun()->GetDetector("DCH"); 
    if(dch)
      geom = dch->GetGeometry(); 
  } 
  if (!geom) ::Error("IlcDCHgeometry::GetGeometry", "Geometry not found"); 
 
  saveDir->cd(); 
  return geom; 
 
} 
    */
//_____________________________________________________________________________ 
/* Bool_t IlcDCHgeometry::ReadGeoMatrices() 
{ 
  // 
  // Read geo matrices from current gGeoManager for each DCH sector 
  // 
 
  if (!gGeoManager) return kFALSE; 
  fMatrixArray = new TObjArray(kNdet);  
  fMatrixCorrectionArray = new TObjArray(kNdet); 
  fMatrixGeo   = new TObjArray(kNdet); 
  IlcAlignObjAngles o; 
 
  for (Int_t iLayer = IlcAlignObj::kTRD1; iLayer <= IlcAlignObj::kTRD6; iLayer++) { 
    for (Int_t iModule = 0; iModule < IlcAlignObj::LayerSize(iLayer); iModule++) { 
      UShort_t volid = IlcAlignObj::LayerToVolUID(iLayer,iModule); 
      const char *path = IlcAlignObj::GetVolPath(volid); 
      if (!gGeoManager->cd(path)) return kFALSE;       
      TGeoHMatrix* m = gGeoManager->GetCurrentMatrix(); 
      Int_t     iLayerDCH    = iLayer-IlcAlignObj::kTRD1; 
      Int_t     isector      = Nsect()-1-(iModule/Ncham()); 
      Int_t     ichamber     = Ncham()-1-(iModule%Ncham()); 
      Int_t     lid          = GetDetector(iLayerDCH,ichamber,isector);     
 
      // 
      // local geo system z-x-y  to x-y--z  
      // 
      fMatrixGeo->AddAt(new TGeoHMatrix(*m),lid); 
       
      TGeoRotation mchange;  
      mchange.RotateY(90); mchange.RotateX(90); 
 
      TGeoHMatrix gMatrix(mchange.Inverse()); 
      gMatrix.MultiplyLeft(m); 
      fMatrixArray->AddAt(new TGeoHMatrix(gMatrix),lid);  
 
      // 
      //  Cluster transformation matrix 
      // 
      TGeoHMatrix  rotMatrix(mchange.Inverse()); 
      rotMatrix.MultiplyLeft(m); 
      Double_t sectorAngle = 20.*(isector%18)+10; 
      TGeoHMatrix  rotSector; 
      rotSector.RotateZ(sectorAngle); 
      rotMatrix.MultiplyLeft(&rotSector);       
 
      fMatrixCorrectionArray->AddAt(new TGeoHMatrix(rotMatrix),lid);        
 
    }     
  } 
 
  return kTRUE; 
 
} 
 */
