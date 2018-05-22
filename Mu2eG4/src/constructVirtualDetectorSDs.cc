//
// Free function to create the make the virtual detectors sensitive
//
// constructVirtualDetectorSDs.cc
// Author: Lisa Goodenough
// 2018/03/12
//

#include "Mu2eG4/inc/constructVirtualDetectorSDs.hh"
#include "Mu2eG4/inc/Mu2eSensitiveDetector.hh"
#include "G4LogicalVolumeStore.hh"

using namespace std;

namespace mu2e {

  // Construct the virtual detectors

  void constructVirtualDetectorSDs(SimpleConfig const & _config,
                                   Mu2eSensitiveDetector* vdSD){
      
      G4LogicalVolumeStore* store = G4LogicalVolumeStore::GetInstance();
      
      for(G4LogicalVolumeStore::iterator pos=store->begin(); pos!=store->end(); pos++){
          
          G4String LVname = (*pos)->GetName();
                    
          //from constructExtMonFNALBoxVirtualDetectors
          //EMFBoxFront, EMFBoxSW, EMFBoxBottom, EMFBoxBack, EMFBoxNE, EMFBoxTop
          if (LVname.find("VirtualDetector_EMFBox") != std::string::npos) {
              store->GetVolume(LVname)->SetSensitiveDetector(vdSD);
          }
          
          //from constructExtMonFNALVirtualDetectors
          //EMFC2Entrance, EMFC2Exit
          if (LVname.find("VirtualDetector_EMFC2") != std::string::npos) {
              store->GetVolume(LVname)->SetSensitiveDetector(vdSD);
          }
          
          //from constructExtMonFNALPlaneStack
          //EMFDetectorUpEntrance, EMFDetectorUpExit, EMFDetectorDnEntrance, EMFDetectorDnExit
          if (LVname.find("VirtualDetector_EMFDetector") != std::string::npos) {
              store->GetVolume(LVname)->SetSensitiveDetector(vdSD);
          }
          
          //***************************************************************
          //***************************************************************
          //everything below here was taken from constructVirtualDetectors()
          //***************************************************************
          //***************************************************************

          
          //Coll1_In, Coll1_Out and Coll1_pBarCollar_In, Coll1_pBarCollar_Out
          if (LVname.find("VirtualDetector_Coll1_") != std::string::npos) {
              store->GetVolume(LVname)->SetSensitiveDetector(vdSD);
          }
          
          if (LVname.find("VirtualDetector_TS2_Bend") != std::string::npos) {
              store->GetVolume(LVname)->SetSensitiveDetector(vdSD);
          }
          
          if (LVname.find("VirtualDetector_TS4_Bend") != std::string::npos) {
              store->GetVolume(LVname)->SetSensitiveDetector(vdSD);
          }
          
          //Coll31_In,Coll31_Out, Coll32_In, Coll32_Out
          if (LVname.find("VirtualDetector_Coll3") != std::string::npos) {
              store->GetVolume(LVname)->SetSensitiveDetector(vdSD);
          }
          
          //Coll5_In, Coll5_Out, Coll5_OutSurf
          if (LVname.find("VirtualDetector_Coll5_") != std::string::npos) {
              store->GetVolume(LVname)->SetSensitiveDetector(vdSD);
          }
          
          if (LVname.find("VirtualDetector_STMUpstream") != std::string::npos) {
              store->GetVolume(LVname)->SetSensitiveDetector(vdSD);
          }
          
          // ST_In, ST_Out
          if (LVname.find("VirtualDetector_ST_") != std::string::npos) {
              store->GetVolume(LVname)->SetSensitiveDetector(vdSD);
          }
          
          //TT_Mid, TT_MidInner, TT_FrontHollow, TT_FrontPA, TT_Back, TT_OutSurf, TT_InSurf
          if (LVname.find("VirtualDetector_TT_") != std::string::npos) {
              store->GetVolume(LVname)->SetSensitiveDetector(vdSD);
          }
          
          //IT_VD_EndCap_Front, IT_VD_EndCap_Back, IT_VD_InSurf
          if (LVname.find("VirtualDetector_IT_VD_") != std::string::npos) {
              store->GetVolume(LVname)->SetSensitiveDetector(vdSD);
          }

          if (LVname.find("VirtualDetector_PS_FrontExit") != std::string::npos) {
              store->GetVolume(LVname)->SetSensitiveDetector(vdSD);
          }
          
          //PT_Front, PT_Back
          if (LVname.find("VirtualDetector_PT_") != std::string::npos) {
              store->GetVolume(LVname)->SetSensitiveDetector(vdSD);
          }
          
          if (LVname.find("VirtualDetector_ProtonBeamDumpCoreFace") != std::string::npos) {
              store->GetVolume(LVname)->SetSensitiveDetector(vdSD);
          }
          
          if (LVname.find("VirtualDetector_DSNeutronShieldExit") != std::string::npos) {
              store->GetVolume(LVname)->SetSensitiveDetector(vdSD);
          }
          
          //EMC_Disk_0_SurfIn, EMC_Disk_0_SurfOut, EMC_Disk_1_SurfIn, EMC_Disk_1_SurfOut
          //EMC_Disk_0_EdgeIn, EMC_Disk_0_EdgeOut, EMC_Disk_1_EdgeIn, EMC_Disk_1_EdgeOut
          if (LVname.find("VirtualDetector_EMC_Disk_") != std::string::npos) {
              store->GetVolume(LVname)->SetSensitiveDetector(vdSD);
          }
          
          //EMC_FEB_0_SurfIn, EMC_FEB_0_SurfOut, EMC_FEB_1_SurfIn, EMC_FEB_1_SurfOut
          //EMC_FEB_0_EdgeIn, EMC_FEB_0_EdgeOut, EMC_FEB_1_EdgeIn, EMC_FEB_1_EdgeOut
          if (LVname.find("VirtualDetector_EMC_FEB_") != std::string::npos) {
              store->GetVolume(LVname)->SetSensitiveDetector(vdSD);
          }
          
          //STM_UpStr, STM_MagDnStr, STM_CollDnStr, STM_Det1UpStr, STM_Det2UpStr
          //STM_FieldOfViewCollDnStr, STM_SpotSizeCollUpStr
          if (LVname.find("VirtualDetector_STM_") != std::string::npos) {
              store->GetVolume(LVname)->SetSensitiveDetector(vdSD);
          }
          
          //PSPbarIn, PSPbarOut
          if (LVname.find("VirtualDetector_PSPbar") != std::string::npos) {
              store->GetVolume(LVname)->SetSensitiveDetector(vdSD);
          }
          
          //CRV_R, CRV_L, CRV_T, CRV_D, CRV_U
          if (LVname.find("VirtualDetector_CRV_") != std::string::npos) {
              store->GetVolume(LVname)->SetSensitiveDetector(vdSD);
          }
   
      }//for G4LogicalVolumeStore::iterator
      
    
  } // constructVirtualDetectorSDs()
}
