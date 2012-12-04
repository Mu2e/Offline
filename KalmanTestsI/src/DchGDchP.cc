#include "KalmanTestsI/inc/DchGDchP.hh"

#include "art/Framework/Services/Registry/ServiceHandle.h"
//#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "ITrackerGeom/inc/ITracker.hh"
#include "MatEnv/RecoMatFactory.hh"
#include "MatEnv/MatDBInfo.hh"
#include "DetectorModel/DetMaterial.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "TMath.h"
#include "G4String.hh"
#include "G4Material.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"

#include "Mu2eG4/inc/WorldMaker.hh"

#include "G4Helper/inc/VolumeInfo.hh"

#include <cmath>
using namespace mu2e;

std::map<std::string,std::string> matnamemap;
bool firstSub;

DchGDchP::DchGDchP(std::string materialfile, bool doDetailedWrSupp, bool useSingleCellElemDescr, bool gas_wire_extWght):
  DchGDch(nLayer)
{
  art::ServiceHandle<GeometryService> geom;
  GeomHandle<ITracker> itr;
  CellGeometryHandle *itwp = itr->getCellGeometryHandle();

  SimpleConfig const& config  = geom->config();
  doDetailedWrSupp = doDetailedWrSupp && config.getBool("itracker.detailedWireSupport",false);

  bool automaterial=true;
  if(materialfile!="") automaterial=false;
  std::cout<<"materialdb "<<materialfile<<" auto re-build "<<automaterial<<std::endl;

  double rz[2][3];
  double inraddensity(0.0), inradtotweigth(0.0), inradradlen(0.0);
  double innerGasGuardWLayer(0.0);
  double insideEndcapMaxR(0.0);
  double sWireRad(0.0), fWireRad(0.0);
  std::map<std::string,double> materialAreaInWll;

  boost::shared_ptr<Wall> rEP, fEP;
  std::multimap<Wall::Walltype,boost::shared_ptr<Wall> >::iterator it;
  for(it=itr->getWalls()->begin();it!=itr->getWalls()->end();it++){
    if((*it).first==Wall::inner||(*it).first==Wall::outer){
      int ind=((*it).first==Wall::inner)?0:1;
      std::cout<<"barrel radius "<<(*it).second->getRmin()<<" "<<(*it).second->getRmax()
	       <<" "<<(*it).second->getDz()<<std::endl;
      rz[ind][0]=(*it).second->getRmin();
      rz[ind][1]=(*it).second->getRmax();
      rz[ind][2]=(*it).second->getDz()*2.0;

      if((*it).first==Wall::inner){
              const std::vector<double>* th = (*it).second->getThicknesses().get();
              const std::vector<std::string>* mats = (*it).second->getMaterialsName().get();
              double thickness=0, tmpDensFact;
              for(unsigned int i=0;i<th->size();i++){
                      G4String matname=(*mats)[i];
                      if (matname.contains("ITGas")) {
                              innerGasGuardWLayer+=(*th)[i];
                              rz[ind][1]-=(*th)[i];
                      } else if (automaterial) {
                              G4Material* mat=mu2e::findMaterialOrThrow(matname);
                              std::cout<<"sum material "<<mat->GetName()<<" thick,mm= "<<(*th)[i]
                                       <<" density,g/cm^3 = "<<mat->GetDensity()/(g/cm/cm/cm)
                                       <<" radlength,g/cm^2 = "<<mat->GetRadlen()/cm*(mat->GetDensity()/(g/cm/cm/cm))<<std::endl;
                              tmpDensFact = (*th)[i]*mat->GetDensity()/(g/cm/cm/cm);
                              inradtotweigth+=tmpDensFact;
                              thickness+=(*th)[i];
                              inradradlen += (*th)[i]/mat->GetRadlen();
                              const G4double* elemfract=mat->GetFractionVector();
                              for(unsigned int i=0;i<mat->GetNumberOfElements();i++){
                                      const G4Element* elem=mat->GetElement(i);
                                      materialAreaInWll[elem->GetName()]+=elemfract[i]*tmpDensFact;
                              }
                      }
              }
              inraddensity=inradtotweigth/thickness;
              inradradlen = thickness/inradradlen/cm*inraddensity;
              std::cout<<"average density = "<<inraddensity
                              <<" g/cm^3, at thickness = "<<thickness<<" mm"<<" average radlength,g/cm^2 = "<< inradradlen <<std::endl;
      } else {
              insideEndcapMaxR = config.getDouble("itracker.rOut",0.0)-(*it).second->getTotalThickness();
      }
    }
    if ((*it).first==Wall::endcap) {
            if ((*it).second->getPos().zz()<0) {
                    rEP = (*it).second;
            } else {
                    fEP = (*it).second;
            }
    }
  }

  matnamemap["G4_Mo"]="Molybdenum";
  matnamemap["G4_Ag"]="Silver";
  matnamemap["G4_Al"]="Aluminum";
  matnamemap["G4_W"]="Tungsten";
  matnamemap["G4_Au"]="Gold";
  matnamemap["G4_N"]="Nitrogen";
  matnamemap["G4_Cu"]="IT-Copper";
  matnamemap["G4_H"]="Hydrogen";
  matnamemap["G4_C"]="Carbon";
  matnamemap["G4_O"]="Oxygen";
  matnamemap["H"]="Hydrogen";
  matnamemap["C"]="Carbon";
  matnamemap["F"]="Fluorine";
  matnamemap["He"]="Helium";
  matnamemap["O"]="Oxygen";
  matnamemap["N"]="Nitrogen";
  matnamemap["Cu"]="Copper";
  matnamemap["WAGVacuum"]="Hydrogen";
  matnamemap["G4_KAPTON"]="Kapton";
  matnamemap["G4_POLYCARBONATE"]="ITPolycarbonateAuto";
  matnamemap["CarbonFiber"]="ITCFibAuto";//"CFiber";
  matnamemap["G10_FR4"]="IT-G10-FR4";//"IT-G10_FR4Auto";  //additional element needed to describe G10 in details  Si Ca Al Mg B Ti Na K Fe
  matnamemap["KptFoam_030"]="ITKptFoam_030";
  matnamemap["MatAluminum"]="IT-Aluminum";
  matnamemap["ITGasVacuum"]="vacuum";//"Hydrogen";


  _name="Mu2e ITracker";    
  initWires();
  _version=999;
  _nLayers=0;

  double activeAreaLength = rz[0][2];
  std::cout<<"gas leng "<<activeAreaLength<<" "<<rz[0][2]<<" "<<rEP->getTotalThickness()<<" "<<fEP->getTotalThickness()<<std::endl;

  if (useSingleCellElemDescr) {
          _GasCyl=DchSimpleCyl(rz[0][1],rz[0][1]+0.1,activeAreaLength);
  } else {
          _GasCyl=DchSimpleCyl(rz[0][1]+innerGasGuardWLayer,rz[1][0],activeAreaLength);
  }
  _GasPos[0]=0;_GasPos[1]=0;_GasPos[2]=0;
  _InnerCyl=DchSimpleCyl(rz[0][0],rz[0][1],activeAreaLength);
  _ICPos[0]=0;_ICPos[1]=0;_ICPos[2]=0;
  _InnerWGrdCyl = DchSimpleCyl(rz[0][1],rz[0][1]+innerGasGuardWLayer,activeAreaLength);
  _IWGrdPos[0]=0;_IWGrdPos[1]=0;_IWGrdPos[2]=0;
  _OuterCyl=DchSimpleCyl(rz[1][0],rz[1][1],rz[1][2]);
  _OCPos[0]=0;_OCPos[1]=0;_OCPos[2]=0;

  //double EPtotLength = itr->maxEndCapDim()-rz[0][2]/2.0;
  //double zEPpos = itr->maxEndCapDim()-0.5*EPtotLength;
  _RearCyl=DchPhiSegmCyl(rz[0][0],rz[1][0],rEP->getTotalThickness());
  _REPPos[0]=0;_REPPos[1]=0;_REPPos[2]=rEP->getPos().dz();
  //_REPMaterialSubs.push_back("");
  //double tmpEPpos = -itr->maxEndCapDim();
  double tmpEPpos = -activeAreaLength/2.0;
  double shellThick;
  double oringWidth  = config.getDouble("itracker.oringWidth",0.0);
  double oringHeight = config.getDouble("itracker.oringHeight",0.0);
  double tmpEPInWlShlMaxR = rz[0][0];
  double tmpEPInWlShlLeng = itr->maxEndCapDim()+tmpEPpos-oringWidth; //maxEndCapDim - innerWall Dz
  double tmpEPInWlShlZpos = tmpEPpos-oringWidth-tmpEPInWlShlLeng/2.0; //_REPPos[2];
  double cableFractionOnREP = config.getDouble("itracker.cableFractionOnREP",0.0);
  if (cableFractionOnREP>1.0) { cableFractionOnREP=1.0; }
  bool solidWall = false;
  firstSub=true;

  int iShl=0;
  double skipEPshlLeng(0.0);
  if (doDetailedWrSupp) {
          skipEPshlLeng = rEP->getThicknesses()->at(0)+rEP->getThicknesses()->at(1);
          iShl=2;
          makeDetailedWrSupp(_RearCylSubs,
                          _REPPosSubs,
                          _REPMaterialSubs,
                          tmpEPpos, skipEPshlLeng, rz[0][1], insideEndcapMaxR, false);
          tmpEPpos-=skipEPshlLeng;
  }
//  for (int iShl=rEP->getNShells()-1; iShl>-1; --iShl){
  for ( ; iShl<rEP->getNShells(); ++iShl){
          solidWall = false;

          shellThick = rEP->getThicknesses()->at(iShl);
          tmpEPpos-=0.5*shellThick;
          if (firstSub) {
                  _RearCylSubs.back() = DchPhiSegmCyl(rz[0][0],rz[1][0],shellThick);
                  _REPPosSubs.back()[2]=tmpEPpos;
                  if ( rEP->getMaterialsName()->at(iShl).find("ITGas")!=std::string::npos ) {
                          _REPMaterialSubs.back()="ITgasAuto";
                  } else {
                          _REPMaterialSubs.back()=matnamemap[rEP->getMaterialsName()->at(iShl)];
                          solidWall = true;
                  }
                  firstSub=false;
          } else {
                  if (iShl==2 && cableFractionOnREP>0) {
                          makeInEPCables(_RearCylSubs,
                                      _REPPosSubs,
                                      _REPMaterialSubs,
                                      tmpEPpos+0.5*shellThick, shellThick, rz[0][0], rz[1][0], true);
                  } else {
                          _RearCylSubs.push_back( DchPhiSegmCyl(rz[0][0],rz[1][0],shellThick) );
                          _REPPosSubs.push_back(new double[3]);
                          _REPPosSubs.back()[0]=_REPPosSubs.back()[1]=0.0;
                          _REPPosSubs.back()[2]=tmpEPpos;
                          if ( rEP->getMaterialsName()->at(iShl).find("ITGas")!=std::string::npos ) {
                                  _REPMaterialSubs.push_back("ITgasAuto");
                          } else {
                                  _REPMaterialSubs.push_back(matnamemap[rEP->getMaterialsName()->at(iShl)]);
                                  solidWall = true;
                          }
                  }
          }
          tmpEPpos-=0.5*shellThick;

          if (solidWall) {
                  //closing the inner part of the endacp with the solid materials
                  _RearCylSubs.push_back( DchPhiSegmCyl(tmpEPInWlShlMaxR-shellThick,tmpEPInWlShlMaxR,tmpEPInWlShlLeng) );
                  _REPPosSubs.push_back(new double[3]);
                  _REPPosSubs.back()[0]=_REPPosSubs.back()[1]=0.0;
                  _REPPosSubs.back()[2]=tmpEPInWlShlZpos;
                  _REPMaterialSubs.push_back(matnamemap[rEP->getMaterialsName()->at(iShl)]);
                  std::cout<<"--- inner REP "<< " --- IR "<<_RearCylSubs.back().getInnerRadius()<<" OR "<<_RearCylSubs.back().getOuterRadius()<<" Leng "<<_RearCylSubs.back().getLength()<<" center zpos "<<_REPPosSubs.back()[2]<<" material "<<_REPMaterialSubs.back()<<std::endl;
                  tmpEPInWlShlMaxR-=shellThick;
          }
  }
  if (oringWidth>0.0) {
          _RearCylSubs.push_back( DchPhiSegmCyl(rz[0][0]-oringHeight,rz[0][0],oringWidth) );
          _REPPosSubs.push_back(new double[3]);
          _REPPosSubs.back()[0]=_REPPosSubs.back()[1]=0.0;
          _REPPosSubs.back()[2]=-(activeAreaLength+oringWidth)/2.0;
          _REPMaterialSubs.push_back(matnamemap[config.getString("itracker.oringMaterial","ITGasVacuum")]);
          std::cout<<"--- oring REP "<< " --- IR "<<_RearCylSubs.back().getInnerRadius()<<" OR "<<_RearCylSubs.back().getOuterRadius()<<" Leng "<<_RearCylSubs.back().getLength()<<" center zpos "<<_REPPosSubs.back()[2]<<" material "<<_REPMaterialSubs.back()<<std::endl;
  }


  _ForwCyl=DchPhiSegmCyl(rz[0][0],rz[1][0],fEP->getTotalThickness());
  _FEPPos[0]=0;_FEPPos[1]=0;_FEPPos[2]=fEP->getPos().dz();
  //_FEPMaterialSubs.push_back("");
  //tmpEPpos = itr->maxEndCapDim();
  tmpEPpos = activeAreaLength/2.0;
  tmpEPInWlShlMaxR = rz[0][0];
  tmpEPInWlShlLeng = itr->maxEndCapDim()-tmpEPpos-oringWidth; //maxEndCapDim - innerWall Dz
  tmpEPInWlShlZpos = tmpEPpos+oringWidth+tmpEPInWlShlLeng/2.0; //_FEPPos[2];
  firstSub=true;

  iShl=0;
  if (doDetailedWrSupp) {
          skipEPshlLeng = rEP->getThicknesses()->at(0)+rEP->getThicknesses()->at(1);
          iShl=2;
          makeDetailedWrSupp(_ForwCylSubs,
                          _FEPPosSubs,
                          _FEPMaterialSubs,
                          tmpEPpos, skipEPshlLeng, rz[0][1], insideEndcapMaxR, true);
          tmpEPpos+=skipEPshlLeng;
  }
//  for (int iShl=fEP->getNShells()-1; iShl>-1; --iShl){
  for ( ; iShl<fEP->getNShells(); ++iShl){
          solidWall = false;

          shellThick = fEP->getThicknesses()->at(iShl);
          tmpEPpos+=0.5*shellThick;
          if (firstSub) {
                  _ForwCylSubs.back() = DchPhiSegmCyl(rz[0][0],rz[1][0],shellThick);
                  _FEPPosSubs.back()[2]=tmpEPpos;
                  if ( fEP->getMaterialsName()->at(iShl).find("ITGas")!=std::string::npos ) {
                          _FEPMaterialSubs.back()="ITgasAuto";
                  } else {
                          _FEPMaterialSubs.back()=matnamemap[fEP->getMaterialsName()->at(iShl)];
                          solidWall = true;
                  }
                  firstSub=false;
          } else {
                  if (iShl==2) {
                          makeInEPCables(_ForwCylSubs,
                                      _FEPPosSubs,
                                      _FEPMaterialSubs,
                                      tmpEPpos-0.5*shellThick, shellThick, rz[0][0], rz[1][0], false);
                  } else {
                          _ForwCylSubs.push_back( DchPhiSegmCyl(rz[0][0],rz[1][0],shellThick) );
                          _FEPPosSubs.push_back(new double[3]);
                          _FEPPosSubs.back()[0]=_FEPPosSubs.back()[1]=0.0;
                          _FEPPosSubs.back()[2]=tmpEPpos;
                          if ( fEP->getMaterialsName()->at(iShl).find("ITGas")!=std::string::npos ) {
                                  _FEPMaterialSubs.push_back("ITgasAuto");
                          } else {
                                  _FEPMaterialSubs.push_back(matnamemap[fEP->getMaterialsName()->at(iShl)]);
                                  solidWall = true;
                          }
                  }
          }
          tmpEPpos+=0.5*shellThick;

          if (solidWall) {
                  //closing the inner part of the endacp with the solid materials
                  _ForwCylSubs.push_back( DchPhiSegmCyl(tmpEPInWlShlMaxR-shellThick,tmpEPInWlShlMaxR,tmpEPInWlShlLeng) );
                  _FEPPosSubs.push_back(new double[3]);
                  _FEPPosSubs.back()[0]=_FEPPosSubs.back()[1]=0.0;
                  _FEPPosSubs.back()[2]=tmpEPInWlShlZpos;
                  _FEPMaterialSubs.push_back(matnamemap[rEP->getMaterialsName()->at(iShl)]);
                  std::cout<<"--- inner FEP "<< " --- IR "<<_ForwCylSubs.back().getInnerRadius()<<" OR "<<_ForwCylSubs.back().getOuterRadius()<<" Leng "<<_ForwCylSubs.back().getLength()<<" center zpos "<<_FEPPosSubs.back()[2]<<" material "<<_FEPMaterialSubs.back()<<std::endl;
                  tmpEPInWlShlMaxR-=shellThick;
          }
  }
  if (oringWidth>0.0) {
          _ForwCylSubs.push_back( DchPhiSegmCyl(rz[0][0]-oringHeight,rz[0][0],oringWidth) );
          _FEPPosSubs.push_back(new double[3]);
          _FEPPosSubs.back()[0]=_FEPPosSubs.back()[1]=0.0;
          _FEPPosSubs.back()[2]=(activeAreaLength+oringWidth)/2.0;
          _FEPMaterialSubs.push_back(matnamemap[config.getString("itracker.oringMaterial","ITGasVacuum")]);
          std::cout<<"--- oring FEP "<< " --- IR "<<_ForwCylSubs.back().getInnerRadius()<<" OR "<<_ForwCylSubs.back().getOuterRadius()<<" Leng "<<_ForwCylSubs.back().getLength()<<" center zpos "<<_FEPPosSubs.back()[2]<<" material "<<_FEPMaterialSubs.back()<<std::endl;
  }

  for(int i=0;i<itr->nSuperLayers();i++){
    boost::shared_ptr<ITLayer> ily;
    for (int il=0; il<itr->getSuperLayer(i)->nLayers(); il++) {
       ily = itr->getSuperLayer(i)->getLayer(il);
       if (ily->nCells()!=0){
         //std::cout<<"DchGDchP SL "<<i<<" layer "<<ily->getCell(0)->Id().getLayer()<<" nCells "<<ily->nCells()<<std::endl;
         _sWires[i].nCells=ily->nCells();
         itwp->SelectCell(i,ily->getCell(0)->Id().getLayer(),0);
         Float_t xup[3],xdown[3];
         itwp->WirePosAtEndcap(xup,xdown);
         _sWires[i].rEnd=TMath::Hypot(xdown[0],xdown[1]);
         _sWires[i].rForw=TMath::Hypot(xup[0],xup[1]);
         _sWires[i].phiEnd=TMath::ATan2(xdown[1],xdown[0]);
         _sWires[i].phiForw=TMath::ATan2(xup[1],xup[0]);
         _sWires[i].twist=_sWires[i].phiForw-_sWires[i].phiEnd;
         _sWires[i].stereo=-TMath::Tan(itwp->GetWireEpsilon());
         _sWires[i].deltaPhi=TMath::TwoPi()/_sWires[i].nCells;
         _sWires[i].zR=xdown[2];
         _sWires[i].zFw=xup[2];
         _sWires[i].sag=0;
         _sWires[i].diameter=2.0*itwp->GetITCell()->getWire()->getDetail()->outerRadius();//0.005;
         _sWires[i].volts=2000;

         //         std::cout<<i<<" "<<itwp->GetWireCenter()<<" n="<<itwp->GetWireDirection()<<std::endl;
         std::cout<<i<<" r "<<itwp->GetWireCenter().rho()<<" phi0 "<<itwp->GetWireCenter().phi()
                  <<" dphi "<<_sWires[i].deltaPhi
                  <<" sangle "<<-itwp->GetWireEpsilon()
                  <<" z "<<fabs(xdown[2])<<" nw "<<_sWires[i].nCells<<" n= "<<itwp->GetWireDirection()<<std::endl;
         _cellStruc[i].nofFWires=0;
         _rOffset[i]=0;
         //break;
         ++_nLayers;
       }
    }
  }
  _zoffset=0;
  _zlen=0;
  _cellHeight=0;

  RecoMatFactory* matFactory = RecoMatFactory::getInstance();

  _ICMaterial = "ITInBarrelAuto";
  _GasMaterial = "ITgaswireAuto";
  _IWGrdMaterial = "ITgasAuto";

  if(!automaterial){
    mu2e::ConfigFileLookupPolicy findFile;
    std::string fullPath = findFile(materialfile);
    MatMaterialList* mtrList = new MatMaterialList(fullPath);
    ((MatMtrDictionary*)matFactory->materialDictionary())->FillMtrDict(mtrList);
    std::cout<<"My materials"<<std::endl;
    MatDBInfo* mtdbinfo=new MatDBInfo;
    mtdbinfo->findDetMaterial(_ICMaterial.c_str())->printAll(std::cout);
    mtdbinfo->findDetMaterial(_GasMaterial.c_str())->printAll(std::cout);
    mtdbinfo->findDetMaterial("ITCFibAuto")->printAll(std::cout);
    mtdbinfo->findDetMaterial("ITKptFoam_030")->printAll(std::cout);
    mtdbinfo->findDetMaterial("ITPolycarbonateAuto")->printAll(std::cout);
    mtdbinfo->findDetMaterial("ITFWireAuto")->printAll(std::cout);
    mtdbinfo->findDetMaterial("ITSWireAuto")->printAll(std::cout);
    mtdbinfo->findDetMaterial("ITgasAuto")->printAll(std::cout);
    mtdbinfo->findDetMaterial("IT-G10-FR4")->printAll(std::cout);
    mtdbinfo->findDetMaterial("ITSignalCableAuto")->printAll(std::cout);
    //mtdbinfo->findDetMaterial("ITSCblREPAuto")->printAll(std::cout);
    //mtdbinfo->findDetMaterial("ITSCblFEPAuto")->printAll(std::cout);

    return;
  }

  //create average materials from geometry description 
  std::map<std::string*, MatMaterialObj*, PtrLess>* matobj=matFactory->materialDictionary();

  double zeff = 0., aeff = 0., intlen = -20., refindex = -30.; //auto produced

  std::string *ptrname=new std::string(_ICMaterial.c_str());
  double density = inraddensity; // g/cm^3
  double radlen = inradradlen;      // g/cm^2
  double temperature = 20.;  // deg
  double pressure = 1.;      // atm
  std::string state="solid"; // solid or gas
  int nbrcomp=materialAreaInWll.size(); // + plus by weigth; - by atoms
  std::vector<int> Compflag; // 0 - element, 1 - material
  std::vector<double> Compweight;
  std::vector<std::string> Compname;
  for (std::map<std::string,double>::iterator itm = materialAreaInWll.begin(); itm != materialAreaInWll.end(); ++itm) {
          Compflag.push_back(0);
          Compweight.push_back(itm->second/inradtotweigth);
          std::cout<<"G4 replace "<<itm->first<<" by MatEnv "<<matnamemap[itm->first]<<" with weight "<<Compweight.back()<<std::endl;
          Compname.push_back(matnamemap[itm->first]);
  }

  MatMaterialObj* Obj = new MatMaterialObj(*ptrname, density, zeff,
                                           aeff, nbrcomp, &Compflag, &Compweight,
                                           &Compname, radlen, intlen, refindex,
                                           temperature, pressure, state);
  (*matobj)[ptrname]=Obj;

  int nsenseWires=0,nfieldWires=0;
  WireDetail* senseWire=0,*fieldWire=0;
  for(int i=0;i<itr->nSuperLayers();i++){
    SuperLayer *sl=itr->getSuperLayer(i);
    for(int j=0;j<sl->nLayers();j++){
      boost::shared_ptr<ITLayer> l=sl->getLayer(j);
      nsenseWires+=l->nCells();
      nfieldWires+=l->nFieldWires();
      if(!senseWire&&l->nCells())
	senseWire=l->getCell(0)->getWire()->getDetail().get();
      if(!fieldWire&&l->nFieldWires())
	fieldWire=l->getFWire(0)->getDetail().get();
      
    }
  }


  G4Material* mat=mu2e::findMaterialOrThrow("CarbonFiber");
  std::cout<<"cfib "<<mat->GetName()<<" density,g/cm^3 = "<<mat->GetDensity()/(g/cm/cm/cm)<<std::endl;
  std::string *ptrcfibname=new std::string("ITCFibAuto");
  density = mat->GetDensity()/(g/cm/cm/cm); // g/cm^3
  radlen = -10;      // g/cm^2
  temperature = 20.;  // deg
  pressure = mat->GetPressure()/atmosphere;      // atm
  state="solid"; // solid or gas
  nbrcomp=mat->GetElementVector()->size(); // + plus by weigth; - by atoms
  Compflag.clear();Compweight.clear();Compname.clear();
  for(unsigned int i=0;i<mat->GetNumberOfElements();i++){
          const G4Element* elem=mat->GetElement(i);
          Compflag.push_back(0);
          Compweight.push_back(mat->GetFractionVector()[i]);
          std::cout<<"G4 replace "<<elem->GetName()<<" by MatEnv "<<matnamemap[elem->GetName()]<<std::endl;
          Compname.push_back(matnamemap[elem->GetName()]);
  }
  Obj = new MatMaterialObj(*ptrcfibname, density, zeff,
                           aeff, nbrcomp, &Compflag, &Compweight,
                           &Compname, radlen, intlen, refindex,
                           temperature, pressure, state);
  (*matobj)[ptrcfibname]=Obj;

  mat=mu2e::findMaterialOrThrow("KptFoam_030");
  std::cout<<"kfoam "<<mat->GetName()<<" density,g/cm^3 = "<<mat->GetDensity()/(g/cm/cm/cm)<<std::endl;
  std::string *ptrkfoamname=new std::string("ITKptFoam_030");
  density = mat->GetDensity()/(g/cm/cm/cm); // g/cm^3
  radlen = -10;      // g/cm^2
  temperature = 20.;  // deg
  pressure = mat->GetPressure()/atmosphere;      // atm
  state="solid"; // solid or gas
  nbrcomp=mat->GetElementVector()->size(); // + plus by weigth; - by atoms
  Compflag.clear();Compweight.clear();Compname.clear();
  for(unsigned int i=0;i<mat->GetNumberOfElements();i++){
          const G4Element* elem=mat->GetElement(i);
          Compflag.push_back(0);
          Compweight.push_back(mat->GetFractionVector()[i]);
          std::cout<<"G4 replace "<<elem->GetName()<<" by MatEnv "<<matnamemap[elem->GetName()]<<std::endl;
          Compname.push_back(matnamemap[elem->GetName()]);
  }
  Obj = new MatMaterialObj(*ptrkfoamname, density, zeff,
                           aeff, nbrcomp, &Compflag, &Compweight,
                           &Compname, radlen, intlen, refindex,
                           temperature, pressure, state);
  (*matobj)[ptrkfoamname]=Obj;

  mat=mu2e::findMaterialOrThrow("G4_POLYCARBONATE");
  std::cout<<"ployc "<<mat->GetName()<<" density,g/cm^3 = "<<mat->GetDensity()/(g/cm/cm/cm)<<std::endl;
  std::string *ptrploycname=new std::string("ITPolycarbonateAuto");
  density = mat->GetDensity()/(g/cm/cm/cm); // g/cm^3
  radlen = -10;      // g/cm^2
  temperature = 20.;  // deg
  pressure = mat->GetPressure()/atmosphere;      // atm
  state="solid"; // solid or gas
  nbrcomp=mat->GetElementVector()->size(); // + plus by weigth; - by atoms
  Compflag.clear();Compweight.clear();Compname.clear();
  for(unsigned int i=0;i<mat->GetNumberOfElements();i++){
          const G4Element* elem=mat->GetElement(i);
          Compflag.push_back(0);
          Compweight.push_back(mat->GetFractionVector()[i]);
          std::cout<<"G4 replace "<<elem->GetName()<<" by MatEnv "<<matnamemap[elem->GetName()]<<std::endl;
          Compname.push_back(matnamemap[elem->GetName()]);
  }
  Obj = new MatMaterialObj(*ptrploycname, density, zeff,
                           aeff, nbrcomp, &Compflag, &Compweight,
                           &Compname, radlen, intlen, refindex,
                           temperature, pressure, state);
  (*matobj)[ptrploycname]=Obj;

  /*
  mat=mu2e::findMaterialOrThrow("G10_FR4");
  std::cout<<"g10 "<<mat->GetName()<<" density,g/cm^3 = "<<mat->GetDensity()/(g/cm/cm/cm)<<std::endl;
  std::string *ptrg10name=new std::string("IT-G10_FR4Auto");
  density = mat->GetDensity()/(g/cm/cm/cm); // g/cm^3
  radlen = -10;      // g/cm^2
  temperature = 20.;  // deg
  pressure = mat->GetPressure()/atmosphere;      // atm
  state="solid"; // solid or gas
  nbrcomp=mat->GetElementVector()->size(); // + plus by weigth; - by atoms
  Compflag.clear();Compweight.clear();Compname.clear();
  for(unsigned int i=0;i<mat->GetNumberOfElements();i++){
          const G4Element* elem=mat->GetElement(i);
          Compflag.push_back(0);
          Compweight.push_back(mat->GetFractionVector()[i]);
          std::cout<<"G4 replace "<<elem->GetName()<<" by MatEnv "<<matnamemap[elem->GetName()]<<std::endl;
          Compname.push_back(matnamemap[elem->GetName()]);
  }
  Obj = new MatMaterialObj(*ptrg10name, density, zeff,
                           aeff, nbrcomp, &Compflag, &Compweight,
                           &Compname, radlen, intlen, refindex,
                           temperature, pressure, state);
  (*matobj)[ptrg10name]=Obj;
  */

  std::cout<<"number of sense wires = "<<nsenseWires<<" field wires = "<<nfieldWires<<std::endl;
  std::map<std::string,double> materialArea;

  std::map<std::string,double> materialAreaSW;
  std::vector<std::string> matnameSW=senseWire->materialNames();
  std::vector<double> thick=senseWire->shellsThicknesses();
  double radius=0.0;
  for(unsigned int i=0;i<matnameSW.size();i++){
    std::cout<<"sense material name = "<<matnameSW[i]<<" dr,mm = "<<thick[i]<<std::endl;
    materialArea[matnameSW[i]]+=TMath::Pi()*((thick[i]+radius)*(thick[i]+radius)-radius*radius)*nsenseWires;//*0.05;
    materialAreaSW[matnameSW[i]]+=TMath::Pi()*((thick[i]+radius)*(thick[i]+radius)-radius*radius);
    radius+=thick[i];
  }
  sWireRad = senseWire->outerRadius();

  std::map<std::string,double> materialAreaFW;
  std::vector<std::string> matnameFW=fieldWire->materialNames();
  thick=fieldWire->shellsThicknesses();
  radius=0.0;
  for(unsigned int i=0;i<matnameFW.size();i++){
    std::cout<<"field material name = "<<matnameFW[i]<<" dr,mm = "<<thick[i]<<std::endl;
    materialArea[matnameFW[i]]+=TMath::Pi()*((thick[i]+radius)*(thick[i]+radius)-radius*radius)*nfieldWires;//*0.5;
    materialAreaFW[matnameFW[i]]+=TMath::Pi()*((thick[i]+radius)*(thick[i]+radius)-radius*radius);
    radius+=thick[i];
  }
  fWireRad = fieldWire->outerRadius();

  double totweigth=0.0, totalArea=0.0;
  std::map<std::string,double>::iterator itm=materialArea.begin();
  for(;itm!=materialArea.end();itm++){
    mat=mu2e::findMaterialOrThrow((*itm).first);
    (*itm).second*=mat->GetDensity()/(g/cm/cm/cm);
    totweigth+=(*itm).second;
    std::cout<<"material "<<mat->GetName()<<" density,g/cm^3 = "<<mat->GetDensity()/(g/cm/cm/cm)<<std::endl;
  }

  double gasarea=TMath::Pi()*(rz[1][0]*rz[1][0]-(rz[0][1]+innerGasGuardWLayer)*(rz[0][1]+innerGasGuardWLayer));
  std::string gasmat=itr->getSuperLayer(0)->getLayer(0)->getDetail()->materialName();
  mat=mu2e::findMaterialOrThrow(gasmat);
  const G4double* elemfract=mat->GetFractionVector();
  double gasPress = mat->GetPressure()/atmosphere;
  
  std::cout<<"gas "<<mat->GetName()<<" density,g/cm^3 = "<<mat->GetDensity()/(g/cm/cm/cm)<<" pressure,atm = "<<gasPress<<std::endl;
 
  double gastotweigth=0.0;
  std::vector<std::string> matnameGas;
  for(unsigned int i=0;i<mat->GetNumberOfElements();i++){
    const G4Element* elem=mat->GetElement(i);
    materialArea[elem->GetName()]=elemfract[i]*gasarea*mat->GetDensity()/(g/cm/cm/cm);//*mat->GetPressure()/(1.0*CLHEP::atmosphere);
    totweigth+=materialArea[elem->GetName()];
    gastotweigth+=materialArea[elem->GetName()];
    matnameGas.push_back( std::string(elem->GetName().c_str()) );
  }
  
  itm=materialArea.begin();
  for(;itm!=materialArea.end();itm++){
    std::cout<<"material "<<(*itm).first<<" has weigth fraction "<<(*itm).second/totweigth<<std::endl;
  }
  std::cout<<"total average gas+wire density,g/cm^3 = "<<totweigth/gasarea<<std::endl;

  double gasDensity, swDensity, fwDensity;

  std::string *ptrgasname=new std::string("ITgasAuto");
  density = gastotweigth/gasarea; // g/cm^3
  radlen = -10;      // g/cm^2
  temperature = 20.;  // deg
  pressure = gasPress;      // atm
  state="gas"; // solid or gas
  nbrcomp=matnameGas.size(); // + plus by weigth; - by atoms
  Compflag.clear();Compweight.clear();Compname.clear();
  for (std::vector<std::string>::iterator imName = matnameGas.begin(); imName != matnameGas.end(); ++imName) {
          Compflag.push_back(0);
          Compweight.push_back(materialArea[*imName]/gastotweigth);
          std::cout<<"G4 replace "<<*imName<<" by MatEnv "<<matnamemap[*imName]<<std::endl;
          Compname.push_back(matnamemap[*imName]);
  }
  Obj = new MatMaterialObj(*ptrgasname, density, zeff,
                           aeff, nbrcomp, &Compflag, &Compweight,
                           &Compname, radlen, intlen, refindex,
                           temperature, pressure, state);
  (*matobj)[ptrgasname]=Obj;
  gasDensity=density;

//  std::string *ptrgaswirename=new std::string(_GasMaterial.c_str());
//  density = totweigth/gasarea; // g/cm^3
//  radlen = -10;      // g/cm^2
//  temperature = 20.;  // deg
//  pressure = gasPress;      // atm
//  state="gas"; // solid or gas
//  nbrcomp=materialArea.size(); // + plus by weigth; - by atoms
//  Compflag.clear();Compweight.clear();Compname.clear();
//  itm=materialArea.begin();
//  for(;itm!=materialArea.end();itm++){
//          Compflag.push_back(0);
//          Compweight.push_back((*itm).second/totweigth);
//          std::cout<<"G4 replace "<<(*itm).first<<" by MatEnv "<<matnamemap[(*itm).first]<<std::endl;
//          Compname.push_back(matnamemap[(*itm).first]);
//  }
//  Obj = new MatMaterialObj(*ptrgaswirename, density, zeff,
//			   aeff, nbrcomp, &Compflag, &Compweight,
//			   &Compname, radlen, intlen, refindex,
//			   temperature, pressure, state);
//  (*matobj)[ptrgaswirename]=Obj;

  std::string *ptrfwrname=0, *ptrswrname=0;

  std::cout<<"for ITFWireAuto"<<std::endl;
  ptrfwrname=new std::string("ITFWireAuto");
  totweigth=0.0, totalArea=0.0;
  radlen = 0.0;
  for (std::vector<std::string>::iterator imName = matnameFW.begin(); imName != matnameFW.end(); ++imName) {
          G4Material* mat=mu2e::findMaterialOrThrow(*imName);
          //std::cout<<"mat->GetRadlen "<<mat->GetRadlen()<<" mat->GetDensity "<< mat->GetDensity()<<" materialAreaFW "<<materialAreaFW[*imName]<<std::endl;
          totalArea+=materialAreaFW[*imName];
          radlen += materialAreaFW[*imName]/mat->GetRadlen();
          materialAreaFW[*imName]*=mat->GetDensity()/(g/cm/cm/cm);
          totweigth+=materialAreaFW[*imName];
          std::cout<<"material "<<mat->GetName()<<" density,g/cm^3 = "<<mat->GetDensity()/(g/cm/cm/cm)<<std::endl;
  }
  density = totweigth/totalArea; // g/cm^3
  radlen = totalArea/radlen/cm*density;      // g/cm^2
  temperature = 20.;  // deg
  pressure = 1.;      // atm
  state="solid"; // solid or gas
  nbrcomp=matnameFW.size(); // + plus by weigth; - by atoms
  Compflag.clear();Compweight.clear();Compname.clear();
  for (std::vector<std::string>::iterator imName = matnameFW.begin(); imName != matnameFW.end(); ++imName) {
          Compflag.push_back(0);
          Compweight.push_back(materialAreaFW[*imName]/totweigth);
          std::cout<<"G4 replace "<<*imName<<" by MatEnv "<<matnamemap[*imName]<<std::endl;
          Compname.push_back(matnamemap[*imName]);
  }

  Obj = new MatMaterialObj(*ptrfwrname, density, zeff,
                  aeff, nbrcomp, &Compflag, &Compweight,
                  &Compname, radlen, intlen, refindex,
                  temperature, pressure, state);
  (*matobj)[ptrfwrname]=Obj;
  fwDensity=density;

  std::cout<<"for ITSWireAuto"<<std::endl;
  ptrswrname=new std::string("ITSWireAuto");
  totweigth=0.0, totalArea=0.0;
  radlen = 0.0;
  for (std::vector<std::string>::iterator imName = matnameSW.begin(); imName != matnameSW.end(); ++imName) {
          G4Material* mat=mu2e::findMaterialOrThrow(*imName);
          totalArea+=materialAreaSW[*imName];
          radlen += materialAreaSW[*imName]/mat->GetRadlen();
          materialAreaSW[*imName]*=mat->GetDensity()/(g/cm/cm/cm);
          totweigth+=materialAreaSW[*imName];
          std::cout<<"material "<<mat->GetName()<<" density,g/cm^3 = "<<mat->GetDensity()/(g/cm/cm/cm)<<std::endl;
  }
  density = totweigth/totalArea; // g/cm^3
  radlen = totalArea/radlen/cm*density;      // g/cm^2
  temperature = 20.;  // deg
  pressure = 1.;      // atm
  state="solid"; // solid or gas
  nbrcomp=matnameSW.size(); // + plus by weigth; - by atoms
  Compflag.clear();Compweight.clear();Compname.clear();
  for (std::vector<std::string>::iterator imName = matnameSW.begin(); imName != matnameSW.end(); ++imName) {
          Compflag.push_back(0);
          Compweight.push_back(materialAreaSW[*imName]/totweigth);
          std::cout<<"G4 replace "<<*imName<<" by MatEnv "<<matnamemap[*imName]<<std::endl;
          Compname.push_back(matnamemap[*imName]);
  }

  Obj = new MatMaterialObj(*ptrswrname, density, zeff,
                  aeff, nbrcomp, &Compflag, &Compweight,
                  &Compname, radlen, intlen, refindex,
                  temperature, pressure, state);
  (*matobj)[ptrswrname]=Obj;
  swDensity=density;

  std::cout<<"nFWires "<<nfieldWires<<" FWRad "<<fWireRad<<"nSWires "<<nsenseWires<<" SWRad "<<sWireRad<<std::endl;
  double fwPerc(0.0), swPerc(0.0);
  //gas_wire_extWght=true;
  if (gas_wire_extWght) {
          fwPerc =8.6e-5;
          swPerc =3.7e-6;
  } else {
          fwPerc = CLHEP::pi*nfieldWires*fWireRad*fWireRad/(gasarea-CLHEP::pi*(nfieldWires*fWireRad*fWireRad+nsenseWires*sWireRad*sWireRad));
          swPerc = CLHEP::pi*nsenseWires*sWireRad*sWireRad/(gasarea-CLHEP::pi*(nfieldWires*fWireRad*fWireRad+nsenseWires*sWireRad*sWireRad));
  }
  double gasPerc = 1.0-fwPerc-swPerc;
  std::cout<<"fw perc "<<fwPerc<<" sw perc "<<swPerc<<std::endl;

  std::string *ptrgaswirename=new std::string(_GasMaterial.c_str());
  density = gasPerc*gasDensity + fwPerc*fwDensity + swPerc*swDensity; // g/cm^3
  radlen = -10;      // g/cm^2
  temperature = 20.;  // deg
  pressure = gasPress;      // atm
  state="gas"; // solid or gas
  nbrcomp=3; // + plus by weigth; - by atoms
  Compflag.clear();Compweight.clear();Compname.clear();
  Compflag.push_back(1);
  Compweight.push_back(gasPerc*gasDensity/density);
  Compname.push_back("ITgasAuto");
  Compflag.push_back(1);
  Compweight.push_back(fwPerc*fwDensity/density);
  Compname.push_back("ITFWireAuto");
  Compflag.push_back(1);
  Compweight.push_back(swPerc*swDensity/density);
  Compname.push_back("ITSWireAuto");

  Obj = new MatMaterialObj(*ptrgaswirename, density, zeff,
                   aeff, nbrcomp, &Compflag, &Compweight,
                   &Compname, radlen, intlen, refindex,
                   temperature, pressure, state);
  (*matobj)[ptrgaswirename]=Obj;

  std::cout<<"for ITSignalCableAuto"<<std::endl;
  std::string *ptrsgncblname=new std::string("ITSignalCableAuto");
  totweigth=0.0, totalArea=0.0;
  radlen = 0.0;

  std::vector<std::string> matnameSgnCbl;
  matnameSgnCbl.push_back( config.getString("itracker.cableDielMaterial") );
  matnameSgnCbl.push_back( config.getString("itracker.cableWireMaterial") );
  std::vector<G4Material *> matnameSgnCblG4Mat;
  std::map<std::string,double> materialAreaSgnCbl;
  double cableDielWidth = config.getDouble("itracker.cableDielWidth");
  double cableWireDiameter = config.getDouble("itracker.cableWireDiameter");
  double cableWirePitch = config.getDouble("itracker.cableWirePitch");
  double nWiresInSgnCbl = std::floor(cableDielWidth/cableWirePitch);
  materialAreaSgnCbl[matnameSgnCbl.at(0)] = config.getDouble("itracker.cableDielThickness")*cableDielWidth;
  materialAreaSgnCbl[matnameSgnCbl.at(1)] = CLHEP::pi * 0.25*cableWireDiameter*cableWireDiameter * nWiresInSgnCbl;

  for (std::vector<std::string>::iterator imName = matnameSgnCbl.begin(); imName != matnameSgnCbl.end(); ++imName) {
          G4Material* mat=mu2e::findMaterialOrThrow(*imName);
          matnameSgnCblG4Mat.push_back(mat);
          //std::cout<<"mat->GetRadlen "<<mat->GetRadlen()<<" mat->GetDensity "<< mat->GetDensity()<<" materialAreaSgnCbl "<<materialAreaSgnCbl[*imName]<<std::endl;
          totalArea+=materialAreaSgnCbl[*imName];
          radlen += materialAreaSgnCbl[*imName]/mat->GetRadlen();
          materialAreaSgnCbl[*imName]*=mat->GetDensity()/(g/cm/cm/cm);
          totweigth+=materialAreaSgnCbl[*imName];
          std::cout<<"material "<<mat->GetName()<<" density,g/cm^3 = "<<mat->GetDensity()/(g/cm/cm/cm)<<std::endl;
  }
//  density = totweigth/totalArea; // g/cm^3
//  radlen = totalArea/radlen/cm*density;      // g/cm^2
//  temperature = 20.;  // deg
//  pressure = 1.;      // atm
//  state="solid"; // solid or gas
//  nbrcomp=matnameSgnCbl.size(); // + plus by weigth; - by atoms
//  Compflag.clear();Compweight.clear();Compname.clear();
//  for (std::vector<std::string>::iterator imName = matnameSgnCbl.begin(); imName != matnameSgnCbl.end(); ++imName) {
//          Compflag.push_back(1);
//          Compweight.push_back(materialAreaSgnCbl[*imName]/totweigth);
//          std::cout<<"G4 replace "<<*imName<<" by MatEnv "<<matnamemap[*imName]<<std::endl;
//          Compname.push_back(matnamemap[*imName]);
//  }

  double wrOfsgnCblPerc(0.0);
  //wrOfsgnCblPerc = nWiresInSgnCbl*cableWireDiameter/cableDielWidth;
  double dltOfsgnCblPerc = 1.0-wrOfsgnCblPerc;

  double dltDensity = matnameSgnCblG4Mat.at(0)->GetDensity()/(g/cm/cm/cm);
  double wrOfsgnCblDensity = matnameSgnCblG4Mat.at(1)->GetDensity()/(g/cm/cm/cm);

  density = dltOfsgnCblPerc*dltDensity + wrOfsgnCblPerc*wrOfsgnCblDensity; // g/cm^3
  radlen = -10;      // g/cm^2
  temperature = 20.;  // deg
  pressure = 1.;      // atm
  state="solid"; // solid or gas
  nbrcomp=1;//2; // + plus by weigth; - by atoms
  Compflag.clear();Compweight.clear();Compname.clear();
  Compflag.push_back(1);
  Compweight.push_back(dltOfsgnCblPerc*dltDensity/density);
  Compname.push_back(matnamemap[matnameSgnCbl.at(0)]);
//  Compflag.push_back(1);
//  Compweight.push_back(wrOfsgnCblPerc*wrOfsgnCblDensity/density);
//  Compname.push_back(matnamemap[matnameSgnCbl.at(1)]);

  Obj = new MatMaterialObj(*ptrsgncblname, density, zeff,
                   aeff, nbrcomp, &Compflag, &Compweight,
                   &Compname, radlen, intlen, refindex,
                   temperature, pressure, state);
  (*matobj)[ptrsgncblname]=Obj;

//  double sgnCblDensity = density;
//  double sgnCblTotalArea = totalArea;
//  double sgnCblwPerc(0.0);
//  sgnCblwPerc = ((double)_nLayers)*sgnCblTotalArea/(cableDielWidth*rEP->getThicknesses()->at(2));
//  gasPerc = 1.0-sgnCblwPerc;
//
//  std::string *ptrsgncblrepname=new std::string("ITSCblFEPAuto");
//  density = gasPerc*gasDensity + sgnCblwPerc*sgnCblDensity; // g/cm^3
//  radlen = -10;      // g/cm^2
//  temperature = 20.;  // deg
//  pressure = 1.;      // atm
//  state="solid"; // solid or gas
//  nbrcomp=2; // + plus by weigth; - by atoms
//  Compflag.clear();Compweight.clear();Compname.clear();
//  Compflag.push_back(1);
//  Compweight.push_back(gasPerc*gasDensity/density);
//  Compname.push_back("ITgasAuto");
//  Compflag.push_back(1);
//  Compweight.push_back(sgnCblwPerc*sgnCblDensity/density);
//  Compname.push_back("ITSignalCableAuto");
//
//  Obj = new MatMaterialObj(*ptrsgncblrepname, density, zeff,
//                   aeff, nbrcomp, &Compflag, &Compweight,
//                   &Compname, radlen, intlen, refindex,
//                   temperature, pressure, state);
//  (*matobj)[ptrsgncblrepname]=Obj;
//
//  sgnCblwPerc *= cableFractionOnREP;
//  gasPerc = 1.0-sgnCblwPerc;
//
//  std::string *ptrsgncblfepname=new std::string("ITSCblREPAuto");
//  density = gasPerc*gasDensity + sgnCblwPerc*sgnCblDensity; // g/cm^3
//  radlen = -10;      // g/cm^2
//  temperature = 20.;  // deg
//  pressure = 1.;      // atm
//  state="solid"; // solid or gas
//  nbrcomp=2; // + plus by weigth; - by atoms
//  Compflag.clear();Compweight.clear();Compname.clear();
//  Compflag.push_back(1);
//  Compweight.push_back(gasPerc*gasDensity/density);
//  Compname.push_back("ITgasAuto");
//  Compflag.push_back(1);
//  Compweight.push_back(sgnCblwPerc*sgnCblDensity/density);
//  Compname.push_back("ITSignalCableAuto");
//
//  Obj = new MatMaterialObj(*ptrsgncblfepname, density, zeff,
//                   aeff, nbrcomp, &Compflag, &Compweight,
//                   &Compname, radlen, intlen, refindex,
//                   temperature, pressure, state);
//  (*matobj)[ptrsgncblfepname]=Obj;

  std::cout<<"My materials"<<std::endl;
  MatDBInfo* mtdbinfo=new MatDBInfo;
  mtdbinfo->findDetMaterial(*ptrname)->printAll(std::cout);
  mtdbinfo->findDetMaterial(*ptrgasname)->printAll(std::cout);
  mtdbinfo->findDetMaterial(*ptrgaswirename)->printAll(std::cout);
  mtdbinfo->findDetMaterial(*ptrcfibname)->printAll(std::cout);
  mtdbinfo->findDetMaterial(*ptrkfoamname)->printAll(std::cout);
  mtdbinfo->findDetMaterial(*ptrploycname)->printAll(std::cout);
  //mtdbinfo->findDetMaterial(*ptrg10name)->printAll(std::cout);
  mtdbinfo->findDetMaterial(*ptrfwrname)->printAll(std::cout);
  mtdbinfo->findDetMaterial(*ptrswrname)->printAll(std::cout);
  mtdbinfo->findDetMaterial(*ptrsgncblname)->printAll(std::cout);
  //mtdbinfo->findDetMaterial(*ptrsgncblrepname)->printAll(std::cout);
  //mtdbinfo->findDetMaterial(*ptrsgncblfepname)->printAll(std::cout);

  //print in format of babar material file
  std::vector<MatMaterialObj*> vecmat;
  vecmat.push_back((*matobj)[ptrname]);
  vecmat.push_back((*matobj)[ptrgasname]);
  vecmat.push_back((*matobj)[ptrgaswirename]);
  vecmat.push_back((*matobj)[ptrcfibname]);
  vecmat.push_back((*matobj)[ptrkfoamname]);
  vecmat.push_back((*matobj)[ptrploycname]);
  //vecmat.push_back((*matobj)[ptrg10name]);
  vecmat.push_back((*matobj)[ptrfwrname]);
  vecmat.push_back((*matobj)[ptrswrname]);
  vecmat.push_back((*matobj)[ptrsgncblname]);
  //vecmat.push_back((*matobj)[ptrsgncblrepname]);
  //vecmat.push_back((*matobj)[ptrsgncblfepname]);

  MatMaterialList matlist(vecmat);
  matlist.print();
  matlist.getMaterialVector()->clear();//to protect from deleting

  // mtdbinfo->findDetMaterial("IT-Swire")->printAll(std::cout);
  // mtdbinfo->findDetMaterial("IT-Fwire")->printAll(std::cout);
  // mtdbinfo->findDetMaterial("IT-gas0")->printAll(std::cout);
  // mtdbinfo->findDetMaterial("CFiber")->printAll(std::cout);
}

void DchGDchP::makeDetailedWrSupp(std::vector<DchPhiSegmCyl> &epCylSubs,
                std::vector<double*> &ePPosSubs,
                std::vector<std::string> &ePMaterialSubs,
                double minZ, double totWidth, double rMin, double rMax, bool isFEP) {

        GeomHandle<ITracker> itracker;

        art::ServiceHandle<GeometryService> geom;
        SimpleConfig const& config  = geom->config();

//        double ncellLayers = itracker->nSuperLayers()*itracker->nRing();

        double zSign(-1);
        if (isFEP) zSign=1;

        double spdWebBaseWidth = config.getDouble("itracker.spdWebBaseWidth",0.0);
        double spdWebBaseThickness = config.getDouble("itracker.spdWebBaseThickness",0.0);

        if (firstSub) {
                epCylSubs.back() = DchPhiSegmCyl(rMin,rMin+spdWebBaseThickness,spdWebBaseWidth);
                ePPosSubs.back()[2]=minZ+zSign*spdWebBaseWidth*0.5;
                ePMaterialSubs.back() = matnamemap[config.getString("itracker.spdWebSpokeMaterial","ITGasVacuum")];
                firstSub=false;
        } else {
                epCylSubs.push_back( DchPhiSegmCyl(rMin,rMin+spdWebBaseThickness,spdWebBaseWidth) );
                ePPosSubs.push_back(new double[3]);
                ePPosSubs.back()[0]=ePPosSubs.back()[1]=0.0;
                ePPosSubs.back()[2]=minZ+zSign*spdWebBaseWidth*0.5;
                ePMaterialSubs.push_back(matnamemap[config.getString("itracker.spdWebSpokeMaterial","ITGasVacuum")]);
        }
        int isub=-1;
        std::cout<<"--- spdWebBase "<<++isub<< " --- IR "<<epCylSubs.back().getInnerRadius()<<" OR "<<epCylSubs.back().getOuterRadius()<<" Leng "<<epCylSubs.back().getLength()<<" center zpos "<<ePPosSubs.back()[2]<<" material "<<ePMaterialSubs.back()<<std::endl;

        epCylSubs.push_back( DchPhiSegmCyl(rMin+spdWebBaseThickness,rMax,spdWebBaseWidth) );
        ePPosSubs.push_back(new double[3]);
        ePPosSubs.back()[0]=ePPosSubs.back()[1]=0.0;
        ePPosSubs.back()[2]=minZ+zSign*spdWebBaseWidth*0.5;
        ePMaterialSubs.push_back("ITgasAuto");
        std::cout<<"--- spdWebGas "<<++isub<< " --- IR "<<epCylSubs.back().getInnerRadius()<<" OR "<<epCylSubs.back().getOuterRadius()<<" Leng "<<epCylSubs.back().getLength()<<" center zpos "<<ePPosSubs.back()[2]<<" material "<<ePMaterialSubs.back()<<std::endl;

        minZ+=zSign*spdWebBaseWidth;

        double maxSubWidth(0.0);

        //parameters for Field wire Boards
        double FwBoardWidth, FwBoardThick;
        FwBoardWidth = config.getDouble("itracker.fwBoardWidth");
        FwBoardThick = config.getDouble("itracker.fwBoardThickness");
        std::string FwBoardMat = config.getString("itracker.fwBoardMaterial","ITGasVacuum");
        double fwRelPos=minZ+zSign*FwBoardWidth*0.5;
        if (FwBoardWidth>maxSubWidth) maxSubWidth=FwBoardWidth;

        //parameters for Sense wire Boards
        double SwBoardWidth, SwBoardThick;
        SwBoardWidth = config.getDouble("itracker.swBoardWidth");
        SwBoardThick = config.getDouble("itracker.swBoardThickness");
        std::string SwBoardMat = config.getString("itracker.swBoardMaterial","ITGasVacuum");
        double swRelPos=minZ+zSign*SwBoardWidth*0.5;
        if (SwBoardWidth>maxSubWidth) maxSubWidth=SwBoardWidth;

        //parameters for spacers
        double startingSpacerRmin, startingSpacerRmax, finishingSpacerRmin, finishingSpacerRmax;
        startingSpacerRmax=0.0;
        double spacerRmin, spacerRmax, spacerBaseWidth, spacerBaseThick, spacerCoreThick;
        double spacerCoreSurfFrac, spokeAngle;
        int    spacerNumCores, spacerNumSpokesPerCore;
        spacerBaseWidth = config.getDouble("itracker.spacerBaseWidth");
        spacerBaseThick = config.getDouble("itracker.spacerBaseThickness");
        spacerCoreThick = config.getDouble("itracker.spacerCoreThickness");
        spacerCoreSurfFrac = config.getDouble("itracker.spacerCoreSurfFrac");
        spacerNumCores = config.getInt("itracker.spacerNumCores");
        spacerNumSpokesPerCore = config.getInt("itracker.spacerNumSpokesPerCore");
        std::string spacerMat = config.getString("itracker.spacerMaterial","ITGasVacuum");
        int spacerCoreNumSpokes = spacerNumSpokesPerCore*spacerNumCores;
        double spetSpokeRot = CLHEP::twopi/((double)spacerCoreNumSpokes);
        spokeAngle = spetSpokeRot*spacerCoreSurfFrac;
        double hollowSpacerAngle = spetSpokeRot-spokeAngle;
        //std::cout<<"------------------ spokeAngle "<<spokeAngle<<" hollowSpacerAngle "<<hollowSpacerAngle<<std::endl;

        double nCellPerLayer;
        //double componContWidth = motherDz-spacerBaseWidth;
        double componContRmin, /*compomContHeight, componContRmax,*/ compStepAngle;
        double componContRelPos = swRelPos+zSign*spacerBaseWidth*0.5;//minZ+zSign*(SwBoardWidth-(SwBoardWidth-spacerBaseWidth)*0.5);
        //parameters for  HV Capacitance
        std::string hvCapMat =  config.getString("itracker.hvCapMaterial","ITGasVacuum");
        double hvCapLength = 0.5*config.getDouble("itracker.hvCapLength");
        double hvCapWidth  = config.getDouble("itracker.hvCapWidth");
        double hvCapHeight = config.getDouble("itracker.hvCapHeight");
        double hvCapAngle;
        if (hvCapMat=="G4_Al") {hvCapMat="MatAluminum";}

        //if (hvCapLength>componContWidth) { throw cet::exception("GEOM") <<"The HV Capacitance exceeds its mother volume dim\n"; }
        //compomContHeight = hvCapHeight;

        //parameters for  Termination Resistance
        std::string termResMat =  config.getString("itracker.trmResMaterial","ITGasVacuum");
        double termResLength = 0.5*config.getDouble("itracker.trmResLength");
        double termResWidth  = config.getDouble("itracker.trmResWidth");
        double termResHeight = config.getDouble("itracker.trmResHeight");
        double termResAngle;
        if (termResMat=="G4_Al") {termResMat="MatAluminum";}

        //if (termResLength>componContWidth) { throw cet::exception("GEOM") <<"The termination Resistance exceeds its mother volume dim\n"; }
        //if (termResHeight>compomContHeight) { compomContHeight = termResHeight; }

        //parameters for  Hv Resistance
        /*std::string hvResMat =  config.getString("itracker.hvResMaterial","ITGasVacuum");
        double hvResLength = 0.5*config.getDouble("itracker.hvResLength");
        double hvResWidth  = config.getDouble("itracker.hvResWidth");
        double hvResHeight = config.getDouble("itracker.hvResHeight");
        double hvResAngle;
        if (hvResMat=="G4_Al") {hvResMat="MatAluminum";}*/


        double spacerBaseRelPos=minZ+zSign*spacerBaseWidth*0.5;
        if (spacerBaseWidth>maxSubWidth) maxSubWidth=spacerBaseWidth;

        if (maxSubWidth>totWidth) { return; }
        //double gasBtwSpacerRelPos=minZ+zSign*maxSubWidth*0.5;

        double minR, maxR(0.0), swR, tmpZincr;
        //int superlayer, iring, icellLayer=0;
        bool firstLayer=true;
        for (int iSl = 0; iSl < itracker->nSuperLayers(); iSl++){

                SuperLayer *SLayer = itracker->getSuperLayer(iSl);

                for (int iLy=0; iLy < SLayer->nLayers(); iLy++ ){

                        boost::shared_ptr<ITLayer> ily = SLayer->getLayer(iLy);
                        //superlayer = ily->Id().getSuperLayer();
                        //iring = ily->Id().getLayer();
                        if (ily->nCells()>0) {
                                minR = ily->getDetail()->centerInnerRadiusRing();
                                tmpZincr = ily->getDetail()->halfLength()*std::tan(ily->getDetail()->stereoAngleInnerRing());
                                minR = std::sqrt( minR*minR + tmpZincr*tmpZincr);
                                if(firstLayer) {
                                        startingSpacerRmax = minR;
                                        firstLayer = false;
                                }

                                maxR = ily->getDetail()->centerOuterRadiusRing();
                                tmpZincr = ily->getDetail()->halfLength()*std::tan(ily->getDetail()->stereoAngleOuterRing());
                                maxR = std::sqrt( maxR*maxR + tmpZincr*tmpZincr);

                                swR = (maxR+minR)*0.5;

                                //------------- start wires boards -------------
                                epCylSubs.push_back( DchPhiSegmCyl(minR,minR+FwBoardThick,FwBoardWidth) );
                                ePPosSubs.push_back(new double[3]);
                                ePPosSubs.back()[0]=ePPosSubs.back()[1]=0.0;
                                ePPosSubs.back()[2]=fwRelPos;
                                ePMaterialSubs.push_back(matnamemap[FwBoardMat]);
                                std::cout<<"--- FwBoard Dwn layer "<<iLy<<" sub "<<++isub<< " --- IR "<<epCylSubs.back().getInnerRadius()<<" OR "<<epCylSubs.back().getOuterRadius()<<" Leng "<<epCylSubs.back().getLength()<<" center zpos "<<ePPosSubs.back()[2]<<" material "<<ePMaterialSubs.back()<<std::endl;

                                epCylSubs.push_back( DchPhiSegmCyl(maxR-FwBoardThick,maxR,FwBoardWidth) );
                                ePPosSubs.push_back(new double[3]);
                                ePPosSubs.back()[0]=ePPosSubs.back()[1]=0.0;
                                ePPosSubs.back()[2]=fwRelPos;
                                ePMaterialSubs.push_back(matnamemap[FwBoardMat]);
                                std::cout<<"--- FwBoard Up layer "<<iLy<<" sub "<<++isub<< " --- IR "<<epCylSubs.back().getInnerRadius()<<" OR "<<epCylSubs.back().getOuterRadius()<<" Leng "<<epCylSubs.back().getLength()<<" center zpos "<<ePPosSubs.back()[2]<<" material "<<ePMaterialSubs.back()<<std::endl;

                                epCylSubs.push_back( DchPhiSegmCyl(swR-SwBoardThick,swR,SwBoardWidth) );
                                ePPosSubs.push_back(new double[3]);
                                ePPosSubs.back()[0]=ePPosSubs.back()[1]=0.0;
                                ePPosSubs.back()[2]=swRelPos;
                                ePMaterialSubs.push_back(matnamemap[SwBoardMat]);
                                std::cout<<"--- SwBoard layer "<<iLy<<" sub "<<++isub<< " --- IR "<<epCylSubs.back().getInnerRadius()<<" OR "<<epCylSubs.back().getOuterRadius()<<" Leng "<<epCylSubs.back().getLength()<<" center zpos "<<ePPosSubs.back()[2]<<" material "<<ePMaterialSubs.back()<<std::endl;
                                //------------- end wires boards -------------

                                //------------- start spacers -------------
                                //start down spacers
                                spacerRmin = minR+FwBoardThick;
                                spacerRmax = swR-SwBoardThick;

                                epCylSubs.push_back( DchPhiSegmCyl(spacerRmin,spacerRmin+spacerBaseThick,spacerBaseWidth) );
                                ePPosSubs.push_back(new double[3]);
                                ePPosSubs.back()[0]=ePPosSubs.back()[1]=0.0;
                                ePPosSubs.back()[2]=spacerBaseRelPos;
                                ePMaterialSubs.push_back(matnamemap[spacerMat]);
                                std::cout<<"--- Dwn spacer Dwn layer "<<iLy<<" sub "<<++isub<< " --- IR "<<epCylSubs.back().getInnerRadius()<<" OR "<<epCylSubs.back().getOuterRadius()<<" Leng "<<epCylSubs.back().getLength()<<" center zpos "<<ePPosSubs.back()[2]<<" material "<<ePMaterialSubs.back()<<std::endl;

                                ////fill space between spacer bases gas layer
                                //epCylSubs.push_back( DchPhiSegmCyl(spacerRmin+spacerBaseThick,spacerRmax-spacerBaseThick,maxSubWidth) );
                                //ePPosSubs.push_back(new double[3]);
                                //ePPosSubs.back()[0]=ePPosSubs.back()[1]=0.0;
                                //ePPosSubs.back()[2]=gasBtwSpacerRelPos;
                                //ePMaterialSubs.push_back("ITgasAuto");
                                ////end fill
                                //spacers spokes
                                epCylSubs.push_back( DchPhiSegmCyl(spacerRmin+spacerBaseThick,spacerRmax-spacerBaseThick,spacerCoreThick,-0.5*spokeAngle,spokeAngle,hollowSpacerAngle) );
                                ePPosSubs.push_back(new double[3]);
                                ePPosSubs.back()[0]=ePPosSubs.back()[1]=0.0;
                                ePPosSubs.back()[2]=spacerBaseRelPos;
                                ePMaterialSubs.push_back(matnamemap[spacerMat]);
                                std::cout<<"--- Dwn spacer Spoke layer "<<iLy<<" sub "<<++isub<< " --- IR "<<epCylSubs.back().getInnerRadius()<<" OR "<<epCylSubs.back().getOuterRadius()<<" Leng "<<epCylSubs.back().getLength()<<" center zpos "<<ePPosSubs.back()[2]<<" material "<<ePMaterialSubs.back()<<std::endl;
                                //end spokes


                                epCylSubs.push_back( DchPhiSegmCyl(spacerRmax-spacerBaseThick,spacerRmax,spacerBaseWidth) );
                                ePPosSubs.push_back(new double[3]);
                                ePPosSubs.back()[0]=ePPosSubs.back()[1]=0.0;
                                ePPosSubs.back()[2]=spacerBaseRelPos;
                                ePMaterialSubs.push_back(matnamemap[spacerMat]);
                                std::cout<<"--- Dwn spacer Up layer "<<iLy<<" sub "<<++isub<< " --- IR "<<epCylSubs.back().getInnerRadius()<<" OR "<<epCylSubs.back().getOuterRadius()<<" Leng "<<epCylSubs.back().getLength()<<" center zpos "<<ePPosSubs.back()[2]<<" material "<<ePMaterialSubs.back()<<std::endl;
                                //end down spacers

                                //start up spacers
                                spacerRmin = swR;
                                spacerRmax = maxR-FwBoardThick;

                                epCylSubs.push_back( DchPhiSegmCyl(spacerRmin,spacerRmin+spacerBaseThick,spacerBaseWidth) );
                                ePPosSubs.push_back(new double[3]);
                                ePPosSubs.back()[0]=ePPosSubs.back()[1]=0.0;
                                ePPosSubs.back()[2]=spacerBaseRelPos;
                                ePMaterialSubs.push_back(matnamemap[spacerMat]);
                                std::cout<<"--- Up spacer Dwn layer "<<iLy<<" sub "<<++isub<< " --- IR "<<epCylSubs.back().getInnerRadius()<<" OR "<<epCylSubs.back().getOuterRadius()<<" Leng "<<epCylSubs.back().getLength()<<" center zpos "<<ePPosSubs.back()[2]<<" material "<<ePMaterialSubs.back()<<std::endl;

                                ////fill space between spacer bases gas layer
                                //epCylSubs.push_back( DchPhiSegmCyl(spacerRmin+spacerBaseThick,spacerRmax-spacerBaseThick,maxSubWidth) );
                                //ePPosSubs.push_back(new double[3]);
                                //ePPosSubs.back()[0]=ePPosSubs.back()[1]=0.0;
                                //ePPosSubs.back()[2]=gasBtwSpacerRelPos;
                                //ePMaterialSubs.push_back("ITgasAuto");
                                ////end fill
                                //spacers spokes
                                epCylSubs.push_back( DchPhiSegmCyl(spacerRmin+spacerBaseThick,spacerRmax-spacerBaseThick,spacerCoreThick,-0.5*spokeAngle,spokeAngle,hollowSpacerAngle) );
                                ePPosSubs.push_back(new double[3]);
                                ePPosSubs.back()[0]=ePPosSubs.back()[1]=0.0;
                                ePPosSubs.back()[2]=spacerBaseRelPos;
                                ePMaterialSubs.push_back(matnamemap[spacerMat]);
                                std::cout<<"--- Up spacer Spoke layer "<<iLy<<" sub "<<++isub<< " --- IR "<<epCylSubs.back().getInnerRadius()<<" OR "<<epCylSubs.back().getOuterRadius()<<" Leng "<<epCylSubs.back().getLength()<<" center zpos "<<ePPosSubs.back()[2]<<" material "<<ePMaterialSubs.back()<<std::endl;
                                //end spokes

                                epCylSubs.push_back( DchPhiSegmCyl(spacerRmax-spacerBaseThick,spacerRmax,spacerBaseWidth) );
                                ePPosSubs.push_back(new double[3]);
                                ePPosSubs.back()[0]=ePPosSubs.back()[1]=0.0;
                                ePPosSubs.back()[2]=spacerBaseRelPos;
                                ePMaterialSubs.push_back(matnamemap[spacerMat]);
                                std::cout<<"--- Dwn spacer Up layer "<<iLy<<" sub "<<++isub<< " --- IR "<<epCylSubs.back().getInnerRadius()<<" OR "<<epCylSubs.back().getOuterRadius()<<" Leng "<<epCylSubs.back().getLength()<<" center zpos "<<ePPosSubs.back()[2]<<" material "<<ePMaterialSubs.back()<<std::endl;
                               //end up spacers
                                //------------- end spacers -------------

                                nCellPerLayer = ily->nCells();
                                compStepAngle = CLHEP::twopi/nCellPerLayer;

                                //------------- start Term components -------------
                                componContRmin = swR;
                                //componContRmax = componContRmin+compomContHeight;

                                //start HV capacitance
                                hvCapAngle = hvCapWidth/componContRmin;
                                epCylSubs.push_back( DchPhiSegmCyl(componContRmin,componContRmin+hvCapHeight,hvCapLength,-0.5*hvCapAngle,hvCapAngle,compStepAngle-hvCapAngle) );
                                ePPosSubs.push_back(new double[3]);
                                ePPosSubs.back()[0]=ePPosSubs.back()[1]=0.0;
                                ePPosSubs.back()[2]=componContRelPos;
                                ePMaterialSubs.push_back(matnamemap[hvCapMat]);
                                std::cout<<"--- HV cap layer "<<iLy<<" sub "<<++isub<< " --- IR "<<epCylSubs.back().getInnerRadius()<<" OR "<<epCylSubs.back().getOuterRadius()<<" Leng "<<epCylSubs.back().getLength()<<" center zpos "<<ePPosSubs.back()[2]<<" material "<<ePMaterialSubs.back()<<std::endl;

                                //end HV capacitance

                                //start Term Resistance
                                termResAngle = termResWidth/componContRmin;
                                epCylSubs.push_back( DchPhiSegmCyl(componContRmin,componContRmin+termResHeight,termResLength,0.5*(compStepAngle-termResAngle),termResAngle,compStepAngle-termResAngle) );
                                ePPosSubs.push_back(new double[3]);
                                ePPosSubs.back()[0]=ePPosSubs.back()[1]=0.0;
                                ePPosSubs.back()[2]=componContRelPos;
                                ePMaterialSubs.push_back(matnamemap[termResMat]);
                                std::cout<<"--- Term res layer "<<iLy<<" sub "<<++isub<< " --- IR "<<epCylSubs.back().getInnerRadius()<<" OR "<<epCylSubs.back().getOuterRadius()<<" Leng "<<epCylSubs.back().getLength()<<" center zpos "<<ePPosSubs.back()[2]<<" material "<<ePMaterialSubs.back()<<std::endl;
                                //end Term Resistance
                                //------------- end Term components -------------

                                //++icellLayer;
                        }

                }

        }

        minR = rMin+spdWebBaseThickness;//((G4Tubs *)localMother->GetSolid())->GetInnerRadius();
        startingSpacerRmin = minR;

        //------------- start starting spacer for guard field wires -------------
        epCylSubs.push_back( DchPhiSegmCyl(startingSpacerRmax-FwBoardThick,startingSpacerRmax,FwBoardWidth) );
        ePPosSubs.push_back(new double[3]);
        ePPosSubs.back()[0]=ePPosSubs.back()[1]=0.0;
        ePPosSubs.back()[2]=fwRelPos;
        ePMaterialSubs.push_back(matnamemap[FwBoardMat]);
        std::cout<<"--- FwBoard Up layer -1 sub "<<++isub<< " --- IR "<<epCylSubs.back().getInnerRadius()<<" OR "<<epCylSubs.back().getOuterRadius()<<" Leng "<<epCylSubs.back().getLength()<<" center zpos "<<ePPosSubs.back()[2]<<" material "<<ePMaterialSubs.back()<<std::endl;

        //start down spacers
        spacerRmin = startingSpacerRmin;
        spacerRmax = startingSpacerRmax-FwBoardThick;

        epCylSubs.push_back( DchPhiSegmCyl(spacerRmin,spacerRmin+spacerBaseThick,spacerBaseWidth) );
        ePPosSubs.push_back(new double[3]);
        ePPosSubs.back()[0]=ePPosSubs.back()[1]=0.0;
        ePPosSubs.back()[2]=spacerBaseRelPos;
        ePMaterialSubs.push_back(matnamemap[spacerMat]);
        std::cout<<"--- spacer Dwn layer -1 sub "<<++isub<< " --- IR "<<epCylSubs.back().getInnerRadius()<<" OR "<<epCylSubs.back().getOuterRadius()<<" Leng "<<epCylSubs.back().getLength()<<" center zpos "<<ePPosSubs.back()[2]<<" material "<<ePMaterialSubs.back()<<std::endl;

        ////fill space between spacer bases gas layer
        //epCylSubs.push_back( DchPhiSegmCyl(spacerRmin+spacerBaseThick,spacerRmax-spacerBaseThick,maxSubWidth) );
        //ePPosSubs.push_back(new double[3]);
        //ePPosSubs.back()[0]=ePPosSubs.back()[1]=0.0;
        //ePPosSubs.back()[2]=gasBtwSpacerRelPos;
        //ePMaterialSubs.push_back("ITgasAuto");
        ////end fill
        //spacers spokes
        epCylSubs.push_back( DchPhiSegmCyl(spacerRmin+spacerBaseThick,spacerRmax-spacerBaseThick,spacerCoreThick,-0.5*spokeAngle,spokeAngle,hollowSpacerAngle) );
        ePPosSubs.push_back(new double[3]);
        ePPosSubs.back()[0]=ePPosSubs.back()[1]=0.0;
        ePPosSubs.back()[2]=spacerBaseRelPos;
        ePMaterialSubs.push_back(matnamemap[spacerMat]);
        std::cout<<"--- spacer Spoke layer -1 sub "<<++isub<< " --- IR "<<epCylSubs.back().getInnerRadius()<<" OR "<<epCylSubs.back().getOuterRadius()<<" Leng "<<epCylSubs.back().getLength()<<" center zpos "<<ePPosSubs.back()[2]<<" material "<<ePMaterialSubs.back()<<std::endl;
        //end spokes

        epCylSubs.push_back( DchPhiSegmCyl(spacerRmax-spacerBaseThick,spacerRmax,spacerBaseWidth) );
        ePPosSubs.push_back(new double[3]);
        ePPosSubs.back()[0]=ePPosSubs.back()[1]=0.0;
        ePPosSubs.back()[2]=spacerBaseRelPos;
        ePMaterialSubs.push_back(matnamemap[spacerMat]);
        std::cout<<"--- spacer Up layer -1 sub "<<++isub<< " --- IR "<<epCylSubs.back().getInnerRadius()<<" OR "<<epCylSubs.back().getOuterRadius()<<" Leng "<<epCylSubs.back().getLength()<<" center zpos "<<ePPosSubs.back()[2]<<" material "<<ePMaterialSubs.back()<<std::endl;
        //------------- end starting spacer for guard field wires -------------

        //------------- start finishing spacer for guard field wires -------------
        finishingSpacerRmin = maxR;
        finishingSpacerRmax = rMax;//((G4Tubs *)localMother->GetSolid())->GetOuterRadius();

        epCylSubs.push_back( DchPhiSegmCyl(finishingSpacerRmin,finishingSpacerRmin+FwBoardThick,FwBoardWidth) );
        ePPosSubs.push_back(new double[3]);
        ePPosSubs.back()[0]=ePPosSubs.back()[1]=0.0;
        ePPosSubs.back()[2]=fwRelPos;
        ePMaterialSubs.push_back(matnamemap[FwBoardMat]);
        std::cout<<"--- FwBoard Dwn layer max+1 sub "<<++isub<< " --- IR "<<epCylSubs.back().getInnerRadius()<<" OR "<<epCylSubs.back().getOuterRadius()<<" Leng "<<epCylSubs.back().getLength()<<" center zpos "<<ePPosSubs.back()[2]<<" material "<<ePMaterialSubs.back()<<std::endl;

        //start up spacer
        spacerRmin = finishingSpacerRmin+FwBoardThick;
        spacerRmax = finishingSpacerRmax;

        epCylSubs.push_back( DchPhiSegmCyl(spacerRmin,spacerRmin+spacerBaseThick,spacerBaseWidth) );
        ePPosSubs.push_back(new double[3]);
        ePPosSubs.back()[0]=ePPosSubs.back()[1]=0.0;
        ePPosSubs.back()[2]=spacerBaseRelPos;
        ePMaterialSubs.push_back(matnamemap[spacerMat]);
        std::cout<<"--- spacer Dwn layer max+1 sub "<<++isub<< " --- IR "<<epCylSubs.back().getInnerRadius()<<" OR "<<epCylSubs.back().getOuterRadius()<<" Leng "<<epCylSubs.back().getLength()<<" center zpos "<<ePPosSubs.back()[2]<<" material "<<ePMaterialSubs.back()<<std::endl;

        ////fill space between spacer bases gas layer
        //epCylSubs.push_back( DchPhiSegmCyl(spacerRmin+spacerBaseThick,spacerRmax-spacerBaseThick,maxSubWidth) );
        //ePPosSubs.push_back(new double[3]);
        //ePPosSubs.back()[0]=ePPosSubs.back()[1]=0.0;
        //ePPosSubs.back()[2]=gasBtwSpacerRelPos;
        //ePMaterialSubs.push_back("ITgasAuto");
        ////end fill
        //spacers spokes
        epCylSubs.push_back( DchPhiSegmCyl(spacerRmin+spacerBaseThick,spacerRmax-spacerBaseThick,spacerCoreThick,-0.5*spokeAngle,spokeAngle,hollowSpacerAngle) );
        ePPosSubs.push_back(new double[3]);
        ePPosSubs.back()[0]=ePPosSubs.back()[1]=0.0;
        ePPosSubs.back()[2]=spacerBaseRelPos;
        ePMaterialSubs.push_back(matnamemap[spacerMat]);
        std::cout<<"--- spacer Spoke layer max+1 sub "<<++isub<< " --- IR "<<epCylSubs.back().getInnerRadius()<<" OR "<<epCylSubs.back().getOuterRadius()<<" Leng "<<epCylSubs.back().getLength()<<" center zpos "<<ePPosSubs.back()[2]<<" material "<<ePMaterialSubs.back()<<std::endl;
        //end spokes

        epCylSubs.push_back( DchPhiSegmCyl(spacerRmax-spacerBaseThick,spacerRmax,spacerBaseWidth) );
        ePPosSubs.push_back(new double[3]);
        ePPosSubs.back()[0]=ePPosSubs.back()[1]=0.0;
        ePPosSubs.back()[2]=spacerBaseRelPos;
        ePMaterialSubs.push_back(matnamemap[spacerMat]);
        std::cout<<"--- spacer Up layer max+1 sub "<<++isub<< " --- IR "<<epCylSubs.back().getInnerRadius()<<" OR "<<epCylSubs.back().getOuterRadius()<<" Leng "<<epCylSubs.back().getLength()<<" center zpos "<<ePPosSubs.back()[2]<<" material "<<ePMaterialSubs.back()<<std::endl;

        //end up spacer

        //remaining gas layer
        epCylSubs.push_back( DchPhiSegmCyl(rMin,rMax,totWidth-maxSubWidth) );
        ePPosSubs.push_back(new double[3]);
        ePPosSubs.back()[0]=ePPosSubs.back()[1]=0.0;
        ePPosSubs.back()[2]=(minZ-zSign*spdWebBaseWidth)+zSign*(totWidth+maxSubWidth)*0.5;//minZ+zSign*totWidth-zSign*(totWidth-maxSubWidth)*0.5
        ePMaterialSubs.push_back("ITgasAuto");
        std::cout<<"--- remaining gas layer "<<++isub<< " --- IR "<<epCylSubs.back().getInnerRadius()<<" OR "<<epCylSubs.back().getOuterRadius()<<" Leng "<<epCylSubs.back().getLength()<<" center zpos "<<ePPosSubs.back()[2]<<" material "<<ePMaterialSubs.back()<<std::endl;

}

void DchGDchP::makeInEPCables(std::vector<DchPhiSegmCyl> &epCylSubs,
                std::vector<double*> &ePPosSubs,
                std::vector<std::string> &ePMaterialSubs,
                double minZ, double totWidth, double rMin, double rMax, bool isREP) {

        GeomHandle<ITracker> itracker;

        art::ServiceHandle<GeometryService> geom;
        SimpleConfig const& config  = geom->config();

        double zRelPos = minZ;
        if(isREP) { zRelPos-=totWidth*0.5; }
        else { zRelPos+=totWidth*0.5; }

        double ncellLayers = itracker->nSuperLayers()*itracker->nRing();
        double minR, maxR, swR(0.0), tmpZincr;
        bool foundFirstLayer=false;
        for (int iSl = 0; iSl < itracker->nSuperLayers(); iSl++){

                SuperLayer *SLayer = itracker->getSuperLayer(iSl);

                for (int iLy=0; iLy < SLayer->nLayers(); iLy++ ){

                        boost::shared_ptr<ITLayer> ily = SLayer->getLayer(iLy);
                        if (ily->nCells()>0) {
                                minR = ily->getDetail()->centerInnerRadiusRing();
                                tmpZincr = ily->getDetail()->halfLength()*std::tan(ily->getDetail()->stereoAngleInnerRing());
                                minR = std::sqrt( minR*minR + tmpZincr*tmpZincr);

                                maxR = ily->getDetail()->centerOuterRadiusRing();
                                tmpZincr = ily->getDetail()->halfLength()*std::tan(ily->getDetail()->stereoAngleOuterRing());
                                maxR = std::sqrt( maxR*maxR + tmpZincr*tmpZincr);

                                swR = (maxR+minR)*0.5;

                                foundFirstLayer=true;
                                break;
                        }
                }
                if (foundFirstLayer) { break;}
        }

        //start signal cables
        double cableDielWidth = config.getDouble("itracker.cableDielWidth");
        double cableWireDiameter = config.getDouble("itracker.cableWireDiameter");
        double cableDielThickness = config.getDouble("itracker.cableDielThickness");
        //double cableWirePitch = config.getDouble("itracker.cableWirePitch");
        //double nWiresInSgnCbl = std::floor(cableDielWidth/cableWirePitch);

        double sgnCblAngle = cableDielWidth/swR;
        sgnCblAngle += cableDielWidth/rMax;
        sgnCblAngle *= 0.5;
        double spetSpokeRot = CLHEP::twopi/((double)config.getInt("itracker.spdWebSpokesNumber"));
        double hollowSgnCblAngle = spetSpokeRot-sgnCblAngle;
//        epCylSubs.push_back( DchPhiSegmCyl(swR,rMax,totWidth,-0.5*sgnCblAngle,sgnCblAngle,hollowSgnCblAngle) );
//        ePPosSubs.push_back(new double[3]);
//        ePPosSubs.back()[0]=ePPosSubs.back()[1]=0.0;
//        ePPosSubs.back()[2]=zRelPos;
//        if (isREP) { ePMaterialSubs.push_back("ITSCblREPAuto"); }
//        else { ePMaterialSubs.push_back("ITSCblFEPAuto"); }
//        std::cout<<"--- Signal Cables "<< " --- IR "<<epCylSubs.back().getInnerRadius()<<" OR "<<epCylSubs.back().getOuterRadius()<<" Leng "<<epCylSubs.back().getLength()<<" center zpos "<<ePPosSubs.back()[2]<<" material "<<ePMaterialSubs.back()<<std::endl;
//        std::cout<<"--- is sub n "<<epCylSubs.size()<<std::endl;
//        //end signal cables
//
//        //start gas between adiacent wire packs
//        epCylSubs.push_back( DchPhiSegmCyl(swR,rMax,totWidth,-0.5*sgnCblAngle-hollowSgnCblAngle,hollowSgnCblAngle,sgnCblAngle) );
//        //if (isREP) { epCylSubs.push_back( DchPhiSegmCyl(swR,rMax,totWidth,-0.5*sgnCblAngle-hollowSgnCblAngle,hollowSgnCblAngle,sgnCblAngle) ); }
//        //else { epCylSubs.push_back( DchPhiSegmCyl(swR,rMax,totWidth,0.5*sgnCblAngle,hollowSgnCblAngle,sgnCblAngle) ); }
//        ePPosSubs.push_back(new double[3]);
//        ePPosSubs.back()[0]=ePPosSubs.back()[1]=0.0;
//        ePPosSubs.back()[2]=zRelPos;
//        ePMaterialSubs.push_back("ITgasAuto");
//        //end gas between adiacent wire packs

        double cblTotThik = cableDielThickness+cableWireDiameter;
        double eqlCblThik = cableDielThickness;
        double cblSpet = totWidth/ncellLayers;
        double gasBtwCblThik = cblSpet-cblTotThik;
        double iCblZpos(minZ), dzSign(1.0);
        if(isREP) { dzSign=-1.0; }
        iCblZpos+=dzSign*(cblTotThik-eqlCblThik*0.5);
        for (int iCbl=0; iCbl<ncellLayers; ++iCbl) {
                epCylSubs.push_back( DchPhiSegmCyl(swR,rMax,eqlCblThik,-0.5*sgnCblAngle,sgnCblAngle,hollowSgnCblAngle) );
                ePPosSubs.push_back(new double[3]);
                ePPosSubs.back()[0]=ePPosSubs.back()[1]=0.0;
                ePPosSubs.back()[2]=iCblZpos;
                ePMaterialSubs.push_back("ITSignalCableAuto");
                //if (isREP) { ePMaterialSubs.push_back("ITSCblREPAuto"); }
                //else { ePMaterialSubs.push_back("ITSCblFEPAuto"); }
                std::cout<<"--- Signal Cables "<< " --- IR "<<epCylSubs.back().getInnerRadius()<<" OR "<<epCylSubs.back().getOuterRadius()<<" Leng "<<epCylSubs.back().getLength()<<" center zpos "<<ePPosSubs.back()[2]<<" material "<<ePMaterialSubs.back()<<std::endl;
                std::cout<<"--- is sub n "<<epCylSubs.size()<<std::endl;

                //if (iCbl<(ncellLayers-1)) {
                        //start gas between adiacent wire packs
                        epCylSubs.push_back( DchPhiSegmCyl(swR,rMax,gasBtwCblThik) );
                        ePPosSubs.push_back(new double[3]);
                        ePPosSubs.back()[0]=ePPosSubs.back()[1]=0.0;
                        ePPosSubs.back()[2]=iCblZpos+dzSign*gasBtwCblThik*0.5;
                        ePMaterialSubs.push_back("ITgasAuto");
                        std::cout<<"--- gas between Signal Cables "<< " --- IR "<<epCylSubs.back().getInnerRadius()<<" OR "<<epCylSubs.back().getOuterRadius()<<" Leng "<<epCylSubs.back().getLength()<<" center zpos "<<ePPosSubs.back()[2]<<" material "<<ePMaterialSubs.back()<<std::endl;
                        std::cout<<"--- is sub n "<<epCylSubs.size()<<std::endl;
                        //end gas between adiacent wire packs
                //}
                iCblZpos+=dzSign*cblSpet;
        }
        //end signal cables

        //start remaining gas
        epCylSubs.push_back( DchPhiSegmCyl(rMin,swR,totWidth) );
        ePPosSubs.push_back(new double[3]);
        ePPosSubs.back()[0]=ePPosSubs.back()[1]=0.0;
        ePPosSubs.back()[2]=zRelPos;
        ePMaterialSubs.push_back("ITgasAuto");
        //end remaining gas

}
