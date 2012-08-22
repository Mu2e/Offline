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

using namespace mu2e;

DchGDchP::DchGDchP(std::string materialfile):
  DchGDch(nLayer)
{
  art::ServiceHandle<GeometryService> geom;
  GeomHandle<ITracker> itr;
  CellGeometryHandle *itwp = itr->getCellGeometryHandle();

  bool automaterial=true;
  if(materialfile!="") automaterial=false;
  std::cout<<"materialdb "<<materialfile<<" auto re-build "<<automaterial<<std::endl;
  //std::cout<<"--------010--------"<<std::endl<<std::flush;
  //WorldMaker* allMu2e    = new WorldMaker();
  ///*G4Material* mat=*/mu2e::findMaterialOrThrow("ITGasHe_90Isob_10");//only to fix a bug in running the Kalman as analyzer FIXME
  //std::cout<<"--------020--------"<<std::endl<<std::flush;

  double rz[2][3];
  double inraddensity=0;
  
  std::multimap<Wall::Walltype,boost::shared_ptr<Wall> >::iterator it;
  for(it=itr->getWalls()->begin();it!=itr->getWalls()->end();it++){
    if((*it).first==Wall::inner||(*it).first==Wall::outer){
      int ind=((*it).first==Wall::inner)?0:1;
      std::cout<<"barrel radius "<<(*it).second->getRmin()<<" "<<(*it).second->getRmax()
	       <<" "<<(*it).second->getDz()<<std::endl;
      rz[ind][0]=(*it).second->getRmin();
      rz[ind][1]=(*it).second->getRmax();
      rz[ind][2]=(*it).second->getDz();

      if((*it).first==Wall::inner){
	const std::vector<double>* th = (*it).second->getThicknesses().get();
	const std::vector<std::string>* mats = (*it).second->getMaterialsName().get();
	double thickness=0;
	for(unsigned int i=0;i<th->size();i++){
	  G4String matname=(*mats)[i];
	  if (matname.contains("ITGas")) {
	          rz[ind][1]-=(*th)[i];
	  } else if (automaterial) {
	          G4Material* mat=mu2e::findMaterialOrThrow(matname);
	          std::cout<<"sum material "<<mat->GetName()<<" thick,mm= "<<(*th)[i]
                       <<" density,g/cm^3 = "<<mat->GetDensity()/(g/cm/cm/cm)<<std::endl;
	          inraddensity+=(*th)[i]*mat->GetDensity()/(g/cm/cm/cm);
	          thickness+=(*th)[i];
      }
	}
	inraddensity/=thickness;
	std::cout<<"average density = "<<inraddensity
		 <<" g/cm^3, at thickness = "<<thickness<<" mm"<<std::endl;
      }
    }
  }


  _name="Mu2e ITracker";    
  initWires();
  _version=999;
  _nLayers=itr->nSuperLayers();
  _GasCyl=DchSimpleCyl(rz[0][1],rz[1][0],itr->zHalfLength()*2);
  _GasPos[0]=0;_GasPos[1]=0;_GasPos[2]=0;
  _InnerCyl=DchSimpleCyl(rz[0][0],rz[0][1],itr->zHalfLength()*2);
  _ICPos[0]=0;_ICPos[1]=0;_ICPos[2]=0;
  _OuterCyl=DchSimpleCyl(rz[1][0],rz[1][1],itr->zHalfLength()*2);
  _OCPos[0]=0;_OCPos[1]=0;_OCPos[2]=0;
  _RearCyl=DchSimpleCyl(itr->r0()-0.1,itr->rOut()+0.1,0.1);
  _REPPos[0]=0;_REPPos[1]=0;_REPPos[2]=itr->zHalfLength()+0.1;
  _ForwCyl=DchSimpleCyl(itr->r0()-0.1,itr->rOut()+0.1,0.1);
  _FEPPos[0]=0;_FEPPos[1]=0;_FEPPos[2]=itr->zHalfLength()+0.1;
  for(int i=0;i<_nLayers;i++){
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
         _sWires[i].diameter=0.005;
         _sWires[i].volts=2000;
         std::cout<<i<<" "<<itwp->GetWireCenter()<<" n="<<itwp->GetWireDirection()<<std::endl;

         _cellStruc[i].nofFWires=0;
         _rOffset[i]=0;
         break;
       }
    }
  }
  _zoffset=0;
  _zlen=0;
  _cellHeight=0;

  RecoMatFactory* matFactory = RecoMatFactory::getInstance();

  if(!automaterial){
    mu2e::ConfigFileLookupPolicy findFile;
    std::string fullPath = findFile(materialfile);
    MatMaterialList* mtrList = new MatMaterialList(fullPath);
    RecoMatFactory* matFactory = RecoMatFactory::getInstance();
    ((MatMtrDictionary*)matFactory->materialDictionary())->FillMtrDict(mtrList);
    std::cout<<"My materials"<<std::endl;
    MatDBInfo* mtdbinfo=new MatDBInfo;
    mtdbinfo->findDetMaterial("ITInBarrelAuto")->printAll(std::cout);
    mtdbinfo->findDetMaterial("ITgasAuto")->printAll(std::cout);
    return;
  }

  //create average materials from geometry description 
  std::map<std::string*, MatMaterialObj*, PtrLess>* matobj=matFactory->materialDictionary();

  double zeff = 0., aeff = 0., intlen = 0., refindex = 0.; //auto produced

  std::string *ptrname=new std::string("ITInBarrelAuto");
  double density = inraddensity; // g/cm^3
  double radlen = 42.7;      // g/cm^2
  double temperature = 20.;  // deg
  double pressure = 1.;      // atm
  std::string state="solid"; // solid or gas
  int nbrcomp=1; // + plus by weigth; - by atoms 
  std::vector<int> Compflag;Compflag.push_back(0); // 0 - element, 1 - material
  std::vector<double> Compweight;Compweight.push_back(1.); 
  std::vector<std::string> Compname;Compname.push_back("Carbon"); 
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
      //  std::cout<<i<<" "<<j<<" "<<l->nCells()<<" field "<<l->nFieldWires()<<std::endl;
      if(!senseWire&&l->nCells())
	senseWire=l->getCell(0)->getWire()->getDetail().get();
      if(!fieldWire&&l->nFieldWires())
	fieldWire=l->getFWire(0)->getDetail().get();
      
    }
  }
  std::cout<<"number of sense wires = "<<nsenseWires<<" field wires = "<<nfieldWires<<std::endl;
  std::map<std::string,double> materialArea;
  std::vector<std::string> matname=senseWire->materialNames();
  std::vector<double> thick=senseWire->shellsThicknesses();
  double radius=0;
  for(unsigned int i=0;i<matname.size();i++){
    std::cout<<"sense material name = "<<matname[i]<<" dr,mm = "<<thick[i]<<std::endl;
    materialArea[matname[i]]+=TMath::Pi()*((thick[i]+radius)*(thick[i]+radius)-radius*radius)*nsenseWires;
    radius+=thick[i];
  }
  matname=fieldWire->materialNames();
  thick=fieldWire->shellsThicknesses();
  radius=0;
  for(unsigned int i=0;i<matname.size();i++){
    std::cout<<"field material name = "<<matname[i]<<" dr,mm = "<<thick[i]<<std::endl;
    materialArea[matname[i]]+=TMath::Pi()*((thick[i]+radius)*(thick[i]+radius)-radius*radius)*nfieldWires;
    radius+=thick[i];
  }

  double totweigth=0;
  std::map<std::string,double>::iterator itm=materialArea.begin();
  for(;itm!=materialArea.end();itm++){
    G4Material* mat=mu2e::findMaterialOrThrow((*itm).first);
    (*itm).second*=mat->GetDensity()/(g/cm/cm/cm);
    totweigth+=(*itm).second;
    std::cout<<"material "<<mat->GetName()<<" density,g/cm^3 = "<<mat->GetDensity()/(g/cm/cm/cm)<<std::endl;
  }

  double gasarea=TMath::Pi()*(rz[1][0]*rz[1][0]-rz[0][1]*rz[0][1]);
  std::string gasmat=itr->getSuperLayer(0)->getLayer(0)->getDetail()->materialName();
  G4Material* mat=mu2e::findMaterialOrThrow(gasmat);
  const G4double* elemfract=mat->GetFractionVector();
  
  std::cout<<"gas "<<mat->GetName()<<" density,g/cm^3 = "<<mat->GetDensity()/(g/cm/cm/cm)<<std::endl;
 
  for(unsigned int i=0;i<mat->GetNumberOfElements();i++){
    const G4Element* elem=mat->GetElement(i);
    materialArea[elem->GetName()]=elemfract[i]*gasarea*mat->GetDensity()/(g/cm/cm/cm);//*mat->GetPressure()/(1.0*CLHEP::atmosphere);
    totweigth+=materialArea[elem->GetName()];
  }
  
  itm=materialArea.begin();
  for(;itm!=materialArea.end();itm++){
    std::cout<<"material "<<(*itm).first<<" has weigth fraction "<<(*itm).second/totweigth<<std::endl;
  }
  std::cout<<"total average gas+wire density,g/cm^3 = "<<totweigth/gasarea<<std::endl;
  std::map<std::string,std::string> matnamemap;
  matnamemap["G4_Mo"]="Molybdenum";
  matnamemap["G4_Ag"]="Silver";
  matnamemap["G4_Al"]="Aluminum";
  matnamemap["G4_W"]="Tungsten";
  matnamemap["G4_Au"]="Gold";
  matnamemap["H"]="Hydrogen";
  matnamemap["C"]="Carbon";
  matnamemap["F"]="Fluorine";
  matnamemap["He"]="Helium";
  matnamemap["WAGVacuum"]="Hydrogen";

  std::string *ptrgasname=new std::string("ITgasAuto");
  density = totweigth/gasarea; // g/cm^3
  radlen = -10;      // g/cm^2
  temperature = 20.;  // deg
  pressure = 1.;      // atm
  state="gas"; // solid or gas
  nbrcomp=materialArea.size(); // + plus by weigth; - by atoms 
  Compflag.clear();Compweight.clear();Compname.clear();
  itm=materialArea.begin();
  for(;itm!=materialArea.end();itm++){
    Compflag.push_back(0);
    Compweight.push_back((*itm).second/totweigth);
    std::cout<<"G4 replace "<<(*itm).first<<" by MatEnv "<<matnamemap[(*itm).first]<<std::endl;
    Compname.push_back(matnamemap[(*itm).first]);
  }
  Obj = new MatMaterialObj(*ptrgasname, density, zeff, 
			   aeff, nbrcomp, &Compflag, &Compweight, 
			   &Compname, radlen, intlen, refindex, 
			   temperature, pressure, state);
  (*matobj)[ptrgasname]=Obj;
  std::cout<<"My materials"<<std::endl;
  MatDBInfo* mtdbinfo=new MatDBInfo;
  mtdbinfo->findDetMaterial(*ptrname)->printAll(std::cout);
  mtdbinfo->findDetMaterial(*ptrgasname)->printAll(std::cout);

  //print in format of babar material file
  std::vector<MatMaterialObj*> vecmat;
  vecmat.push_back((*matobj)[ptrname]);
  vecmat.push_back((*matobj)[ptrgasname]);
  MatMaterialList matlist(vecmat);
  matlist.print();
  matlist.getMaterialVector()->clear();//to protect from deleting

  // mtdbinfo->findDetMaterial("IT-Swire")->printAll(std::cout);
  // mtdbinfo->findDetMaterial("IT-Fwire")->printAll(std::cout);
  // mtdbinfo->findDetMaterial("IT-gas0")->printAll(std::cout);
}
