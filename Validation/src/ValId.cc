
#include "Validation/inc/ValId.hh"

int mu2e::ValId::declare(art::TFileDirectory tfs, 
			 std::string name, std::string title) {
  _hid = tfs.make<TH1D>( name.c_str(), title.c_str(), 121, -60.5, 60.5);
  return 0;
}

int mu2e::ValId::fill(int id) {

  int idc = compress(id);
  _hid->Fill(idc); 

  return idc;
}

int mu2e::ValId::compress(int id) {

  //int jj = id%10;
  int q1 = (id/10)%10;
  int q2 = (id/100)%10;
  int q3 = (id/1000)%10;
  int rr = id/10000;
  int code = 50; // unknown/other

  int aid = abs(id);
  int sign = (aid==0 ? 1 : id/aid);
  if(aid<30) {
    code = id;
  } else if(aid==111) { // pi0
    code = 30;
  } else if(aid==211) { // pi+/-
    code = sign*31;
  } else if(aid==130) { // K_L^0
    code = sign*32;
  } else if(aid==310) { // K_S^0
    code = sign*33;
  } else if(aid==321) { // K+/-
    code = sign*34;
    } else if(aid==2112) { // n
    code = 40;
  } else if(aid==2212) { // p+/-
    code = sign*41;
  } else if(aid>1000000000) { // excited nuclei
    code = sign*51;
  } else if(q1>0 && q2>0 && q3==0 && rr<1000) { // meson
    int i = q2;
    if(i<1) i = 1;
    if(i>5) i = 6;
    code = sign*(31+i);
  } else if(q1>0 && q2>0 && q3>0 && rr<1000) { // baryon
    int i = q3;
    if(i<1) i = 1;
    if(i>5) i = 6;
    code = sign*(41+i);
  }
  return code;
}
