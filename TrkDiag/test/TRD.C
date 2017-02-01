#include "../../Scripts/FillChain.C"
#include "TrkDiag/test/TrkRecoDiag.C+"
void TRD() {
  TChain* mtce(0);
  FillChain(mtce,"MTCE.txt","TrkRecoDiag/trdiag");
  TChain* mtbk(0);
  FillChain(mtbk,"bkgfiles.txt","TrkRecoDiag/trdiag");
  TrkRecoDiag trdce(mtce);
  trdce.Loop();
  trdce.drawHistos();
  TrkRecoDiag trdbk(mtbk);
  trdbk.Loop();
  trdbk.drawHistos();




}
