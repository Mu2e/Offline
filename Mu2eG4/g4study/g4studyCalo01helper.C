{
TTree*   _sttree = 0x0;
_file0->GetObject("checkhits/tstep", _sttree);
_sttree->Process("Mu2eG4/g4study/g4studyCalo01Selector.C+g");
}
