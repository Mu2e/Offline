# root -q 'plotAllSpectra.C({"Stage1/S1EleVD.root"},  {"Stage1/S1MuVD.root", "Stage1/S1Mu3VD.root",  "Stage1/S11809VD.root"},  "Stage1virtualdetector/ttree", "virtualdetector", 42)'
# root -q 'plotAllSpectra.C({"Stage2/S2EleVD.root"},  {"Stage2/S2MuVD.root",  "Stage2/S21809VD.root"},  "Stage2virtualdetector/ttree", "virtualdetector")'
root -q 'plotAllSpectra.C({"Stage2/S2EleDet.root"}, {"Stage2/S2MuDet.root", "Stage2/S21809Det.root"}, "Stage2HPGe/ttree",        "HPGe")'
