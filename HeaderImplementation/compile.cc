{
	gROOT->ProcessLine(".L FindPeakBaseRoot.cc+");
	gROOT->ProcessLine(".L FitModel.cc+");
	gROOT->ProcessLine(".L FitModelRoot.cc+");
	gROOT->ProcessLine(".L FindMultiplePeak.cc+");

	gROOT->ProcessLine(".L FindPeakBase.hh+");
	gROOT->ProcessLine(".L FindPeakBaseRoot.hh+");
	gROOT->ProcessLine(".L FitModel.hh+");
	gROOT->ProcessLine(".L FitModelRoot.hh+");
	gROOT->ProcessLine(".L FindMultiplePeak.hh+");
	gROOT->ProcessLine(".L config.hh+");
	gROOT->ProcessLine(".L ParamStructs.hh+");

}