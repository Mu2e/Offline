{
  gROOT->LoadMacro("TrkPatRec/test/Timing.C+");
  std::vector<std::string> names;
  names.push_back("MakeStereoHits");
  names.push_back("FlagStrawHits");
  names.push_back("FlagBkgHits");
  names.push_back("TPRDownstreameMinus");
  Timing(names);
}

