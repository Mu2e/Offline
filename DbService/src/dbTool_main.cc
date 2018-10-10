
#include <string>
#include <vector>
#include <iostream>
#include "DbService/inc/DbTool.hh"

int main(int argc, char**argv) {

  std::vector<std::string> words;
  for(size_t i=1; i<size_t(argc);++i) {
    //std::cout << argv[i] <<std::endl;
    words.emplace_back(argv[i]);
  }

  mu2e::DbTool tool;
  int rc;
  rc = tool.setArgs(words);
  if(rc!=0) return rc;
  rc = tool.init();
  if(rc!=0) return rc;
  rc = tool.run();
  if(rc!=0) return rc;

  return 0;
}
