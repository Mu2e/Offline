//
// A utility program that does two things:
// 1) reads the scons dependency dump and summarizes  
//    dependcies between directories (a.k.a packages), write out the summary
//    The output looks like:
// HDR BFieldGeom BFieldGeom GeneralUtilities Mu2eInterfaces 
// LIB BFieldGeom 
// PRD BFieldGeom boost cetlib cetlib_except clhep gsl messagefacility 
// ...
//     where HDR means Mu2e headers, LIB is Mu2e libraries, and
//     PRD means products dependencies (header or lib).
//     
// 2) read in the above summary file and make a list of packages
//    that depend on a given packge.  This is a what you need to 
//    warn the developer about in the partial checkout
//
// deps -h  shows help
//

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <numeric>
#include <boost/algorithm/string.hpp>


using namespace std;
enum type{header,lib,product,ntype};


class gnode {
public:
  gnode(string x):name(x) {}
  string name;
  // things pointing in, parents, these depend on this node
  set<string> in; 
  // things pointing out, children, node depends on these
  set<string> out; 
};
vector<gnode> gnodes;


class package {
public:
  string name;
  set<string> headers;
  set<string> libs;
  set<string> products;
  void clear() {
    name.clear();
    libs.clear();
    headers.clear();
    products.clear();
  }
  void print() const {
    cout << "HDR " << name << " " ;
    for(auto const& x : headers) cout << x << " ";
    cout << endl;
    cout << "LIB " << name << " " ;
    for(auto const& x : libs) cout << x << " ";
    cout << endl;
    cout << "PRD " << name << " " ;
    for(auto const& x : products) cout << x << " ";
    cout << endl;
  }
};


vector<package> packages;

void usage() {
  cout << 
    "\n"
    "    deps [OPTIONS]\n"
    "    -i  make tree summary, reading from cin, print to cout\n"
    "    -h PKG  print packages that depend on PKG via headers\n"
    "    -f FILENAME  the input dependency summary\n"
    "    \n"
    "     For example, to make a dependency summary file:\n"
    "        scons -Q --tree=prune | deps -i > deps.txt\n"
    "     Then to read that summary and print what packages depend on\n"
    "     a given package:\n"
    "        deps -h ConfigTools -f deps.txt\n"
    "\n";
}

/********************************************************************/
// parse the scons dependency tree
// dependencies are saved in vector "packages"
//
// The tree looks like
  /*
  | +-Alignment/src
  |   +-Alignment/src/AlignmentMap.cc
  |   +-Alignment/src/AlignmentMap.os
  |   | +-Alignment/src/AlignmentMap.cc
  |   | +-ConfigTools/inc/SimpleConfig.hh
  |   | +-/cvmfs/mu2e.opensciencegrid.org/artexternals/cetlib_except/v1_01_07/in
clude/cetlib_except/exception.h
  | | | +-lib/libmu2e_GeneralUtilities.so
  | +-lib/libmu2e_Alignment_AlignmentService_service.so
  | | +-[lib/libmu2e_Alignment.so]
  */

int maketree(istream& inp) {

  string line;
  int n = 0;
  package pkg;
  string currPkg;
  vector<string> strs;
  int ilib = -1;

  while(getline(inp,line)) {
    size_t i = line.find("+-"); // sense the line depth
    size_t bracket = line.find("[");
    if( bracket!=string::npos ) line.erase(bracket,1);
    bracket = line.find("]");
    if( bracket!=string::npos ) line.erase(bracket,1);

    if(line.find("scons:")!=string::npos) { //skip first line
    } else if (i==0) { // skip the root of the tree
    } else if (i==2) { // top level dir
      string pkgName = line.substr(4,line.size()-4);
      bool skip = currPkg=="bin" || currPkg=="tmp" || 
	currPkg=="lib" || currPkg=="SConstruct";
      if(!pkg.name.empty() && !skip) {
	packages.push_back(pkg);
      }
      pkg.clear();
      pkg.name = pkgName;
      currPkg = pkgName;
      //cout << pkgName << endl;
    } else if (currPkg=="bin") { // skipping
    } else if (currPkg=="tmp") { // skipping
    } else if (currPkg=="lib") { // note library level dependecies of a dir
      if(i==4) { // start new library
	boost::split(strs,line,boost::is_any_of("_."));
	string libPkg = strs[1];
	bool found = false;
	for(size_t i = 0; i< packages.size(); i++) {
	  if(packages[i].name==libPkg) {
	    ilib = i;
	    found = true;
	  }
	}
	if(!found) {
	  cout << "*************  failed to find lib "<< libPkg<< endl;
	  return 1;
	}
      } else { // amongst the dependencies for a library
	// if a package dependence
	if(line.find("/cvmfs")!=string::npos) {
	  boost::split(strs,line,boost::is_any_of("/"));
	  //cout << "insert " << strs[4] << " " << pkg << endl;
	  if(strs[4]!="gcc") {
	    packages[ilib].libs.insert(strs[4]);
	  }
	}
	// if a library dependence
	if(line.find("lib/libmu2e")!=string::npos) {
	  size_t i = line.find("_");
	  size_t j = line.find(".");
	  string libn = line.substr(i+1,j-i-1);
	  packages[ilib].libs.insert(libn);
	}
      }
    } else if (i>2) { // depth is below the main directory level
      auto ss = line.substr(i+2,line.size()-i-2);
      if(ss.find("/cvmfs")==0) { // must be product
	ss = ss.substr(45,ss.size()-45);
	size_t j = ss.find("/");
	if(j==string::npos) j = ss.size();
	ss = ss.substr(0,j);
	if(ss!="gcc") pkg.products.insert(ss);
      } else if(ss.find("lib")==0) { // alibrary dependence
	ss = ss.substr(4,ss.size()-4);
	pkg.libs.insert(ss);
      } else if(ss.find("tmp")==0) {  // a temp dependence - ignore
      } else if(ss.find("/usr")==0) { // system dep - ignore
      } else { // a subdirectory - a mu2e package dep
	size_t j = ss.find("/");
	if(j==string::npos) j = ss.size();
	ss = ss.substr(0,j);
	pkg.headers.insert(ss);
      }
    } else {
      cout << "bad line" << line << endl;
      return 1;
    }
    n++;
  }

  bool skip = pkg.name=="bin" || pkg.name=="tmp" || pkg.name=="lib";
  if(!skip) packages.push_back(pkg);

  // this prints the results to cout
  // good for saving a deps.txt
  for(auto const& x : packages) x.print();

  return 0;
}


/**********************************************************/
// find a node in the tree
// return maxuint if failed
size_t findNode(string ss) {

  for(size_t i=0; i<gnodes.size(); i++) {
    if(gnodes[i].name==ss) {
      return i;
    }
  }
  return std::numeric_limits<size_t>::max();
}

/**********************************************************/
// insert if needed, return node number
size_t insertNode(string ss) {
  size_t np = findNode(ss);
  if(np!=std::numeric_limits<size_t>::max()) {
    return np;
  }
  gnodes.emplace_back(ss);
  return gnodes.size()-1;
}

/**********************************************************/
// print the tree
int printNodes() {
  for(auto const& nn: gnodes) {
    cout << nn.name << " : " ;
    for(auto const& oo: nn.out) {
      cout << oo << " " ;
    }
    cout << endl;
    cout << nn.name << " :: " ;
    for(auto const& ii: nn.in) {
      cout << ii << " " ;
    }
    cout << endl;
  }
  return 0;
}

/**********************************************************/
// read deps.txt, the format written by "maketree"
// only saves the header dependencies, but could 
// also save library and product deps

int readGraph( string fn ) {

  ifstream inp(fn);
  string line;
  vector<string> strs;

  while(getline(inp,line)) {
    boost::split(strs,line,boost::is_any_of(" "));
    if(strs[0]=="HDR") { // if this line lists header dependecies
      size_t np = insertNode(strs[1]); // np = the main directory
      for(size_t i=2; i<strs.size(); i++) { // the dirs main dir depends on
	size_t mp = insertNode(strs[i]); // insert if neeed
	// add it to this node's dep
	if (mp!=np) { // don't put in a self-reference
	  gnodes[np].out.insert(strs[i]); // out means np depends on it
	  gnodes[mp].in.insert(strs[1]);
	} // not the same
      } // loop over dependencies in a line
    } // if HDR
  } // input lines

  return 0;
}

/**********************************************************/
// print what depends on pkgName in type tt
// note that the tree is not reduced,
// every node has a full list of its direct dependencies, 
// both in and out, but nodes do not have links
// to their transitive dependencies

int printDep(type tt, string pkgName) {

  size_t np = findNode(pkgName);
  if(np<0) {
    cout << "package "<<pkgName<<" is not in the dependency list " << endl;
    return 1;
  }

  vector<string> ss;
  for(auto name: gnodes[np].in) {
    ss.emplace_back(name);
  } // loop over in

  for(auto s : ss) {
    cout << s << " ";
  }
  if(ss.size()>0) cout << endl;

  return 0;
}


/**********************************************************/

int main( int argc, char** argv ) {
  
   int opt;
   string fn;
   string pkgName;
   ifstream inp;
   type tt = header;

   while ((opt = getopt(argc, argv, "is:f:h:l:p:")) != -1) {
     switch (opt) {
     case 'i':
       return maketree(cin);
     case 's':
       inp = ifstream(optarg);
       return maketree(inp);
     case 'f':
       fn = optarg;
       break;
     case 'h':
       pkgName = optarg;
       tt = header;
       break;
     case 'l':
       pkgName = optarg;
       tt = lib;
       break;
     case 'p':
       pkgName = optarg;
       tt = product;
       break;
     default: /* '?' */
       usage();
       return 1;
     }
   }

   if(fn=="") {
     cout << "ERROR - no deps file specified" << endl;
     return 1;
   }

   int rc;
   rc = readGraph(fn);
   if(rc!=0) return rc;
   rc = printDep(tt,pkgName);
   return rc;
}
