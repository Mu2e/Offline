//
// A utility program that does two things:
// 1) reads the scons dependecny dump and summarizes it 
//    and dependcies between directories (packages), write out the summary
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
#include <map>
#include <unistd.h>
#include <algorithm>
#include <numeric>
#include <boost/algorithm/string.hpp>
//#include <boost/graph/adjacency_list.hpp>
//#include <boost/graph/transitive_reduction.hpp>

using namespace std;
enum type{header,lib,product,ntype};


class gnode {
public:
  gnode(string x):name(x) {}
  string name;
  set<string> in; // things pointing in, parents
  set<string> out; // thing spointing out, chlidren
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
    "    -h PKG  print packages that depend on PKG\n"
    "    -f FILENAME  the input dependency summary\n"
    "    \n"
    "     For example, to make a dependency summary file:\n"
    "        scons -Q --tree=prune | deps -i > deps.txt\n"
    "     Then to read that summary and print what packages depend on\n"
    "     a given package:\n"
    "        deps -h ConfigTools -f deps.txt\n"
    "\n";
}

int maketree(istream& inp) {
  //cout << "maketree " << std::endl;
  //auto inp = cin;
  //if(!fn.empty()) inp = ifstream(fn);

  string line;
  int n = 0;
  package pkg;
  string currPkg;
  vector<string> strs;
  int ilib = -1;

  while(getline(inp,line)) {
    //if(n<10) {
    // cout << line << endl;
    //}
    size_t i = line.find("+-");
    //cout << i << endl;
    size_t bracket = line.find("[");
    if( bracket!=string::npos ) line.erase(bracket,1);
    bracket = line.find("]");
    if( bracket!=string::npos ) line.erase(bracket,1);

    if(line.find("scons:")!=string::npos) {
    } else if (i==0) {
    } else if (i==2) {
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
    } else if (currPkg=="bin") {
    } else if (currPkg=="tmp") {
    } else if (currPkg=="lib") {
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
	}
	//cout << "new lib "<<libPkg << " " <<found << endl;
      } else { // in a stanza for a library
	if(line.find("/cvmfs")!=string::npos) {
	  boost::split(strs,line,boost::is_any_of("/"));
	  //cout << "insert " << strs[4] << " " << pkg << endl;
	  if(strs[4]!="gcc") {
	    packages[ilib].libs.insert(strs[4]);
	  }
	}
	if(line.find("lib/libmu2e")!=string::npos) {
	  size_t i = line.find("_");
	  size_t j = line.find(".");
	  string libn = line.substr(i+1,j-i-1);
	  packages[ilib].libs.insert(libn);
	}
      }
    } else if (i>2) {
      auto ss = line.substr(i+2,line.size()-i-2);
      if(ss.find("/cvmfs")==0) { // must be product
	ss = ss.substr(45,ss.size()-45);
	size_t j = ss.find("/");
	if(j==string::npos) j = ss.size();
	ss = ss.substr(0,j);
	if(ss!="gcc") pkg.products.insert(ss);
      } else if(ss.find("lib")==0) {
	cout << "******************** lib in code" <<endl;
	ss = ss.substr(4,ss.size()-4);
	pkg.libs.insert(ss);
      } else if(ss.find("tmp")==0) {
      } else if(ss.find("/usr")==0) {
      } else { // a subdirectory
	size_t j = ss.find("/");
	//cout << "subdir " << ss<< endl;
	if(j==string::npos) j = ss.size();
	ss = ss.substr(0,j);
	//cout << "subdir " << ss << " " << j<< endl;
	pkg.headers.insert(ss);
      }
    } else {
      cout << "bad line" << line << endl;
    }
    n++;
    //if(n>1000) return 1;
  }

  bool skip = pkg.name=="bin" || pkg.name=="tmp" || pkg.name=="lib";
  if(!skip) packages.push_back(pkg);

  if(true) {
    for(auto const& x : packages) x.print();
  }

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
  return 0;
}


/**********************************************************/
size_t findNode(string ss) {

  for(size_t i=0; i<gnodes.size(); i++) {
    if(gnodes[i].name==ss) {
      return i;
    }
  }
  return std::numeric_limits<size_t>::max();
}

/**********************************************************/
size_t insertNode(string ss) {
  size_t np = findNode(ss);
  if(np!=std::numeric_limits<size_t>::max()) {
    return np;
  }
  gnodes.emplace_back(ss);
  return gnodes.size()-1;
}

/**********************************************************/

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

int readGraph( string fn ) {

  //cout << "open " << fn<<endl;
  ifstream inp(fn);
  string line;
  vector<string> strs;
  //cout << "stat " << inp.is_open() << endl;

  while(getline(inp,line)) {
    //cout << "read " << line<<endl;
    boost::split(strs,line,boost::is_any_of(" "));
    if(strs[0]=="HDR") {
      //cout << "working on " << strs[1]<<endl;
      size_t np = insertNode(strs[1]);
      for(size_t i=2; i<strs.size(); i++) {
	size_t mp = insertNode(strs[i]);
	// add it to this node's dep
	gnodes[np].out.insert(strs[i]);
	gnodes[mp].in.insert(strs[1]);
      } // loop over dependencies in a line
    } // if HDR
  } // input lines

    //print_graph(tr);
  //printNodes();

  return 0;
}

/**********************************************************/
int printDep(type tt,string pkgName) {
  //  cout << "printDep " << tt << " " << pkgName << endl;
  size_t np = findNode(pkgName);
  if(np<0) {
    cout << "package "<<pkgName<<" is not in the dependency list " << endl;
  }

  map<string,int> sr;
  for(auto const gn : gnodes) sr[gn.name] = 0;
  sr[pkgName] = 1;
  bool more = true;
  while(more) {
    more = false;
    for(auto si: sr) {
      if(si.second==1) {
	size_t np = findNode(si.first);
	for(auto const& ii: gnodes[np].in) {
	  if(ii!=si.first) { // skip the one we are working on
	    auto stat = sr[ii];
	    if(stat==0) {
	      sr[ii] = 1;
	      more=true;
	    }
	  }
	} // loop over in
      } // si.second==1
      si.second = 2;
    } // si : sr
  }

  for(auto si: sr) {
    if(si.first!=pkgName) {
      if(si.second>0) cout << si.first << " ";
    }
  }
  cout << endl;
  return 0;
}


///**********************************************************/
//int printDep(type tt,string pkgName) {
//  cout << "printDep " << tt << " " << pkgName << endl;
//  size_t np = findNode(pkgName);
//  if(np<0) {
//    cout << "package "<<pkgName<<" is not in the dependency list " << endl;
//  }
//
//  set<string> parents;
//  parents.insert(pkgname);
//
//  for(auto const gn : gnodes) {
//  }
//  return 0;
//}

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

   readGraph(fn);
   printDep(tt,pkgName);

  return 0;
}
