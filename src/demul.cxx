#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <stdexcept>

#include <boost/program_options.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>

typedef unsigned long long ullong;

using namespace std;

vector<string> load_index_file(string filename) throw (exception) {
  iostream fi(filename.c_str());
  if (!fi.is_open()){
    throw runtime_error(string("cannot open ") + filename);
  }
  while (!fi.eof()) {
    string line;
    list<string> items;
    getline(fi, line);
    boost::regex acgt("^[ACGT]\\b");
    boost::smatch matches;
    boost::split(items, line, boost::is_space());
    if (boost::regex_match(line, matches, acgt)) {
      
    }
  }
  
  fi.close();
}
