#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <stdexcept>
#include <list>
#include <set>
#include <fstream>

#include <boost/program_options.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>

typedef unsigned long long ullong;

int main(int argc, char** argv) {
  try {
    return 0;
  } catch (exception& e) {
    cerr << e.what() << endl;
    return -1;
  }
}
