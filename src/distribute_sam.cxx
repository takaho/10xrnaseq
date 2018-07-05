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

using namespace std;
/*
first_indices = {
    'AACGTCAA':0,#	85457155
    'GCATCTCC':1,#	84432940
    'CTGCGATG':2,#	70279153
    'TGTAAGGT':3,#	54412965

    'GTCCTAAC':4,#	83802878
    'CGGTCCCA':5,#	78867120
    'ACAGAGGT':6,#	73400348
     NCAGAGGT
    'TATAGTTG':7}#	67410698
*/

namespace {
  // string assume_original(const map<string,int>& indices, int num) {
  //   vector<map<char,int>> freq;
  //   for (auto it : indices) {
  //     if ( it.second == num ) {
  // 	string seq = it.first;
  // 	if (seq.size() >= freq.size()) {
  // 	  for (auto i = freq.size(); i < seq.size(); i++) {
  // 	    freq.push_back(map<char,int>());
  // 	  }
  // 	}
  // 	for (auto i = 0; i < (int)seq.size(); i++) {
  // 	  char c = seq.c_str()[i];
  // 	  if (freq[i].find(c) == freq[i].end()) {
  // 	    freq[i][c] = 1;
  // 	  } else {
  // 	    freq[i][c] ++;
  // 	  }
  // 	}
  //     }
  //   }
  //   int N = freq.size();
  //   char* buf = new char[N + 1];
  //   buf[N] = '\0';
  //   for (int i = 0; i < N; i++) {
  //     int maxnum = 0;
  //     char base = 'X';
  //     for (auto it : freq[i]) {
  // 	if (it.second > maxnum) {
  // 	  base = it.first;
  // 	  maxnum = it.second;
  // 	}
  // 	buf[i] = base;
  //     }
  //   }
  //   string tag(buf);
  //   delete[] buf;
  //   return tag;
  // }
  string convert_from_code(ullong code, int length) {
    char *buf = new char[length + 1];
    buf[length] = '\0';
    for (int i = 0; i < length; i++) {
      ullong n = (code >> ((length - 1 - i) * 2)) & 0x03;
      switch (n) {
      case 0:
	buf[i] = 'A'; break;
      case 1:
	buf[i] = 'C'; break;
      case 2:
	buf[i] = 'G'; break;
      default:
	buf[i] = 'T'; break;
      }
    }
    string seq = buf;
    delete[] buf;
    return seq;
  }
  ullong convert_to_code(const string& line, int length=0) {
    const char* ptr = line.c_str();
    if (length == 0) {
      length = line.size();
    }
    ullong code = 0;
    for (int i = 0; i < (int)line.size(); i++) {
      char c = ptr[i];
      ullong c_ = 0;
      switch (c) {
      case 'A': c_ = 0; break;
      case 'C': c_ = 1; break;
      case 'G': c_ = 2; break;
      case 'T': c_ = 3; break;
      default :
	c_ = 0; break;
      }
      code = (code << 2) | c_;
    }
    return code;
  }
  
  
  void _substitute_recursively(const char* seq, char* buffer, int mm, int pos, int span, vector<string>& seeds) {
    if (mm <= 0 || pos >= span) return;
    for (char base : {'A', 'C', 'G', 'T', 'N'}) {
      if (seq[pos] != base) {
	buffer[pos] = base;
	//cout << seq << " " << pos << " " << buffer << endl;
	seeds.push_back(buffer);
	_substitute_recursively(seq, buffer, mm - 1, pos + 1, span, seeds);
	buffer[pos] = seq[pos];
      } else {
	_substitute_recursively(seq, buffer, mm, pos + 1, span, seeds);
      }
    }
  }
  
  vector<string> generate_substituted(const string& seq, int mismatches) {
    vector<string> sub;
    char* buffer = new char[seq.size() + 1];
    memcpy(buffer, seq.c_str(), seq.size() + 1);
    _substitute_recursively(seq.c_str(), buffer, mismatches, 0, seq.size(), sub);
    delete[] buffer;
    return sub;
  }
}

map<string,int> load_well_indices(const string& filename_tag) {
  ifstream fi(filename_tag.c_str());
  map<string,int> codes;
  if (fi.is_open() == false) {
    codes["AACGTCAA"]=0;
    codes["GCATCTCC"]=1;
    codes["CTGCGATG"]=2;
    codes["TGTAAGGT"]=3;
    codes["GTCCTAAC"]=4;
    codes["CGGTCCCA"]=5;
    codes["ACAGAGGT"]=6;
    codes["TATAGTTG"]=7;
    //return codes;
  } else {
    boost::regex acgt("^[ACGTN]{4,}$");
    while (!fi.eof()) {
      string line;
      list<string> items;
      //vector<string> results;
      //boost::match_results<string> results;
      boost::smatch matches;
      getline(fi, line);
      boost::split(items, line, boost::is_any_of("\t ,"));
      if (boost::regex_match(line, matches, acgt)) {
	cout << matches[0] << endl;
	//if (boost::regex_search(line.begin(), line.end(), results, acgt)) {
	// cout << results[0] << endl;
	if (codes.find(matches[0]) == codes.end()) {
	  codes[matches[0]] = codes.size();
	}
      }
    }
  }
  return codes;
}

void expand_indices(map<string,int>& codes, int mismatch_tolerance=1) {
// }
// map<string,int> load_well_indices(const string& filename_tag, int mismatch_tolerance=1) {
  if (mismatch_tolerance > 0) {
    map<string,int> additional;
    set<string> ambiguous;
    for (auto&& it : codes) {
      const string& seq = it.first;
      int num = it.second;
      for (int j = 0; j < (int)seq.size(); j++) {
	//char c = seq.c_str()[j];
	vector<string> seeds = generate_substituted(seq, mismatch_tolerance);
	//cout << it.second << " " << seq << endl;
	for (auto s : seeds) {
	  if (additional.find(s) != additional.end() && additional[s] != num) {
	    ambiguous.insert(s);
	  }
	  //	  if (additional.find(s) != additional.end()) ambiguous.insert(s);
	  additional[s] = num;
	  //cout << s << " => " << num << endl;
	}
      }
    }
    for (auto&& it : codes ) {
      cout << it.first << " " << it.second;
      char delim = '\t';
      for (auto&& itr : additional) {
	if (itr.second == it.second) {
	  cout << delim << itr.first;
	  if (ambiguous.find(itr.first) != ambiguous.end()) {
	    cout << "*";
	  }
	  delim = ',';
	}
      }
      cout << endl;
    }
    
    for (auto&& it : additional) {
      if (ambiguous.find(it.first) == ambiguous.end()) {
	codes[it.first] = it.second;
      }
    }
  }
  // for (map<string,int>::const_iterator seq = codes.begin(); seq != codes.end(); seq++) {
  //   cout << seq->first << " " << seq->second << endl;
  // }
  //return codes;
}

int main(int argc, char** argv) {
  namespace po = boost::program_options;
  namespace io = boost::iostreams;

  po::options_description opt("options");
  string filename_sam;
  string filename_index;
  string filename_tag;
  map<string,int> wellcode;
  //map<string,int> wellcode_strict;
  ullong range_max = 0;
  ullong range_min = 0;
  int lap;
  int shift = 60;
  int barcode_size = 16;
  opt.add_options()
    ("tags,t", po::value<string>(&filename_tag)->default_value(""), "tag file")
    ("samfile,s", po::value<string>(&filename_sam), "sam file")
    ("index_file,I", po::value<string>(&filename_index), "fastq.gz")
    ("range_max", po::value<ullong>(&range_max)->default_value(0), "Index range")
    ("range_min", po::value<ullong>(&range_min)->default_value(0), "Index range")
    ("barcode_size,b", po::value<int>(&barcode_size)->default_value(16), "Barcode length (16)")
    ("lap,L", po::value<int>(&lap)->default_value(0), "log period")
  ("help,h", "display help");

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, opt), vm);
  } catch (exception& e) {
    cerr << opt << endl;
    cerr << e.what() << endl;
    return -1;
  }
  po::notify(vm);
  //cout << opt << endl;
  //cout << vm << endl;
  if (vm.count("help")) {
    cout << opt << endl;
  } else {
    try {
      //cerr << "well barcode\n";
      wellcode = load_well_indices(filename_tag);
      map<int,string> rev_wellcode;
      for (auto it: wellcode) {
	rev_wellcode[it.second] = it.first;
	//cout << it.second << " " << it.first << endl;
      }
      expand_indices(wellcode, 1);
      //cout << "expanded" << endl;
      //cerr << "well barcode\n";
      // for (map<string,int>::const_iterator it = wellcode.begin(); it != wellcode.end();it++) {
      // 	cout << it->first << " " << it->second << endl;
      // }
      
      // cerr << "1\n";
      // const string filename_sam = vm["s"];//.as<string>();
      // cerr << "2\n";
      // string filename_index = vm["I"];//.as<string>();
      // cerr << "3\n";
      //cout << filename_index << endl;
      //      exit(0);
      io::filtering_istream fistr;// io::filtering_istream();
      fistr.push(io::gzip_decompressor());
      fistr.push(io::file_source(filename_index));
      list<string> items;
      //istream fi;//*in = fltr_is;
      map<ullong,int> freq;
      ullong num_lines = 0;
      ullong num_reads = 0;
      int indexcode = -1;
      while (true) {
	//for (int i = 0; i < 10; i++) {
	string line;
	getline(fistr, line);
	//cout << num_lines << " " << line << endl;
	size_t flag = num_lines & 3;
	//cout << num_lines << " " << flag << " " << line << endl;
	if (flag == 0) {
	  //switch ((num_lines & 3)) {
	  //case 0:
	  boost::split(items, line, boost::is_any_of(":"));//boost::is_space());
	  map<string,int>::const_iterator it = wellcode.find(*items.rbegin());
	  if (it != wellcode.end()) {
	    indexcode = it->second;
	    //cout << "INDEX " << indexcode << " " << *(items.rbegin()) << endl;
	  }
	  //	  break;
	  //case 1:
	} else if (flag == 1) {
	  //cout << num_lines << " " << indexcode << endl; // 
	  if (indexcode >= 0) {
	    ullong nuccode = convert_to_code(line, barcode_size) | ((ullong)indexcode << shift);
	    //cout << line << "=> " << nuccode << " " << indexcode << endl;
	    if (range_max <= 0 || (range_min <= nuccode && nuccode < range_max)) {
	      map<ullong,int>::iterator it = freq.find(nuccode);
	      if (it == freq.end()) {
		//cout << "NEW " << nuccode << endl;
		freq[nuccode] = 1;
	      } else {
		it->second++;
		//cout << "ADD " << freq[nuccode] << endl;
	      }
	    }
	    num_reads++;
	  }
	//   break;
	// default:
	//   break;
	} else if (flag == 3) {
	  indexcode = -1;
	}
	num_lines ++;
	if (lap > 0 && num_lines % lap == 0) {
	  vector<ullong> keys;
	  for (auto&& it : freq) {
	    keys.push_back(it.first);
	  }
	  sort(keys.begin(), keys.end(), [&](ullong a, ullong b) { return freq[a] > freq[b]; });//const string& k_, const string& l_) {
	  //});
	  // sort(keys.begin(), keys.end(), [=](const string& k_, const string& l_) {
	  //     return freq[k_] < freq[l_];
	  //   } );
	  //cout << freq.size() << endl;

	  ofstream log("log.txt");
	  log << "#num_reads:" << num_reads << endl;
	  log << "#num_lines:" << num_lines << endl;
	  log << "#filename:" << filename_index << endl;
	  cout << num_lines << " ; " << freq.size() << endl;
	  int num_disp = 0;
	  for (auto&& k_ : keys) {
	    int num = freq[k_];
	    //if (num >= 10) {
	      ullong group = (k_ >> shift) & 0x0f;
	      string tag = convert_from_code((k_ ^ (group << shift)), barcode_size);
	      log << rev_wellcode[group] << "\t" << tag << "\t" << num << endl;
	      if (++num_disp <= 10) {
		cout << k_ << " " << rev_wellcode[group] << " " << tag << " " << num << endl;
	      }
	      //}
	  }
	  log.close();
	}
      }
      //fistr.close();
      //delete fltr_is;
    } catch (exception& e) {
      cerr << e.what() << endl;
    }
  }
  return 0;
}
