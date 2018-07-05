#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <stdexcept>
#include <list>
#include <set>
#include <fstream>
#include <iomanip>

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
  
  pair<ullong,ullong> encode_sequence(const string& seq, int N=0) {
    ullong code = 0;
    ullong mask = 0;
    const char* ptr = seq.c_str();// + code.size() - 1;
    if (N == 0) {
      N = seq.size();
    }
    for (int i = 0; i < N; i++) {
      char c = ptr[i];
      int shift = (N - 1 - i) << 1;
      if (c == 'N') {
	mask |= (3 << shift);
      } else if (c == 'C') {
	code |= (1 << shift);
      } else if (c == 'G') {
	code |= (2 << shift);
      } else if (c == 'T') {
	code |= (3 << shift);
      }
    }
    return make_pair(code, mask);
  }

  vector<ullong> encode_indexes(const vector<string>& seqs, int length=0) {
    vector<ullong> codes;
    //cout << "encode requested " << seqs.size() << ", " << length << endl;
    for (auto it = seqs.begin(); it != seqs.end(); it++) {
      //cout << length << " " << *it << " " << it->size() << endl;
      if (length > 0 && it->size() != length) {
	//cout << "SKIP " << *it << " " << length << " " << it->size() << endl;
	continue;
      }
      pair<ullong,ullong> cm = encode_sequence(*it, length);
      if (cm.second == 0) {
	codes.push_back(cm.first);
      // } else {
      // 	cout << "N included " << *it << endl;
	continue;
      }
    }
    cout << codes.size() << " codes generated\n";
    return codes;
  }
}

namespace {
  map<string,ullong> reduce_cache_size(const map<string,ullong>& cache, ullong capacity) {
    if (capacity > cache.size()) {
      return cache;
    }
    map<string,ullong> reduced;
    ullong* buf = new ullong[cache.size()];
    size_t i = 0;
    for (auto it = cache.begin(); it != cache.end(); it++) {
      buf[i++] = it->second;
    }
    sort(buf, buf + cache.size());
    ullong threshold = buf[capacity];
    //cerr << cache.size() << " => " << threshold << "      ";
    delete[] buf;
    for (auto it = cache.begin(); it != cache.end(); it++) {
      if (it->second > threshold) {
	reduced[it->first] = it->second;
      }
    }

    //cerr << cache.size() << " => " << reduced.size() << "         \r";
    return reduced;
  }

  void display_frequency_stat(const map<string, ullong>& cache, ullong size, ostream& out) {
    out << "display frequency\n";
    vector<string> codes;
    for (auto it = cache.begin(); it != cache.end(); it++) {
      codes.push_back(it->first);
    }
    sort(codes.begin(), codes.end(), [&](const string& a, const string& b){return cache.find(a)->second > cache.find(b)->second;});
    ullong end = size < cache.size() ? size : cache.size();
    //out << end << ", " << codes.size() << endl;
    for (size_t i = 0; i < end; i++) {
      out << codes[i] << "\t" << cache.find(codes[i])->second << "\n";
    }
  }


  pair<map<string,ullong>,map<string,ullong>> load_index_sets(const string& filename_index, int well_barcode_size=8, int cell_barcode_size=16, ullong well_index_capacity=200,ullong cell_index_capacity=16000, ullong lap=1000000, bool verbose=false) throw (exception) {
    namespace io = boost::iostreams;

    map<string,ullong> well_indexes;
    map<string,ullong> cell_indexes;
    // count_freq
    io::filtering_istream fistr;// io::filtering_istream();
    if (filename_index.rfind(".gz") == filename_index.size() - 3) {
      fistr.push(io::gzip_decompressor());
      //cerr << ".gz file\n";
    }
    fistr.push(io::file_source(filename_index));
    fistr.set_auto_close(true);
    list<string> items;
    ullong num_line = 0;
    while (true) {
      string line;
      if (!getline(fistr, line)) {
	break;
      }
      //cout << num_line << " " << line << endl;
      size_t slot = num_line & 3;
      num_line ++;
      if (slot == 0) { // title
	size_t pos = line.rfind(':');
	if (pos != string::npos) {
	  string wellcode = line.substr(pos + 1, well_barcode_size);
	  auto it = well_indexes.find(wellcode);
	  if (it == well_indexes.end()) {
	    //cerr << (well_indexes.size() + 1) << " " << wellcode << " new well index ";
	    if (well_indexes.size() > (well_index_capacity * 4)) {
	      well_indexes = reduce_cache_size(well_indexes, well_index_capacity);
	    }
	    //cout << wellcode << " new well index\n";
	    well_indexes[wellcode] = 1;
	    //cout << well_indexes.size() << endl;
	  } else {
	    //cout << wellcode << " count up\n";
	    it->second ++;
	  }
	}
      } else if (slot == 1) {
	const string& code = line.substr(0, cell_barcode_size);
	//cout << code << " " << cell_barcode_size << " " << line << endl;
	auto it = cell_indexes.find(code);
	if (it == cell_indexes.end()) {
	  //cout << (cell_indexes.size() + 1) << " " << code << " new cell index\n";
	  if (cell_indexes.size() > (cell_index_capacity * 4) ) {
	    //cerr << cell_indexes.size() << " ? " << cell_index_capacity << " ? " << (cell_indexes.size() > (cell_index_capacity * 4)) << endl;
	    cell_indexes = reduce_cache_size(cell_indexes, cell_index_capacity);
	  }
	  cell_indexes[code] = 1;
	} else {
	  it->second++;
	}
      } else if (slot == 3) {
	if (lap > 0 && ((num_line >> 2) % lap) == 0) {
	  cerr << "#Well index : " << well_indexes.size() << " in " << (num_line / 4) << " reads\n";
	  display_frequency_stat(well_indexes, 8, cerr);
	  cerr << "#Cell index : " << cell_indexes.size() << "\n";
	  display_frequency_stat(cell_indexes, 10, cerr);
	  //cerr << "EXIT\n";
	}
      }
    }
    
    io::close(fistr);
    //fistr.close();
    if (verbose) {
      display_frequency_stat(well_indexes, 20, cerr);
      display_frequency_stat(cell_indexes, 20, cerr);
    }
    well_indexes = reduce_cache_size(well_indexes, well_index_capacity);
    cell_indexes = reduce_cache_size(cell_indexes, cell_index_capacity);
    return make_pair(well_indexes, cell_indexes);
  }
  vector<string> obtain_reduced_indexes(const map<string,ullong>& indexes, int num) {
    vector<string> seqs;
    for (auto it = indexes.begin(); it != indexes.end(); it++) {
      if (it->second > num) {
	seqs.push_back(it->first);
      }
    }
    return seqs;
  }

  string decode_index(int length, ullong code) {//, int offset=0) {
    char* buf = new char[length + 1];
    //buf[length] = buf[length + 1] = '\0';
    int index = 0;
    for (int i = 0; i < length; i++) {
      ullong c = (code >> ((length - 1 - i) << 1)) & 0x03;
      char base = 'A';
      if (c == 1) {
	base = 'C';
      } else if (c == 2) {
	base = 'G';
      } else if (c == 3) {
	base = 'T';
      }
      buf[index++] = base;
      // if (offset > 0 && offset == i) {
      // 	  buf[index++] = '_';
      // }
    }
    buf[index] = '\0';
    string name(buf);
    delete[] buf;
    return name;
  }

#define HAVE_X64INTRIN 1
#ifdef HAVE_X64INTRIN
#include <x86intrin.h>
  
  inline uint64_t popcount(uint64_t n) {
    uint64_t res;
    __asm__( "popcntq %1, %0" : "=r"(res) : "r"(n));
    return res;
  }
#else
#error cannot use popcount  
#endif
  void generate_renamed_fastq(const string& filename_index,
			      int well_code_length, const vector<ullong>& well_codes,
			      int cell_code_length, const vector<ullong>& cell_codes,
			      const string& filename_read, int mismatches, boost::iostreams::filtering_ostream& ostr) throw (exception) {
    
    namespace io = boost::iostreams;
    io::filtering_istream fistr_index;// io::filtering_istream();
    if (filename_index.rfind(".gz") == filename_index.size() - 3) {
      fistr_index.push(io::gzip_decompressor());
    }
    fistr_index.push(io::file_source(filename_index));
    fistr_index.set_auto_close(true);
    io::filtering_istream fistr_read;// io::filtering_istream();
    if (filename_read.rfind(".gz") == filename_index.size() - 3) {
      fistr_read.push(io::gzip_decompressor());
    }
    fistr_read.push(io::file_source(filename_read));
    fistr_read.set_auto_close(true);

    size_t num_lines = 0;
    bool resolved = false;
    //map<ullong,vector<pair<string,string>> seqs;
    //map<ullong>,
    string seq;
    ullong code;
    map<ullong,ullong> counts;
    string umi;
    ullong wellcode = 0;
    ullong cellcode = 0;
    //pair<ullong,ullong> well_code_wm; // with mask
    //pair<ullong,ullong> cell_code_wm; // with mask
    //string qual;
    while (true) {
      string li, lr;
      if (!getline(fistr_index, li) || !getline(fistr_read, lr)) {
	break;
      }
      //cout << num_lines << " " << li << endl;
      size_t slot = num_lines & 3;
      if (slot == 0) {
	wellcode = cellcode = 0;
	size_t pos = li.rfind(':');
	if (pos == string::npos || li.size() - pos < well_code_length) {
	  resolved = false;
	} else {
	  pair<ullong,ullong> well_code_wm = encode_sequence(li.substr(pos + 1, well_code_length));
	  int num_N = popcount(well_code_wm.second) >> 1;
	  //cout << li.substr(pos + 1, well_code_length) << " : " << well_code_wm.second << " => " << num_N << endl;
	  if (num_N > mismatches) {
	    resolved = false;
	  } else {
	    ullong wc = well_code_wm.first;
	    ullong mask = 0xffffffffffffffff ^ well_code_wm.second;
	    for (auto it = well_codes.begin(); it != well_codes.end(); it++) {
	      ullong c = *it;
	      ullong a = (c ^ wc) & mask;
	      ullong mm = (a & 0x5555555555555555) | ((a & 0xaaaaaaaaaaaaaaaa) >> 1);
	      int num_mm = popcount(mm) + num_N;
	      //cout << li.substr(pos + 1, well_code_length) << " " << hex << *it << " " << num_mm << endl;
	      if (num_mm <= mismatches) {
		resolved = true;
		wellcode = c;//code = c << (cell_code_length << 1);
		break;
	      }
	    }
	  }
	}
      } else if (slot == 3) {
	if (cellcode != 0) {
	  //if (resolved) {
	  ostr << lr << "\n";
	}
	// cell code
      } else if (slot == 1) {
	//cout << resolved << " " << li.size() << " " << cell_code_length << endl;
	if (resolved && li.size() >= cell_code_length) {
	  pair<ullong,ullong> cell_code_wm = encode_sequence(li.substr(0, cell_code_length));
	  int num_N = popcount(cell_code_wm.second) >> 1;
	  //cout << num_N << " " << mismatches << endl;
	  if (num_N <= mismatches) {
	    ullong cc = cell_code_wm.first;
	    ullong mask = 0xffffffffffffffff ^ cell_code_wm.second;
	    umi = li.substr(cell_code_length, li.size() - cell_code_length);
	    resolved = false;
	    //cout << hex << cc << " " << mask << " " << cell_codes.size() << endl;
	    for (auto it = cell_codes.begin(); it != cell_codes.end(); it++) {
	      ullong c = *it;
	      ullong a = (c ^ cc) & mask;
	      ullong mm = (a & 0x5555555555555555) | ((a & 0xaaaaaaaaaaaaaaaa) >> 1);
	      int num_mm = popcount(mm) + num_N;
	      if (num_mm <= mismatches) {
		//cout << "OK " << num_mm << " " << decode_index(cell_code_length, c) << " " << decode_index(cell_code_length, cc) << ":" << li.substr(0, 16) << " " << umi << endl;
		//		cout << "OK\n";
		//code |= c;
		//cout << hex << code << endl;
		cellcode = c;
		resolved = true;
		break;
	      }
	    }
	  }
	  if (cellcode != 0) {//resolved) {
	    ullong code = (ullong)((ullong)wellcode << (ullong)(cell_code_length << (ullong)1)) | (ullong)cellcode;
	    auto it = counts.find(code);
	    int n;
	    if (it == counts.end()) {
	      n = 1;
	      counts.insert(make_pair(code, 1));
	    } else {
	      n = (it->second++);
	    }
	    //cout << hex << code << "=>" << decode_index(well_code_length + cell_code_length, code) << ", UMI:" << umi << dec << endl;
	    //cout << "@" << decode_index(well_code_length, wellcode) << "_" << decode_index(cell_code_length, cellcode) << ":" << n << ":" << umi << endl;
	    //cout << "@" << decode_index(8, code >> 32) << "_" << decode_index(16, cell_code_length & 0xffffffff) << ":" << n << ":" << umi << endl;
	    //cout << "@" << decode_index(cell_code_length + well_code_length, code, well_code_length) << ":" << n << ":" << umi << "\n" << lr << "\n+\n";
	    //ostr << "@" << decode_index(cell_code_length + well_code_length, code) << ":" << n << ":" << umi << "\n" << lr << "\n+\n";
	    ostr << "@" << decode_index(well_code_length, wellcode) << "_" << decode_index(cell_code_length, cellcode) << ":" << n << ":" << umi << "\n" << lr << "\n+\n";
	  } else {
	    code = 0;
	  }
	}
      }
      num_lines++;
    }

    io::close(ostr);
    //ostr.sync();
    fistr_index.sync();
    fistr_read.sync();
  }
}


int main(int argc, char** argv) {
  namespace io = boost::iostreams;
  namespace po = boost::program_options;
  po::options_description opt("options");
  //string filename_sam;
  string filename_index;
  string filename_out;
  string filename_read;
  map<string,int> wellcode;
  //map<string,int> wellcode_strict;
  ullong range_max = 0;
  ullong range_min = 0;
  int lap;
  //int shift = 60;
  int well_barcode_size = 8;
  int cell_barcode_size = 16;
  bool verbose = false;
  float threshold_well_composition = 0.05f;
  int minimum_reads_per_cell = 500;
  int mismatches = 1;
  opt.add_options()
    ("output,o", po::value<string>(&filename_out)->default_value(""), "output file")
    ("R1,1", po::value<string>(&filename_index), "fastq.gz")
    ("R2,2", po::value<string>(&filename_read), "fastq.gz")
    ("range_max", po::value<ullong>(&range_max)->default_value(0), "Index range")
    ("range_min", po::value<ullong>(&range_min)->default_value(0), "Index range")
    ("barcode_size,b", po::value<int>(&cell_barcode_size)->default_value(16), "Barcode length (16)")
    ("lap,L", po::value<int>(&lap)->default_value(0), "log period")
    ("minimum_reads,m", po::value<int>(&minimum_reads_per_cell)->default_value(500), "Minimum reads per cell (default 500)")
    ("mismatches,k", po::value<int>(&mismatches)->default_value(1), "Mismatches (default 1)")
    ("verbose,V", "verbose mode")
  ("help,h", "display help");
  //("samfile,s", po::value<string>(&filename_sam), "sam file")

  ullong well_index_capacity = 100;
  ullong cell_index_capacity = 0x1fff;

  // if (true) {
  //   for (int i = 0; i < 100; i++) {
  //     cout << i << " " << popcount(i) << endl;
  //   }
  //   exit(0);
  // }
  
  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, opt), vm);
  } catch (exception& e) {
    cerr << opt << endl;
    cerr << e.what() << endl;
    return -1;
  }
  po::notify(vm);
  verbose = vm.count("verbose");
  if (verbose && lap < 1000) {
    lap = 1000;
  }
  lap *= 1000;
  // {
  //   map<ullong,ullong> test;
  //   vector<string> a = {"CAAAAAAA", "AAAAAAAT", "GGGGGGGG", "NNTTTTTT", "NNNNNNNN"}; 
  //     //char*[] a = {
  //   for (auto b : a) {
  //     auto test = encode_sequence(b);
  //     cout << hex << b << " " << test.first << " " << test.second << endl;
  //   }
  //   exit(0);
  // }
  
  //cout << opt << endl;
  //cout << vm << endl;
  if (vm.count("help")) {
    cout << opt << endl;
  } else {
    try {
      pair<map<string,ullong>,map<string,ullong>> index_sets = load_index_sets(filename_index, well_barcode_size, cell_barcode_size, well_index_capacity, cell_index_capacity, lap, verbose);
      size_t num_total = 0;
      for (auto it = index_sets.first.begin(); it != index_sets.first.end(); it++) {
	num_total += it->second;
      }
      
      
      vector<string> well_indexes = obtain_reduced_indexes(index_sets.first, int(num_total * threshold_well_composition));
      vector<string> cell_indexes = obtain_reduced_indexes(index_sets.second, minimum_reads_per_cell);
      //cout << well_indexes.size() << " " << cell_indexes.size() << endl;
      for (int i = 0; i < 1; i++) {
	cout << well_indexes[i] << " " << cell_indexes[i] << endl;
      }

      //cout << "WELL " << well_indexes.size() << "\n";
      vector<ullong> well_codes = encode_indexes(well_indexes, well_barcode_size);
      //cout << "CELL " << cell_indexes.size() << "\n";
      vector<ullong> cell_codes = encode_indexes(cell_indexes, cell_barcode_size);

      // cout << well_codes.size() << endl;
      // cout << cell_codes.size() << endl;
      // cout << minimum_reads_per_cell << endl;
      // exit(0);

      io::filtering_ostream fostr;
      fostr.set_auto_close(true);
      ostream* ost_ = &cout;
      if (filename_out != "") {
	ost_ = new ofstream(filename_out.c_str());
	if (filename_out.rfind(".gz") == filename_out.size() - 3) {
	  fostr.push(io::gzip_compressor());
	}
      }
      fostr.push(*ost_);//filename_out);
      generate_renamed_fastq(filename_index, well_barcode_size, well_codes, cell_barcode_size, cell_codes,
			     filename_read, 
			     mismatches, fostr);

      if (filename_out != "") {
	reinterpret_cast<ofstream*>(ost_)->close();
	delete ost_;
      }
	

      // if (filename_out.rfind(".gz") == filename_out.size() - 3) {
      // 	fostr.push(io::gzip_compressor());
      // }

      
//ostream* ost = &cout;
//      ostream& ost = cout;
//       io::filtering_ostream fostr;
// //      ofstream* fstr = NULL;
//       if (filename_out.size() > 0) {
// 	fostr.push(filename_out);
// 	//fstr = ofstream(filename_out);
// 	if (filename_out.rfind(".gz") == filename_out.size() - 3) {
// 	  fostr.push(io::gzip_compressor());
// 	}
// 	//fstr = ofstream(filename_out);
// 	ost = &fostr;
//       }
//      fostr.push(ost);
      //generate_renamed_fastq(filename_index, filename_read, well_codes, cell_codes, mismatches, fostr);
      
      // vector<string> cell_index = obtain_reduced_indexes(index_sets.second, minimum_reads_per_cell);
      // vector<ullong> well_codes = encode_indexes(well_index, well_barcode_size);
      // vector<ullong> cell_codes = encode_indexes(cell_index, barcode_size);
      
      //process_fastq_files(filename_index, filename_read, filename_out
      
      //pair<vector<ullong>,vector<ullong>> available_index_codes = obtain_reduced_index_sets(filename_index, barcode_size, well_index_capacity, cell_index_capacity, verbose);

      //vector<ullong> available_well_index;
    } catch (exception& e) {
      cerr << e.what() << endl;
    }
  }
  return 0;
}
