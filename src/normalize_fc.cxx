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
#include <string>

#include <boost/program_options.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string.hpp>
//#include <boost/regex.hpp>

typedef unsigned long long ullong;

using namespace std;

int main(int argc, char** argv) {
  namespace io = boost::iostreams;
  namespace po = boost::program_options;

  string filename_out;
  vector<string> filenames_fc;
  
  po::options_description opt("options");
  opt.add_options()
    ("output,o", po::value<string>(&filename_out)->default_value("out"), "output file")
    ("input,i", po::value<vector<string>>(&filenames_fc)->multitoken(), "featureCount output")
    ("verbose,V", "verbose mode")
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
  bool verbose = vm.count("verbose");

  if (filenames_fc.size() == 0) {
    cerr << "no input files given" << endl;
    return -1;
  }

  vector<string> sample_names;
  vector<vector<int>> data;
  vector<vector<string>> info;
  bool first = true;
  for (auto filename : filenames_fc) {
    size_t rpos = filename.rfind("/");
    string stem;
    if (rpos != string::npos) {
      stem = filename.substr(rpos + 1);
    } else {
      stem = filename;
    }
    rpos = stem.rfind(".");
    if (rpos != string::npos) {
      stem = stem.substr(0, rpos);
    }
      
    ifstream fi(filename.c_str());
    if (!fi.is_open()) {
      throw runtime_error(string("cannot open ") + filename);
    }
    //vector<string> samples;
    int num_columns = 0;
    int row_index = 0;
    if (verbose) {
      cerr << "loading " << filename << endl;
    }
    //vector<string> sample_ids;
    while (!fi.eof()) {
      string line;
      list<string> items_;
      getline(fi, line);
      if (line.c_str()[0] == '#') continue;
      if (verbose && (row_index + 1) % 100 == 0) {
	cerr << row_index << " : " << line.substr(0, 60) << endl;
      }
      boost::split(items_, line, boost::is_any_of("\t"));
      if (items_.size() <= 1) break;
      //cout << "SIZE " << items_.size() << endl;
      
      //cout << num_columns << " " << items_.size() << endl;

      if (num_columns == 0) {
	int index = 0;
	for (auto item : items_) {
	  if (index++ < 6) {
	    continue;
	  }
	  string label = stem + string(":") + std::to_string(num_columns + 1);
	  string name = label;
	  int count = 0;
	  while (true) {
	    if (find(sample_names.begin(), sample_names.end(), name)
		!= sample_names.end()) {
	      count ++;
	      name = label + string(".") + std::to_string(count);
	    } else {
	      sample_names.push_back(name);
	      break;
	    }
	  }
	  // if (verbose) {
	  //   cout << " " << name << endl;
	  // }
	  num_columns++;
	}
	if (verbose) cerr << stem << ": " << num_columns << " columns" << endl;
      } else {
	vector<string> items{std::begin(items_), std::end(items_)};
	//cout << items.size() << " columns" << endl;
	if (first) {
	  info.push_back(std::vector<string>(items.begin(), items.begin() + 6));
	  vector<int> row;
	  for (auto it = items.begin() + 6; it != items.end(); it++) {
	    row.push_back(std::stoi(*it));
	  }
	  //cout << row.size() << " / " << num_columns << endl;
	  data.push_back(row);
	  // if (data.size() % 1000 == 0) {
	  //   cerr << row[0] << endl;
	  // }
	  if (row.size() != num_columns) {
	    throw runtime_error("inconsistent column size");
	  }
	} else {
	  //cout << "#########row index " << row_index << " " << items[0] << endl; 
	  if (items[0] != info[row_index][0]) {
	    // for (auto ri : info) {
	    //   cout << ri[0] << endl;
	    // }
	    throw runtime_error(string("invalid index row ") + *(items.begin()) + string(" <= ") + info[row_index][0]);
	  }
	  vector<int>& row = data[row_index];
	  for (auto it = items.begin() + 6; it != items.end(); it++) {
	    row.push_back(std::stoi(*it));
	  }
	  if (items.size() - 6 != num_columns) {
	    throw runtime_error("inconsistent column size");
	  }
	}
	// for (auto row : data) {
	//   cout << ", " << row.size();
	// }
	// cout << endl;
	row_index++;
      }
    }
    fi.close();
    first = false;
  }
  if (verbose) {
    cerr << "integration\n";
  }
  pair<int,int> shape(data.size(), data[0].size());
  vector<int> total;
  for (int i = 0; i < shape.second; i++) total.push_back(0);
  for (auto row : data) {
    for (int i = 0; i < (int)row.size(); i++) {
      total[i] += row[i];
    }
  }
  vector<double> coeff;
  for (auto t: total) {
    if (t == 0) {
      coeff.push_back(0.0);
    } else {
      coeff.push_back(1e6 / t);
    }
  }
  
  ofstream file_count(string(filename_out + ".cnt").c_str());
  ofstream file_tpm(string(filename_out + ".tpm").c_str());
  // header
  file_count << "Geneid";
  file_tpm << "Geneid";
  for (auto item : sample_names) {
    file_count << "\t" << item;
    file_tpm << "\t" << item;
  }
  file_count << "\n";
  file_tpm << "\n";

  // for (int i = 0; i < shape.first; i++) {
  //   file_count << info[i][0];
  //   file_tpm << info[i][0];
  // }
  int row_index = 0;
  for (auto row : data) {
    file_count << info[row_index][0];
    file_tpm << info[row_index][0];
    for (int i = 0; i < shape.second; i++) {
      file_count << "\t" << row[i];
      file_tpm << "\t" << row[i] * coeff[i];
    }
    file_count << "\n";
    file_tpm << "\n";
    row_index++;
  }
  file_count.close();
  file_tpm.close();
}
