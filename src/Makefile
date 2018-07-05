all:normalize_fc

normalize_fc:normalize_fc.o
	c++ normalize_fc.o -lboost_program_options -lboost_iostreams -lboost_regex -std=c++11 -o normalize_fc

normalize_fc.o:normalize_fc.cxx
	c++ normalize_fc.cxx -std=c++11 -c -g
