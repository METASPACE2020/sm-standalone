convert: convert.cpp imzml_reader.cpp imzml_reader.h imzb_writer.h imzb_reader.h
	g++ -std=c++11 -O2 convert.cpp imzml_reader.cpp -lxml2 -lblosc -o convert

cffi:
	g++ -shared -fPIC -std=c++11 -O2 imzb_cffi.cpp -pthread -lblosc -o libimzb.so
