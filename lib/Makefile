all: 
	rootcint -f DelilaEventDict.cxx -c DelilaEvent.h DelilaEventLinkDef.h
	rootcint -f LaBrEventDict.cxx -c LaBrEvent.h LaBrEventLinkDef.h
	rootcint -f ElissaEventDict.cxx -c ElissaEvent.h ElissaEventLinkDef.h
	g++ -shared -o libDelilaEvent.so`root-config --ldflags`  -std=c++11 -fPIC -O3 -I`root-config --incdir` DelilaEventDict.cxx
	g++ -shared -o libLaBrEvent.so`root-config --ldflags`  -std=c++11 -fPIC -O3 -I`root-config --incdir` LaBrEventDict.cxx
	g++ -shared -o libElissaEvent.so`root-config --ldflags`  -std=c++11 -fPIC -O3 -I`root-config --incdir` ElissaEventDict.cxx
	g++ -ggdb -fPIC DelilaEventDict.cxx LaBrEventDict.cxx ElissaEventDict.cxx `root-config --cflags` `root-config --libs` -o fake
clean:
	rm *~ *.so

