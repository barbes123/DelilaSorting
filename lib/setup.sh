#!/bin/bash

rootcint -f DelilaEventDict.cxx -c DelilaEvent.h DelilaEventLinkDef.h
rootcint -f LaBrTreeEventDict.cxx -c LaBrTreeEvent.h LaBrTreeEventLinkDef.h
rootcint -f ElissaTreeEventDict.cxx -c ElissaTreeEvent.h ElissaTreeEventLinkDef.h
#rootcint -f LaBrEventDict.cxx -c LaBrEvent.h LaBrEventLinkDef.h
#rootcint -f ElissaEventDict.cxx -c ElissaEvent.h ElissaEventLinkDef.h
g++ -shared -o libDelilaEvent.so`root-config --ldflags`  -std=c++11 -fPIC -O3 -I`root-config --incdir` DelilaEventDict.cxx
g++ -shared -o libLaBrTreeEvent.so`root-config --ldflags`  -std=c++11 -fPIC -O3 -I`root-config --incdir` LaBrTreeEventDict.cxx
g++ -shared -o libElissaTreeEvent.so`root-config --ldflags`  -std=c++11 -fPIC -O3 -I`root-config --incdir` ElissaTreeEventDict.cxx
#g++ -shared -o libLaBrEvent.so`root-config --ldflags`  -std=c++11 -fPIC -O3 -I`root-config --incdir` LaBrEventDict.cxx
#g++ -shared -o libElissaEvent.so`root-config --ldflags`  -std=c++11 -fPIC -O3 -I`root-config --incdir` ElissaEventDict.cxx


