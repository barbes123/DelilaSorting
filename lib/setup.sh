#!/bin/bash

rootcint -f DelilaEventDict.cxx -c DelilaEvent.h DelilaEventLinkDef.h
rootcint -f HPGeTreeEventDict.cxx -c HPGeTreeEvent.h HPGeTreeEventLinkDef.h
rootcint -f LaBrTreeEventDict.cxx -c LaBrTreeEvent.h LaBrTreeEventLinkDef.h
rootcint -f ElissaTreeEventDict.cxx -c ElissaTreeEvent.h ElissaTreeEventLinkDef.h
g++ -shared -o libDelilaEvent.so`root-config --ldflags`  -std=c++11 -fPIC -O3 -I`root-config --incdir` DelilaEventDict.cxx
g++ -shared -o libHPGeTreeEvent.so`root-config --ldflags`  -std=c++11 -fPIC -O3 -I`root-config --incdir` HPGeTreeEventDict.cxx
g++ -shared -o libLaBrTreeEvent.so`root-config --ldflags`  -std=c++11 -fPIC -O3 -I`root-config --incdir` LaBrTreeEventDict.cxx
g++ -shared -o libElissaTreeEvent.so`root-config --ldflags`  -std=c++11 -fPIC -O3 -I`root-config --incdir` ElissaTreeEventDict.cxx


