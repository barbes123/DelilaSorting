#include <iostream>
#include <vector>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"


#include <sstream>
#include <string>
#include <iostream> 


class TreeData
{
public:
  TreeData(){};
  ~TreeData(){};

  UChar_t Mod;
  UChar_t Ch;
  ULong64_t TimeStamp;
  Double_t FineTS;
  UShort_t ChargeLong;
  UShort_t ChargeShort;
  UInt_t Extras;
  UInt_t RecordLength;
  std::vector<uint16_t> Signal;
  // std::vector<uint16_t> Trace2;
  // std::vector<uint8_t> DTrace1;
  // std::vector<uint8_t> DTrace2;

  void Print()
  {
    std::cout << Int_t(Mod) <<"\t"<< Int_t(Ch) <<"\t"<< TimeStamp <<"\t"<< FineTS <<"\t"
	      << ChargeLong <<"\t"<< ChargeShort <<"\t"<< RecordLength << std::endl;
  };
};

void sort(TString fileName = "run75_720_ssgant1.root")
{
  if(fileName == "") return;

  // Load tree 
  auto file = new TFile(fileName);
  auto tree = (TTree*)file->Get("ELIADE_Tree");
  TreeData event;
  event.Signal.resize(1);// If the signal was recorded, set enough number.  But easy to use all memory space
  
  tree->SetBranchAddress("Mod", &event.Mod);
  tree->SetBranchAddress("Ch", &event.Ch);
  tree->SetBranchAddress("TimeStamp", &event.TimeStamp);
  tree->SetBranchAddress("FineTS", &event.FineTS);
  tree->SetBranchAddress("ChargeLong", &event.ChargeLong);
  tree->SetBranchAddress("ChargeShort", &event.ChargeShort);
  tree->SetBranchAddress("RecordLength", &event.RecordLength);
  tree->SetBranchAddress("Signal", &event.Signal[0]);

  const auto nEvents = tree->GetEntries();
  std::vector<TreeData> dataVec;
  dataVec.reserve(nEvents);
  ULong64_t lastHit = 0;
  for(auto iEve = 0; iEve < nEvents; iEve++) {
    tree->GetEntry(iEve);
    event.FineTS = event.TimeStamp * 1000 + event.FineTS;
    dataVec.push_back(event);
  }
  file->Close();
  delete file;

  // Sort by FineTS
  std::sort(dataVec.begin(), dataVec.end(), [](const TreeData &a, const TreeData &b){
      return a.FineTS < b.FineTS;
    });

  // Output tree
  auto outputName = fileName;
  outputName.ReplaceAll(".root", "_sorted.root");
  auto outputFile = new TFile(outputName, "RECREATE");
  auto sorted = new TTree("ELIADE_Tree", "Sorted data");
  
  UChar_t Mod;
  UChar_t Ch;
  ULong64_t TimeStamp;
  Double_t FineTS;
  UShort_t ChargeLong;
  UShort_t ChargeShort;
  UInt_t Extras;
  UInt_t RecordLength;
  UShort_t Signal[100000]{0};
    
  sorted->Branch("Mod", &Mod, "Mod/b");
  sorted->Branch("Ch", &Ch, "Ch/b");
  sorted->Branch("TimeStamp", &TimeStamp, "TimeStamp/l");
  sorted->Branch("FineTS", &FineTS, "Finets/D");
  sorted->Branch("ChargeLong", &ChargeLong, "ChargeLong/s");
  sorted->Branch("ChargeShort", &ChargeShort, "ChargeShort/s");
  sorted->Branch("RecordLength", &RecordLength, "RecordLength/i");
  sorted->Branch("Signal", Signal, "Signal[RecordLength]/s");
  for(auto iEve = 0; iEve < dataVec.size(); iEve++) {
    Mod = dataVec.at(iEve).Mod;
    Ch = dataVec.at(iEve).Ch;
    TimeStamp = dataVec.at(iEve).TimeStamp;
    FineTS = dataVec.at(iEve).FineTS;
    ChargeLong = dataVec.at(iEve).ChargeLong;
    ChargeShort = dataVec.at(iEve).ChargeShort;
    RecordLength = dataVec.at(iEve).RecordLength;
    if(RecordLength > 0)
      std::copy(&dataVec.at(iEve).Signal[0], &dataVec.at(iEve).Signal[RecordLength], Signal);
    
    sorted->Fill();
  }

  sorted->Write();
  outputFile->Close();
  delete outputFile;
}
