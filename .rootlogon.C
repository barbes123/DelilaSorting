{
  TString home_path = gSystem->GetFromPipe("echo $HOME");
  gROOT->ProcessLine(Form(".L %s/DelilaSorting/lib/libDelilaEvent.so-m64",home_path.Data()));
  gROOT->ProcessLine(Form(".L %s/DelilaSorting/lib/libHPGeTreeEvent.so-m64",home_path.Data()));
  gROOT->ProcessLine(Form(".L %s/DelilaSorting/lib/libLaBrTreeEvent.so-m64",home_path.Data()));
  gROOT->ProcessLine(Form(".L %s/DelilaSorting/lib/libElissaTreeEvent.so-m64",home_path.Data()));
   
  }
