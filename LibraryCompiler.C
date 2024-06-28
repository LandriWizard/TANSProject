#include <Riostream.h>
#include <TROOT.h>
#include <TSystem.h> 

void LibraryCompiler() {

  //Librerie utili per la simulazione
  gROOT->LoadMacro("MyPoint.cxx");
  gROOT->LoadMacro("MyVertex.cxx");




}