#include <Riostream.h>
#include <TROOT.h>
#include <TSystem.h> 

void LibraryCompiler() {

  //Librerie utili per la simulazione
//  gROOT->LoadMacro("MyPoint.cxx");
//  gROOT->LoadMacro("MyVertex.cxx");
//  gROOT->LoadMacro("MyRandom.cxx");
//  gROOT->LoadMacro("MyParticle.cxx");
//  gROOT->LoadMacro("MyPhysics.cxx");

  gSystem->CompileMacro("MyPoint.cxx","kfg");
  gSystem->CompileMacro("MyVertex.cxx","kfg");
  gSystem->CompileMacro("MyRandom.cxx","kfg");
  gSystem->CompileMacro("MyParticle.cxx","kfg");
  gSystem->CompileMacro("MyPhysics.cxx","kfg");

  gSystem->CompileMacro("LampadinaSimulation.C","kfg");
  gSystem->CompileMacro("LampadinaReconstruction.C","kfg");



}