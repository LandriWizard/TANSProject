#include <Riostream.h>
#include <TROOT.h>
#include <TSystem.h> 

void LibraryCompiler(TString compiling_flag="fast") {

  TString Compiling_option;
  if(compiling_flag.Contains("force")) Compiling_option = "kfg";
  else Compiling_option = "kg";

  //Librerie utili per la simulazione

  gSystem->CompileMacro("MyPoint.cxx",Compiling_option.Data());
  gSystem->CompileMacro("MyVertex.cxx",Compiling_option.Data());
  gSystem->CompileMacro("MyRandom.cxx",Compiling_option.Data());
  gSystem->CompileMacro("MyParticle.cxx",Compiling_option.Data());
  gSystem->CompileMacro("MySignal.cxx",Compiling_option.Data());
  gSystem->CompileMacro("MyPhysics.cxx",Compiling_option.Data());
  gSystem->CompileMacro("MyTracklet.cxx",Compiling_option.Data());

  gSystem->CompileMacro("DeltaPhi_Extimation.C",Compiling_option.Data());

  gSystem->CompileMacro("LampadinaSimulation.C",Compiling_option.Data());
  gSystem->CompileMacro("LampadinaReconstruction.C",Compiling_option.Data());



}