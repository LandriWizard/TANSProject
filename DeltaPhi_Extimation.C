#include "Riostream.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TMath.h"
#include "TStopwatch.h"
#include "TTree.h"

#include "MyParticle.h"
#include "MyPhysics.h"
#include "MyPoint.h"
#include "MyRandom.h"
#include "MySignal.h"
#include "MyTracklet.h"
#include "MyVertex.h"

#define FALSE 0
#define TRUE 1
#define DEBUG FALSE

using namespace std;

void DeltaPhi_Extimation(const char* input_file = "simulation.root"){
  //Stopwatch declaration and start
  TStopwatch Clock;
  Clock.Start();

  //physical quantities
  double multiscattering_angle = 0.0012;
  double inner_radius = 4.;
  double outer_radius = 7.;

  double reconstructed_z = 0;;
  double residual_z;

  //Declaring "auxiliary" objects
  MyVertex* Vertex = new MyVertex();
  MyTracklet* Tracklet = new MyTracklet(inner_radius,outer_radius,0.,0.);
  
  //Declaring TClonesArray
  TClonesArray *HitsL1 = new TClonesArray("MySignal",70);
  TClonesArray *HitsL2 = new TClonesArray("MySignal",70); 
 
  //Opening input file
  TFile infile(input_file);
  if(infile.IsZombie()){
    cout << "There is no file named " << input_file << " in the chosen directory" << endl;
    return;
  }

  //Creating the histogram for delta phi
  TH1D* hphi = new TH1D("hphi", "Istogramma delle phi", 100, -0.01, 0.01);

  //Reading the tree
  TTree *tree = (TTree*)infile.Get("T");
  TBranch *b1 = tree->GetBranch("VertMult");
  TBranch *b2 = tree->GetBranch("HitsL1");
  TBranch *b3 = tree->GetBranch("HitsL2");

  //Address definition to read the data
  b1->SetAddress(&Vertex);
  b2->SetAddress(&HitsL1);
  b3->SetAddress(&HitsL2);

  for(int ev = 0; ev < tree->GetEntries(); ev++){
    tree->GetEvent(ev);

    //Print event number every 100000 events to check the progress
    if(ev%100000 == 0) cout << "Event #" << ev << endl;

    //Loop over the hits on the second layer
    for(int i = 0; i < HitsL2->GetEntries(); i++){
      MySignal* InteractionOnLayer2 = (MySignal*)HitsL2->At(i);
      int layer2_flag = InteractionOnLayer2->GetFlag();
      //Loop over the hits on the first layer
      for(int j = 0; j < HitsL1->GetEntries(); j++){
        MySignal* InteractionOnLayer1 = (MySignal*)HitsL1->At(j);
        //Check if the hits are on the same track, if so, compute delta phi and fill the histogram
        if(InteractionOnLayer1->GetFlag() == layer2_flag){
          hphi->Fill((InteractionOnLayer1->GetPhi())-(InteractionOnLayer2->GetPhi()));
        }
      }
    }
    hphi->DrawCopy();
  }

  //Print the time elapsed and end the function
  Clock.Stop();
  Clock.Print();
  return;
}