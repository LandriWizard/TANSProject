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
#include "MyVertex.h"

#define FALSE 0
#define TRUE 1
#define DEBUG FALSE

using namespace std;



void Reconstruction(){

//Stopwatch declaration and start
  TStopwatch Clock;
  Clock.Start();


  double multiscattering_angle = 0.0012; //in mrad
  double delta_phi = TMath::ASin(3./7.*TMath::Sin(multiscattering_angle)); //used in the while loop to "slice" the azimuth angle in many parts for the reconstruction
  MyVertex* Vertex = new MyVertex();
  
  // Dichiarazione TClonesArray
//  TClonesArray *hitsL1 = new TClonesArray("MyPoint",100);
//  TClonesArray *hitsL2 = new TClonesArray("MyPoint",100);
  TClonesArray *HitsL1 = new TClonesArray("MySignal",100);
  TClonesArray *HitsL2 = new TClonesArray("MySignal",100);
  //Apertura file di input
  TFile hfile("simulation.root");
  //Lettura TTree  e branch
  TTree *tree = (TTree*)hfile.Get("T");
  TBranch *b1 = tree->GetBranch("VertMult");
  TBranch *b2 = tree->GetBranch("HitsL1");
  TBranch *b3 = tree->GetBranch("HitsL2");

  // Definizione degli indirizzi per la lettura dei dati su ttree
//  b1->SetAddress(&point.X);
//  b2->SetAddress(*s);
  b1->SetAddress(&Vertex);
  b2->SetAddress(&HitsL1);
  b3->SetAddress(&HitsL2);

  // loop sugli ingressi nel TTree
  for(int ev = 0; ev < tree->GetEntries(); ev++){
    tree->GetEvent(ev);

    #if DEBUG == TRUE
      cout << "Evento " << ev << "; Molteplicita= " << Vertex->GetMult() << endl;
      cout << "X,Y,Z = " << Vertex->GetX() << "; " << Vertex->GetY() << "; " << Vertex->GetZ() << endl;
      int num = hitsL1->GetEntries();
      cout << "Numero di elementi nel TClonesArray " << num << endl;
      for (int j=0; j < num; j++){
//        MyPoint *tst1=(MyPoint*)hitsL1->At(j);
//        MyPoint *tst2=(MyPoint*)hitsL2->At(j);
//        cout << "Hit on L2 # " << j << ") x, y, z = " << tst2->GetX() << "; " << tst2->GetY() << "; " << tst2->GetZ() << endl;
        MySignal *tst1=(MySignal*)HitsL1->At(j);
        MySignal *tst2=(MySignal*)HitsL2->At(j);
        cout << "Hit on L2 # " << j << ") z, phi = " << tst2->GetZ() << "; " << tst2->GetPhi() << "; " << endl;
      }
    #endif

    

  }
 
//Clock stop and time print
  Clock.Stop();
  Clock.Print();

}