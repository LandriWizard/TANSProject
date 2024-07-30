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
#define DEBUG TRUE

using namespace std;


//IMPORTANT: LENGHTS ARE IN CM, ANGLES IN RAD

void Reconstruction(const char* input_file = "simulation.root"){

//Stopwatch declaration and start
  TStopwatch Clock;
  Clock.Start();

  //physical quantities
  double multiscattering_angle = 0.0012;
  double inner_radius = 4.;
  double outer_radius = 7.;

  double delta_phi = TMath::ASin(3./7.*TMath::Sin(multiscattering_angle)); //I would like to use it in the while loop to "slice" the azimuth angle, cannot because its not a divisor of 2Pi
  int slice_number = 2*TMath::Pi()/delta_phi + 1; //number of azimuth angle slices
  double real_delta_phi = 2*TMath::Pi()/(double)slice_number; //actual value used to divide the azimuth angle, is a divisor of 2Pi; differs in the order ~10^-8 from delta_phi
 
  //Declaring "auxiliary" objects
  MyVertex* Vertex = new MyVertex();
  MyTracklet* Tracklet = new MyTracklet(inner_radius,outer_radius,0.,0.);
  
  // Dichiarazione TClonesArray
  TClonesArray *HitsL1 = new TClonesArray("MySignal",100);
  TClonesArray *HitsL2 = new TClonesArray("MySignal",100);
 
  //Apertura file di input
  TFile infile(input_file);
  if(infile.IsZombie()){
		cout << "There is no file named " << input_file << " in the chosen directory" << endl;
		return;
	}
  //Lettura TTree  e branch
  TTree *tree = (TTree*)infile.Get("T");
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
      int num1 = HitsL1->GetEntries();
      int num2 = HitsL2->GetEntries();
      cout << "Numero di elementi nel primo TClonesArray " << num1 << endl;
      cout << "Numero di elementi nel secondo TClonesArray " << num2 << endl;
      for (int j=0; j < num1; j++){
        MySignal *tst1=(MySignal*)HitsL1->At(j);
        cout << "Hit on L1 # " << j << ") z, phi, r = " << tst1->GetZ() << ";\t " << tst1->GetPhi() << ";\t " << tst1->GetR() << ";\t " << endl;
      }
      for (int j=0; j < num2; j++){
        MySignal *tst2=(MySignal*)HitsL2->At(j);
        cout << "Hit on L2 # " << j << ") z, phi, r = " << tst2->GetZ() << ";\t " << tst2->GetPhi() << ";\t " << tst2->GetR() << ";\t " << endl;
      }
    #endif

    for(int i = 0; i < HitsL1->GetEntries(); i++){ //for cycle on 1st layer
      MySignal* InteractionOnLayer1 = (MySignal*)HitsL1->At(i);
      #if DEBUG == TRUE
        cout << "Interatcion on layer 1 #" << i << "; z, phi, r = " << InteractionOnLayer1->GetZ() << ";\t " << InteractionOnLayer1->GetPhi() 
                                                                                                   << ";\t " << InteractionOnLayer1->GetR() << "; " << endl;
      #endif
      Tracklet->SetZ1(InteractionOnLayer1->GetZ()); //insertion of Z1 in the tracklet

      for(int j = 0; j < HitsL2->GetEntries(); j++){  //for cycle on 2nd layer
      MySignal* InteractionOnLayer2 = (MySignal*)HitsL2->At(j);
      #if DEBUG == TRUE
        cout << "Interatcion on layer 2 #" << j << "; z, phi, r = " << InteractionOnLayer2->GetZ() << ";\t " << InteractionOnLayer2->GetPhi()   
                                                                                                   << ";\t " << InteractionOnLayer2->GetR() << "; " << endl;
      #endif

      }
    }
  }
 
//Clock stop and time print
  Clock.Stop();
  Clock.Print();

//delete pointers
  delete Vertex;

}