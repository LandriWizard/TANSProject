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

void Reconstruction(const char* input_file = "simulation.root", const char* log_file = "reconstruction_log.txt"){

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

  double reconstructed_z = 0;;
  double residual_z;

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

  //Opening log file
  ofstream lfile(log_file);

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

    if(ev%100000 == 0) cout << "Event #" << ev << endl;

    #if DEBUG == TRUE
      lfile << "\nConsidering event #" << ev << endl;
    #endif


//    #if DEBUG == TRUE
//      cout << "Evento " << ev << "; Molteplicita= " << Vertex->GetMult() << endl;
//      cout << "X,Y,Z = " << Vertex->GetX() << "; " << Vertex->GetY() << "; " << Vertex->GetZ() << endl;
//      int num1 = HitsL1->GetEntries();
//      int num2 = HitsL2->GetEntries();
//      cout << "Numero di elementi nel primo TClonesArray " << num1 << endl;
//      cout << "Numero di elementi nel secondo TClonesArray " << num2 << endl;
//      for (int j=0; j < num1; j++){
//        MySignal *tst1=(MySignal*)HitsL1->At(j);
//        cout << "Hit on L1 # " << j << ") z, phi, r = " << tst1->GetZ() << ";\t " << tst1->GetPhi() << ";\t " << tst1->GetR() << ";\t " << endl;
//      }
//      for (int j=0; j < num2; j++){
//        MySignal *tst2=(MySignal*)HitsL2->At(j);
//        cout << "Hit on L2 # " << j << ") z, phi, r = " << tst2->GetZ() << ";\t " << tst2->GetPhi() << ";\t " << tst2->GetR() << ";\t " << endl;
//      }
//    #endif

    int vertex_counter = 0;

    for(int i = 0; i < HitsL1->GetEntries(); i++){ //for cycle on 1st layer



      bool signal_vertex_check = 0; //bool to check if we find a vertex for each signal on layer 1

      MySignal* InteractionOnLayer1 = (MySignal*)HitsL1->At(i);
//      #if DEBUG == TRUE
//        lfile << "Interatcion on layer 1 #" << i << "; z, phi, r = " << InteractionOnLayer1->GetZ() << ";\t " << InteractionOnLayer1->GetPhi() 
//                                                                                                   << ";\t " << InteractionOnLayer1->GetR() << "; " << endl;
//      #endif
      Tracklet->SetZ1(InteractionOnLayer1->GetZ()); //insertion of Z1 in the tracklet

      for(int j = 0; j < HitsL2->GetEntries(); j++){  //for cycle on 2nd layer
        MySignal* InteractionOnLayer2 = (MySignal*)HitsL2->At(j);
//        #if DEBUG == TRUE
//          lfile << "Interatcion on layer 2 #" << j << "; z, phi, r = " << InteractionOnLayer2->GetZ() << ";\t " << InteractionOnLayer2->GetPhi()   
//                                                                                                     << ";\t " << InteractionOnLayer2->GetR() << "; " << endl;
//        #endif

        if(TMath::Abs(InteractionOnLayer2->GetPhi() - InteractionOnLayer1->GetPhi()) < 0.004/*2.*real_delta_phi*/){
          Tracklet->SetZ1(InteractionOnLayer1->GetZ());
          Tracklet->SetZ2(InteractionOnLayer2->GetZ());
          reconstructed_z = Tracklet->Intersection();

          signal_vertex_check = 1;

          #if DEBUG == TRUE
            if(signal_vertex_check){
              lfile << "Reconstructed vertex Z = " << reconstructed_z << " thanks to hit #" << i << " on the 1st detector layer and to hit #" 
                                                                                            << j << " on the 2nd detector layer" << endl;
              if(InteractionOnLayer1->GetFlag() == InteractionOnLayer2->GetFlag())  lfile << "Vertex reconstructed using the same particles" << endl;
              else lfile << "Vertex reconstructed using particle " << InteractionOnLayer1->GetFlag() << " on the 1st layer and particle "
                                                                   << InteractionOnLayer2->GetFlag() << " on the second layer" << endl;
            }
          #endif

          vertex_counter++;
        }


      }

      #if DEBUG == TRUE
        if(!signal_vertex_check) lfile << "!! NOT FOUND ANY VERTEX FOR HIT ON LAYER 1 #" << i << endl;
        if(signal_vertex_check){
          lfile << "Real vertex z cohordinate = " << Vertex->GetZ() << endl;
          lfile << "Residual = " << Vertex->GetZ() - reconstructed_z << endl;
        }
      #endif

    }

    #if DEBUG == TRUE
      lfile << "Number of hits on the 1st layer = " << HitsL1->GetEntries() << endl; 
      lfile << "Maximum number of intersections expected, based on the number of hits on the 2nd layer = " << HitsL2->GetEntries() << endl;
      lfile << "Found " << vertex_counter << " intersections" << endl;
    #endif

  }
 

//closing log file
  lfile.close();

//delete pointers
  delete Vertex;
  delete Tracklet;


//Clock stop and time print
  Clock.Stop();
  Clock.Print();

}