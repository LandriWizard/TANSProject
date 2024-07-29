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

#define TRUE 1
#define FALSE 0
#define DEBUG FALSE

using namespace std;

//IMPORTANT: DISTANCES ARE MEASURED IN CENTIMETRES IN THIS SIMULATION
//MULTIPLICITY FLAG VALUES: 1 FOR EXTRACTION FROM GIVEN DISTIBUTION, 2 FOR CONSTANT VALUE, 3 FOR UNIFORM DISTRIBUTION

void Simulation(int N_exp = 1e6, unsigned int seed = 69420, int multiplicity_flag = 1, const char* input_file = "kinem.root", const char* output_file = "simulation.root"){

  MyRandom *RndmPtr = new MyRandom(input_file,seed);
  delete gRandom;
  gRandom = RndmPtr;
  if(RndmPtr->GetFlag()) cout << "There is no file named " << input_file << " in the chosen directory" << endl;
  else cout << "File " << input_file << " was found, beginning the operations" << endl;

//Stopwatch declaration and start
  TStopwatch Clock;
  Clock.Start();

//Output file declaration
  TFile hfile(output_file,"RECREATE");

  double X,Y,Z;
  int mult;

  TTree* Tree = new TTree("T","Tree with 3 branches");
//  TClonesArray* HitsOnL1 = new TClonesArray("MyPoint",100);
//  TClonesArray& L1Hit = *HitsOnL1;
//  TClonesArray* HitsOnL2 = new TClonesArray("MyPoint",100);
//  TClonesArray& L2Hit = *HitsOnL2;
  TClonesArray* HitsOnL1 = new TClonesArray("MySignal",100);
  TClonesArray& L1Hit = *HitsOnL1;
  TClonesArray* HitsOnL2 = new TClonesArray("MySignal",100);
  TClonesArray& L2Hit = *HitsOnL2;

  MyPoint* Point = new MyPoint();
  MyVertex* Vertex = new MyVertex();


  Tree->Branch("VertMult",&Vertex);
  Tree->Branch("HitsL1",&HitsOnL1);
  Tree->Branch("HitsL2",&HitsOnL2);

  Tree->SetAutoSave(0);

  //Generatori di vertice, inizio con la generazione di un vertice con molteplicità estratta da kinem.root
  //Ancora da inserire: molteplicità fissa e estratta da distribuzione uniforme

  MyPoint* Hit = new MyPoint(); //Points used to store the true position of where the particles hit the detectors
  MyParticle* Particle = new MyParticle(); //The particle that will be transported, contains the theta and phi of the trajectory

  MyPhysics BeamPipe(3.,54.);
  MyPhysics Layer1(4.,27.);
  MyPhysics Layer2(7.,27.);

  //Events loop
  for(int i = 0; i < N_exp; i++){

    if(i%100000 == 0) cout << "Vertex #" << i << endl;

    //Vertex generation
    Vertex->SetX(RndmPtr->Gaus(0.,0.01));
    Vertex->SetY(RndmPtr->Gaus(0.,0.01));
    Vertex->SetZ(RndmPtr->Gaus(0.,5.3));
    Vertex->SetMult(RndmPtr->RndmMult());
    mult = Vertex->GetMult();

    #if DEBUG == TRUE
      cout << "Generated vertex #" << i << " = (" << Vertex->GetX() << ", " <<
                                                     Vertex->GetY() << ", " <<
                                                     Vertex->GetZ() << ");" << endl;
      cout << "Random multiplicity = " << mult << ";" << endl;
    #endif

    //Multiplicity loop
    int j1 = 0, j2 = 0; //Indices used to save hits when the z is actually on the detector
    for(int j = 0; j < mult; j++){

      #if DEBUG == TRUE
        cout << "Indeces at the beginning of the for cycle:" << endl;
        printf("j = %d\n",j);
        printf("j1 = %d\n",j1);
        printf("j2 = %d\n",j2);
      #endif


      //Particle generation
      Particle->SetTheta(RndmPtr->RndmTheta());
      Particle->SetPhi(RndmPtr->Uniform(0.,2.*TMath::Pi()));

      //Particle transport
      //Beam pipe interaction
      MyPoint* tmpPoint = new MyPoint(Vertex->GetX(),Vertex->GetY(),Vertex->GetZ());
      *Hit = BeamPipe.Hit(tmpPoint, Particle);  //Particle transport
      #if DEBUG == TRUE
        cout << "Hit position on the beam pipe = (" << Hit->GetX() << ", " <<
                                                       Hit->GetY() << ", " <<
                                                       Hit->GetZ() << "); Radius of the position = " << 
                                                       Hit->GetRadius() << endl;
      #endif


      //First layer interaction
      *Hit = Layer1.Hit(Hit, Particle); //Particle transport

      #if DEBUG == TRUE
        cout << "Hit position on the first detector layer = (" << Hit->GetX() << ", " <<
                                                                  Hit->GetY() << ", " <<
                                                                  Hit->GetZ() << "); Radius of the position = " << 
                                                                  Hit->GetRadius() << endl;
      #endif

      if(Hit->GetZ() > -1.*Layer1.GetH()/2. && Hit->GetZ() < Layer1.GetH()/2.){ //Check if Z is on the detector

        new(L1Hit[j1])MySignal(Hit);

        //Second layer interaction
        *Hit = Layer2.Hit(Hit, Particle);
        #if DEBUG == TRUE
          cout << "Hit position on the second detector layer = (" << Hit->GetX() << ", " <<
                                                                    Hit->GetY() << ", " <<
                                                                    Hit->GetZ() << "); Radius of the position = " << 
                                                                    Hit->GetRadius() << endl;
        #endif

        if(Hit->GetZ() > -1.*Layer2.GetH()/2. && Hit->GetZ() < Layer2.GetH()/2.){ //Check if Z is on the detector

          new(L2Hit[j2])MySignal(Hit);



//        #if DEBUG == TRUE
//          printf("Evento %d - moltepl: %d - interazione: %d\n",i,mult,j+1);
//          printf("x= %f ; y= %f; z= %f \n",Vertex->GetX(),Vertex->GetY(),Vertex->GetZ());
//          printf("Entries nel TClonesArray1: %d\n",HitsOnL1->GetEntries());
//          MySignal *tst1=(MySignal*)HitsOnL1->At(j1);
//          std::cout<< "Hit on L1 " << j1 << ") phi, z = " << tst1->GetPhi() << "; " << tst1->GetZ() << std::endl;
//          printf("Entries nel TClonesArray2: %d\n",HitsOnL2->GetEntries());
//          MySignal *tst2=(MySignal*)HitsOnL2->At(j2);
//          std::cout << "Hit on L2 " << j2 << ") phi, z = " << tst2->GetPhi() << "; " << tst2->GetZ() << std::endl;
//        #endif



          j2++;

        }

        j1++;

      }



      #if DEBUG == TRUE
        cout << "Indeces at the end of the for cycle:" << endl;
        printf("j = %d\n",j);
        printf("j1 = %d\n",j1);
        printf("j2 = %d\n\n\n",j2);
      #endif

    }


    Tree->Fill();

    HitsOnL1->Clear();
    HitsOnL2->Clear();

  }

//Writing and saving the output file
  hfile.Write();
  hfile.Close();

//Stopwatch stop and time print
  Clock.Stop();
  Clock.Print();

//Deallocating pointers
  delete Hit;
  delete Particle;

//NECESSARIO AGGIUNGERE SMEARING E MULTISCATTERING
//AGGIUNGERE FUNTORI PER MOTLEPLICITA' FISSA E DA DISTRIBUZIONE UNIFORME
//OPPORTUNO INIZIARE ANCHE LA RICOSTRUZIONE

//SCATTERING MULTIPLO: 1.2 MILLIRADIANTI PER SOVRASTIMARE IL FENOMENO


}