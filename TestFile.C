#include "Riostream.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TMath.h"

#include "TLeaf.h"

#include "TTree.h"

#include "MyParticle.h"
#include "MyPhysics.h"
#include "MyPoint.h"
#include "MyRandom.h"
#include "MyVertex.h"

#define TRUE 1
#define FALSE 0
#define DEBUG FALSE

using namespace std;

//IMPORTANT: DISTANCES ARE MEASURED IN CENTIMETRES IN THIS SIMULATION

void Test(int N_exp = 1e6, unsigned int seed = 69420, const char* input_file = "kinem.root", const char* output_file = "simulation.root"){

  MyRandom *RndmPtr = new MyRandom(input_file,seed);
  delete gRandom;
  gRandom = RndmPtr;
  if(RndmPtr->GetFlag()) cout << "There is no file named " << input_file << " in the chosen directory" << endl;
  else cout << "File " << input_file << " was found, beginning the operations" << endl;

  TFile hfile(output_file,"RECREATE");

  double X,Y,Z;
  int mult;

  TTree* Tree = new TTree("T","Tree with 3 branches");
  TClonesArray* HitsOnL1 = new TClonesArray("MyPoint",100);
  TClonesArray& L1Hit = *HitsOnL1;
  TClonesArray* HitsOnL2 = new TClonesArray("MyPoint",100);
  TClonesArray& L2Hit = *HitsOnL2;

  MyPoint* Point = new MyPoint();
  MyVertex* Vertex = new MyVertex();


  Tree->Branch("VertMult",&Vertex);
  Tree->Branch("Hits on Layer 1",&HitsOnL1);
  Tree->Branch("Hits on Layer 2",&HitsOnL2);

  Tree->SetAutoSave(0);

  //Generatori di vertice, inizio con la generazione di un vertice con molteplicità estratta da kinem.root
  //Ancora da inserire: molteplicità fissa e estratta da distribuzione uniforme

  MyPoint* Hit = new MyPoint(); //Points used to store the true position of where the particles hit the detectors
  MyParticle* Particle = new MyParticle();
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
      cout << "Generated vertex #" << i << " = (" << Vertex->GetPoint()->GetX() << ", " <<
                                                     Vertex->GetPoint()->GetY() << ", " <<
                                                     Vertex->GetPoint()->GetZ() << ");" << endl;
      cout << "Random multiplicity = " << mult << ";" << endl;
    #endif

    //Multiplicity loop
    for(int j = 0; j < mult; j++){
      //Particle generation
      Particle->SetTheta(RndmPtr->RndmTheta());
      Particle->SetPhi(RndmPtr->Uniform(0.,2.*TMath::Pi()));

      //Particle transport
      //Beam pipe interaction
      MyPoint* tmpPoint = new MyPoint(Vertex->GetX(),Vertex->GetY(),Vertex->GetZ());
//      tmpPoint = Vertex->GetPoint();
      *Hit = BeamPipe.Hit(tmpPoint, Particle);
      #if DEBUG == TRUE
        cout << "Hit position on the beam pipe = (" << Hit->GetX() << ", " <<
                                                       Hit->GetY() << ", " <<
                                                       Hit->GetZ() << "); Radius of the position = " << 
                                                       Hit->GetRadius() << endl;
      #endif


      //First layer interaction
      *Hit = Layer1.Hit(Hit, Particle);
      new(L1Hit[j])MyPoint(Hit->GetX(),Hit->GetY(),Hit->GetZ());
      #if DEBUG == TRUE
        cout << "Hit position on the first detector layer = (" << Hit->GetX() << ", " <<
                                                                  Hit->GetY() << ", " <<
                                                                  Hit->GetZ() << "); Radius of the position = " << 
                                                                  Hit->GetRadius() << endl;
      #endif

      //Second layer interaction
      *Hit = Layer2.Hit(Hit, Particle);
      new(L2Hit[j])MyPoint(Hit->GetX(),Hit->GetY(),Hit->GetZ());
      #if DEBUG == TRUE
        cout << "Hit position on the first detector layer = (" << Hit->GetX() << ", " <<
                                                                  Hit->GetY() << ", " <<
                                                                  Hit->GetZ() << "); Radius of the position = " << 
                                                                  Hit->GetRadius() << endl;
      #endif



      #if DEBUG == TRUE
        printf("Evento %d - moltepl: %d - interazione: %d\n",i,mult,j+1);
        printf("x= %f ; y= %f; z= %f \n",Vertex->GetPoint()->GetX(),Vertex->GetPoint()->GetY(),Vertex->GetPoint()->GetZ());
        printf("Entries nel TClonesArray1: %d\n",HitsOnL1->GetEntries());
        MyPoint *tst1=(MyPoint*)HitsOnL1->At(j);
        std::cout<<"Hit on L1 "<<j<<") x, y, z = "<<tst1->GetX()<<"; "<<tst1->GetY()<<"; "<<tst1->GetZ()<<std::endl;
        printf("Entries nel TClonesArray2: %d\n",HitsOnL2->GetEntries());
        MyPoint *tst2=(MyPoint*)HitsOnL2->At(j);
        std::cout<<"Hit on L2 "<<j<<") x, y, z = "<<tst2->GetX()<<"; "<<tst2->GetY()<<"; "<<tst2->GetZ()<<std::endl;
      #endif


    }

    Tree->Fill();

    HitsOnL1->Clear();
    HitsOnL2->Clear();

  }

  hfile.Write();
  hfile.Close();


//NECESSARIO AGGIUNGERE SMEARING E MULTISCATTERING
//AGGIUNGERE ANCHE MOTLEPLICITA' FISSA E DA DISTRIBUZIONE UNIFORME
//OPPORTUNO INIZIARE ANCHE LA RICOSTRUZIONE E L'ESTRAZIONE DELLE COSE DAL TREE


}