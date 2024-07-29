#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "Riostream.h"

#include "MyParticle.h"
#include "MyPhysics.h"
#include "MyPoint.h"
#include "MyRandom.h"
#include "MyVertex.h"

#define TRUE 1
#define FALSE 0
#define DEBUG TRUE

using namespace std;

//IMPORTANT: DISTANCES ARE MEASURED IN CENTIMETRES IN THIS SIMULATION

void Test(int N_exp = 2, unsigned int seed = 69420, const char* input_file = "kinem.root"){

  double X,Y,Z;
  int mult;
//  const char* input_file = "kinem.root";

  MyRandom *RndmPtr = new MyRandom(input_file,seed);
  delete gRandom;
  gRandom = RndmPtr;
  if(RndmPtr->GetFlag()) cout << "There is no file named " << input_file << " in the chosen directory" << endl;
  else cout << "File " << input_file << " was found, beginning the operations" << endl;

  MyPoint* Point = new MyPoint();

//  #if DEBUG == TRUE
//    Point->SetX(RndmPtr->Gaus(0.,0.053));
//    Point->SetY(RndmPtr->Gaus(0.,0.0001));
//    Point->SetZ(RndmPtr->Gaus(0.,0.0001));
//    cout << "Punto generato casualmente, gettato dal punto = (" << Point->GetX() << ", " 
//                                                                << Point->GetY() << ", " 
//                                                                << Point->GetZ() << ");" << endl;
//  #endif

  MyVertex* Vertex = new MyVertex();

//  #if DEBUG == TRUE
//    Vertex->SetPoint(Point);
//    MyPoint* PointRec = new MyPoint();
//    PointRec = Vertex->GetPoint();
//    cout << PointRec << endl;
//    cout << "Punto generato casualmente, gettato dal vertice = (" << PointRec->GetX() << ", " 
//                                                                  << PointRec->GetY() << ", " 
//                                                                  << PointRec->GetZ() << ");" << endl;
//
//    MyVertex* Vertex2 = new MyVertex();
//    MyPoint* PointRec2 = new MyPoint();
//    Vertex2->SetPoint(Point);
//    PointRec2 = Vertex2->GetPoint();
//    cout << "Punto copiato da vertex a vertex2 = (" << Vertex2->GetPoint()->GetX() << ", " 
//                                                    << Vertex2->GetPoint()->GetY() << ", " 
//                                                    << Vertex2->GetPoint()->GetZ() << ");" << endl;
//  #endif

//Generatori di vertice, inizio con la generazione di un vertice con molteplicità estratta da kinem.root
//Ancora da inserire: molteplicità fissa e estratta da distribuzione uniforme


//  #if DEBUG == TRUE
//    int generated_mult;
//    for(int i = 0; i < 10; i++){
//      generated_mult = RndmPtr->RndmMult();
//      cout << "Multiplicity extracted from the given distribution = " << generated_mult << endl;
//      }
//  #endif

//Generazione delle varie particelle di ogni vertice

//  #if DEBUG == TRUE
//    double generated_theta;
//    for(int i = 0; i < 10; i++){
//      generated_theta = RndmPtr->RndmTheta();
//      cout << "Theta extracted from the given distribution = " << generated_theta << endl;
//    }
//  #endif

  MyPoint* Hit = new MyPoint(); //Points used to store the true position of where the particles hit the detectors
  MyParticle* Particle = new MyParticle();
  MyPhysics BeamPipe(.03,.54);
  MyPhysics Layer1(.04,.27);
  MyPhysics Layer2(.07,.27);

  for(int i = 0; i < N_exp; i++){
    //Vertex generation
    Point->SetX(RndmPtr->Gaus(0.,.0001));
    Point->SetY(RndmPtr->Gaus(0.,.0001));
    Point->SetZ(RndmPtr->Gaus(0.,.053));
    Vertex->SetPoint(Point);
    Vertex->SetMult(RndmPtr->RndmMult());

    #if DEBUG == TRUE
      cout << "Generated vertex #" << i+1 << " = (" << Vertex->GetX() << ", " <<
                                                       Vertex->GetY() << ", " <<
                                                       Vertex->GetZ() << ");" << endl;
      cout << "Random multiplicity = " << Vertex->GetMult() << ";" << endl;
    #endif

    for(int j = 0; j < Vertex->GetMult(); j++){
      //Particle generation
      Particle->SetTheta(RndmPtr->RndmTheta());
      Particle->SetPhi(RndmPtr->Uniform(0.,2.*TMath::Pi()));

      #if DEBUG == TRUE
        cout << "Random theta = " << Particle->GetTheta() << ";" << endl;
        cout << "Random phi = " << Particle->GetPhi() << ";" << endl;
      #endif

      //Particle transport
      //Beam pipe interaction
      MyPoint* tmpPoint = new MyPoint();
      tmpPoint = Vertex->GetPoint();
      tmpPoint->SetZ(0.);
      cout << tmpPoint->GetZ() << endl;
      *Hit = BeamPipe.Hit(Vertex->GetPoint(), Particle);
      #if DEBUG == TRUE
        cout << "Hit position on the beam pipe = (" << Hit->GetX() << ", " <<
                                                       Hit->GetY() << ", " <<
                                                       Hit->GetZ() << "); Radius of the position = " << 
                                                       TMath::Sqrt(Hit->GetX()*Hit->GetX() + 
                                                                   Hit->GetY()*Hit->GetY()) << endl;
      #endif


      //First layer interaction
      *Hit = Layer1.Hit(Hit, Particle);
      #if DEBUG == TRUE
        cout << "Hit position on the first detector layer = (" << Hit->GetX() << ", " <<
                                                                  Hit->GetY() << ", " <<
                                                                  Hit->GetZ() << "); Radius of the position = " << 
                                                                  TMath::Sqrt(Hit->GetX()*Hit->GetX() + 
                                                                              Hit->GetY()*Hit->GetY()) << endl;
      #endif

      //Second layer interaction
      *Hit = Layer2.Hit(Hit, Particle);
      #if DEBUG == TRUE
        cout << "Hit position on the first detector layer = (" << Hit->GetX() << ", " <<
                                                                  Hit->GetY() << ", " <<
                                                                  Hit->GetZ() << "); Radius of the position = " << 
                                                                  TMath::Sqrt(Hit->GetX()*Hit->GetX() + 
                                                                              Hit->GetY()*Hit->GetY()) << endl;
      #endif

//DA FARE: SALVARE GLI HIT SUI DUE LAYER SU UN QUALCHE CONTAINER (TTREE? TCLONESARRAY?) E PROVARE A FARE
//UNA RUDIMENTALE RICOSTRUZIONE, SERVE LO SMEARING!!

    }
  }
}