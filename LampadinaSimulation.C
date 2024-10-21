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
#define MULTISCATTERING_DEBUG FALSE
#define INDEX_DEBUG FALSE
#define NOISE_DEBUG FALSE

using namespace std;

//IMPORTANT: DISTANCES ARE MEASURED IN CENTIMETRES IN THIS SIMULATION
//MULTIPLICITY FLAG VALUES: 1 FOR EXTRACTION FROM GIVEN DISTIBUTION, 2 FOR CONSTANT VALUE, 3 FOR UNIFORM DISTRIBUTION
//ZGEN FLAG VALUES: 1 FOR GAUSS, 2 FOR UNIFORM
//MULTISCATTERING FLAG VALUES: 0 FOR NO SCATTERING, 1 FOR SCATTERING
//SMEARING FLAG VALUES: 0 FOR NO SMEARING, 1 FOR SMEARING
//N_NOISE IS THE NUMBER OF NOISE HITS TO BE GENERATED

void Simulation(int N_exp = 1e6, unsigned int seed = 69420, int zgen_flag = 1, int multiplicity_flag = 1, int multiscattering_flag = 1, int smearing_flag = 1, int N_noise = 0, const char* input_file = "kinem.root", const char* output_file = "simulation.root"){

  MyRandom *RndmPtr = new MyRandom(input_file,seed);
  delete gRandom;
  gRandom = RndmPtr;
  if(RndmPtr->GetFlag()) cout << "There is no file named " << input_file << " in the chosen directory" << endl;
  else cout << "File " << input_file << " was found, beginning the operations" << endl;

//Stopwatch declaration and start
  TStopwatch Clock;
  Clock.Start();

  double X,Y,Z;
  int mult;
  int noise_count = 0;
  int layer1_array_position = 0;
  int layer2_array_position = 0;

//Physic constants
  double multiscattering_angle = 0.0012/TMath::Sqrt(2); //mrad
  double smearing_z = 0.012; //cm 
  double smearing_rphi = 0.003; //cm

//Auxialiary objects
  MyPoint* Point = new MyPoint();
  MyVertex* Vertex = new MyVertex();
  MyPoint* Hit = new MyPoint(); //Points used to store the true position of where the particles hit the detectors
  MyParticle* Particle = new MyParticle(); //The particle that will be transported, contains the theta and phi of the trajectory

  MyPhysics BeamPipe(3.,54.,multiscattering_angle,smearing_z,smearing_rphi);
  MyPhysics Layer1(4.,27.,multiscattering_angle,smearing_z,smearing_rphi);
  MyPhysics Layer2(7.,27.,multiscattering_angle,smearing_z,smearing_rphi);

//Object used to store the multiplicity/z generation method, used to differentiate the histograms in the Reconstruction
  TObject Multiplicity_Generation;
  TObject Z_Generation;
  Z_Generation.SetUniqueID(zgen_flag);

  //Functors definition
  //Multiplicity extraction functor
  int dim = 0;  //Dimension of the TClonesArray storing the hits
  int N;
  int (MyRandom::*RndmMult) (int);
  //MULTIPLICITY FLAG VALUES: 1 FOR EXTRACTION FROM GIVEN DISTIBUTION, 2 FOR CONSTANT VALUE, 3 FOR UNIFORM DISTRIBUTION
  switch (multiplicity_flag){
  case 1:
    cout << "Extracting the multiplicity from a given distribution" << endl;
    RndmMult = &MyRandom::RndmMult_FromDistribution;
    dim = 70;
    Multiplicity_Generation.SetUniqueID(0);
    break;
  case 2:
    cout << "How many particles do you want for each vertex? ";
    cin >> N;
    cout << "You chose to have " << N << " particles for each vertex" << endl;
    RndmMult = &MyRandom::RndmMult_Constant;
    dim = N;
    Multiplicity_Generation.SetUniqueID(N);
    break;
  case 3: 
    cout << "What is the wanted maximum number of particles for each vertex? ";
    cin >> N;
    cout << "You chose to have a maximum of " << N << " particles for each vertex" << endl;
    RndmMult = &MyRandom::RndmMult_Uniform;
    dim = N/2 + 1;
    Multiplicity_Generation.SetUniqueID(-N);
    break;
  default:
    cout << "multiplicity_flag value choice not valid. 1 for extraction from given distribution, 2 for constant value, 3 for uniform distribution" << endl;
    cout << "Default choice: extracting the multiplicity from a given distribution" << endl;
    RndmMult = &MyRandom::RndmMult_FromDistribution;
    dim = 70;
    Multiplicity_Generation.SetUniqueID(0);
    break;
  }
  //Multiscattering functor
  MyParticle (MyPhysics::*ScatteringFunc)(MyParticle*);
  switch (multiscattering_flag)
  {
  case 0:
    cout << "Multiple scattering not activated" << endl;
    ScatteringFunc = &MyPhysics::NoScattering;
    break;
    case 1:
    cout << "Multiple scattering activated" << endl;
    ScatteringFunc = &MyPhysics::MultipleScattering;
    break;
  default:
    cout << "multiscattering_flag value choice not valid. 0 for no scattering, 1 for multiple scattering" << endl;
    cout << "Default choice: multiple scattering activated" << endl;
    ScatteringFunc = &MyPhysics::MultipleScattering;
    break;
  }
  //Smearing functor
  MySignal (MyPhysics::*SmearingFunc)(MySignal*);
  switch (smearing_flag)
  {
  case 0:
    cout << "Smearing not activated" << endl;
    SmearingFunc = &MyPhysics::NoSmearing;
    break;
  case 1:
    cout << "Smearing activated" << endl;
    SmearingFunc = &MyPhysics::SmearingOn;
    break;
  default:
    cout << "smearing_flag value choice not valid. 0 for no smearing, 1 for smearing" << endl;
    cout << "Default choice: smearing activated" << endl;
    SmearingFunc = &MyPhysics::SmearingOn;
    break;
  }

//Output file declaration
  TFile outfile(output_file,"RECREATE");
  Multiplicity_Generation.Write("Multiplicity_Generation");
  Z_Generation.Write("Z_Generation");

  TTree* Tree = new TTree("T","Tree with 3 branches");
  TClonesArray* HitsOnL1 = new TClonesArray("MySignal",dim);
  TClonesArray& L1Hit = *HitsOnL1;
  TClonesArray* HitsOnL2 = new TClonesArray("MySignal",dim);
  TClonesArray& L2Hit = *HitsOnL2;

  Tree->Branch("VertMult",&Vertex);
  Tree->Branch("HitsL1",&HitsOnL1);
  Tree->Branch("HitsL2",&HitsOnL2);

  Tree->SetAutoSave(0);

  //Events loop
  for(int i = 0; i < N_exp; i++){

    if(i%100000 == 0) cout << "Vertex #" << i << endl;

    //Vertex cohordinates generation
    Vertex->SetX(RndmPtr->Gaus(0.,0.01));
    Vertex->SetY(RndmPtr->Gaus(0.,0.01));

    switch(zgen_flag){
      case 1:
        if(i == 0)  cout << "Extracting the vertex z from a gaussian distribution" << endl;
        Vertex->SetZ(RndmPtr->Gaus(0.,5.3));
        break;
      case 2:
        if(i == 0)  cout << "Extracting the vertex z from a uniform distribution" << endl;
        Vertex->SetZ(RndmPtr->Uniform(-13.5, 13.5));
        break;
      default:
        if(i == 0){
          cout << "zgen_flag value choice not valid. 1 for gaussian distribution, 2 for uniform distribution" << endl;
          cout << "Default choice: exctacting the vertex z from the gaussian distribution" << endl;
        }
        Vertex->SetZ(RndmPtr->Gaus(0.,5.3));
        break;
    }
    //Vertex multiplicity generation
    Vertex->SetMult((RndmPtr->*RndmMult)(N));

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

      #if INDEX_DEBUG == TRUE
        cout << "Indeces at the beginning of the for cycle:" << endl;
        printf("j = %d\n",j);
        printf("j1 = %d\n",j1);
        printf("j2 = %d\n",j2);
      #endif


      //Particle generation
      Particle->SetTheta(RndmPtr->RndmTheta());
      Particle->SetPhi(RndmPtr->Uniform(0.,2.*TMath::Pi()));

      #if MULTISCATTERING_DEBUG == TRUE
        double first_theta = Particle->GetTheta();
        double first_phi = Particle->GetPhi();
      #endif

      //Particle transport
      //Beam pipe interaction
      MyPoint* tmpPoint = new MyPoint(Vertex->GetX(),Vertex->GetY(),Vertex->GetZ());
      *Hit = BeamPipe.Transport(tmpPoint, Particle);  //Particle transport
      #if DEBUG == TRUE
        cout << "Hit position on the beam pipe = (" << Hit->GetX() << ", " <<
                                                       Hit->GetY() << ", " <<
                                                       Hit->GetZ() << "); Radius of the position = " << 
                                                       Hit->GetRadiusXY() << endl;
      #endif
      *Particle = (BeamPipe.*ScatteringFunc)(Particle);
      #if MULTISCATTERING_DEBUG == TRUE
        double after_bp_theta = Particle->GetTheta();
        double after_bp_phi = Particle->GetPhi();
        cout << "Difference in theta after beam pipe = " << after_bp_theta - first_theta << endl;
        cout << "Difference in phi after beam pipe = " << after_bp_phi - first_phi << endl;
      #endif

      //First layer interaction
      *Hit = Layer1.Transport(Hit, Particle); //Particle transport
      #if DEBUG == TRUE
        cout << "Hit position on the first detector layer = (" << Hit->GetX() << ", " <<
                                                                  Hit->GetY() << ", " <<
                                                                  Hit->GetZ() << "); Radius of the position = " << 
                                                                  Hit->GetRadiusXY() << endl;
      #endif

      if(Hit->GetZ() > -1.*Layer1.GetH()/2. && Hit->GetZ() < Layer1.GetH()/2.){ //Check if Z is on the detector

        *Particle = (Layer1.*ScatteringFunc)(Particle);
        #if MULTISCATTERING_DEBUG == TRUE
        double after_l1_theta = Particle->GetTheta();
        double after_l1_phi = Particle->GetPhi();
        cout << "Difference in theta after layer 1 = " << after_l1_theta - after_bp_theta << endl;
        cout << "Difference in phi after layer 1 = " << after_l1_phi - after_bp_phi << endl;
        #endif

        MySignal* tempSignal1 = new MySignal(Hit,j);
        *tempSignal1 = (Layer1.*SmearingFunc)(tempSignal1);
        new(L1Hit[j1]) MySignal(*tempSignal1);

        //Second layer interaction
        *Hit = Layer2.Transport(Hit, Particle);
        #if DEBUG == TRUE
          cout << "Hit position on the second detector layer = (" << Hit->GetX() << ", " <<
                                                                    Hit->GetY() << ", " <<
                                                                    Hit->GetZ() << "); Radius of the position = " << 
                                                                    Hit->GetRadiusXY() << endl;
        #endif

        if(Hit->GetZ() > -1.*Layer2.GetH()/2. && Hit->GetZ() < Layer2.GetH()/2.){ //Check if Z is on the detector

          MySignal* tempSignal2 = new MySignal(Hit,j);
          *tempSignal2 = (Layer2.*SmearingFunc)(tempSignal2);
          new(L2Hit[j2]) MySignal(*tempSignal2);


          j2++;
          delete tempSignal2;

        }

        j1++;

        delete tempSignal1;

      }

      #if INDEX_DEBUG == TRUE
        cout << "Indeces at the end of the for cycle:" << endl;
        printf("j = %d\n",j);
        printf("j1 = %d\n",j1);
        printf("j2 = %d\n\n\n",j2);
      #endif

      delete tmpPoint;

    }//end of the multiplicity loop

    //Noise generation

    for(int i = 0; i < N_noise; i++){
      double radius = Layer1.GetR();
      double z = -(Layer1.GetH())/2. + (RndmPtr->Rndm()) * Layer1.GetH();
      double phi = RndmPtr->Rndm() * 2. * TMath::Pi();
      new(L1Hit[j1]) MySignal(radius, z, phi, -(i+1));
      #if NOISE_DEBUG == TRUE
        cout << "Noise hit number " << i << " generated on layer 1" << endl;
        cout << "Generated noise hit position on layer 1 r, z, phi = (" << radius << ", " <<
                                                                                z << ", " <<
                                                                                phi << ");" << endl;
      #endif
      j1++;
      radius = Layer2.GetR();
      z = -(Layer2.GetH())/2. + (RndmPtr->Rndm()) * Layer2.GetH();
      phi = RndmPtr->Rndm() * 2. * TMath::Pi();
      new(L2Hit[j2]) MySignal(radius, z, phi, -(i+1));
      #if NOISE_DEBUG == TRUE
        cout << "Noise hit number " << i << " generated on layer 2" << endl;
        cout << "Generated noise hit position on layer 2 = (" << radius << ", " <<
                                                                      z << ", " <<
                                                                      phi << ");" << endl;
      #endif
      j2++;
    }
    Tree->Fill();

    HitsOnL1->Clear();
    HitsOnL2->Clear();

  }//end of the events loop
  cout << "End of the events loop" << endl;

//Writing and saving the output file
  outfile.Write();
  outfile.Close();


//Deallocating pointers
  delete Point;
  delete Vertex;
  delete Hit;
  delete Particle;

//Stopwatch stop and time print
  Clock.Stop();
  Clock.Print();

}
