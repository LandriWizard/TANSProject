#include "Riostream.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TMath.h"
#include "TStopwatch.h"
#include "TTree.h"

#include "MyParticle.h"
#include "MyPhysics.h"
#include "MyPoint.h"
#include "MyRandom.h"
#include "MyRunningWindow.h"
#include "MySignal.h"
#include "MyTracklet.h"
#include "MyVertex.h"

#define FALSE 0
#define TRUE 1
#define PRINT_EVENT FALSE
#define DEBUG FALSE
#define LOGGING FALSE
#define SORTING_DEBUG FALSE

using namespace std;
using std::vector;

//IMPORTANT: LENGHTS ARE IN CM, ANGLES IN RAD

void Reconstruction(double window_size = 0.35, double window_step = 0.175, const char* input_file = "simulation.root", const char* log_file = "reconstruction_log.txt"){

//Stopwatch declaration and start
  TStopwatch Clock;
  Clock.Start();

  //physical quantities
  double multiscattering_angle = 0.0012;
  double inner_radius = 4.;
  double outer_radius = 7.;

//Not actually used?
//  double delta_phi = TMath::ASin(3./7.*TMath::Sin(multiscattering_angle)); //I would like to use it in the while loop to "slice" the azimuth angle, cannot because its not a divisor of 2Pi
//  int slice_number = 2*TMath::Pi()/delta_phi + 1; //number of azimuth angle slices
//  double real_delta_phi = 2*TMath::Pi()/(double)slice_number; //actual value used to divide the azimuth angle, is a divisor of 2Pi; differs in the order ~10^-8 from delta_phi

  double reconstructed_z = 0;;
  double residual_z;

  //Declaring "auxiliary" objects
  MyVertex* Vertex = new MyVertex();
  MyTracklet* Tracklet = new MyTracklet(inner_radius,outer_radius,0.,0.);
  MyRunningWindow* RunningWindow = new MyRunningWindow(window_size,window_step);
  
  //Declaring TClonesArray
  TClonesArray *HitsL1 = new TClonesArray("MySignal",70);
  TClonesArray *HitsL2 = new TClonesArray("MySignal",70);
 
  //Declaring histograms
  TH1D* hResidual = new TH1D("Residual","Residual distribution",200,-1000.,1000.);
  hResidual->GetXaxis()->SetTitle("Residual [#mum]");
  hResidual->GetYaxis()->SetTitle("Counts");
  hResidual->SetLineColor(kBlack);

  //Opening input file
  TFile infile(input_file);
  if(infile.IsZombie()){
		cout << "There is no file named " << input_file << " in the chosen directory" << endl;
		return;
	}

  #if LOGGING
    //Opening log file
    ofstream lfile(log_file);
  #endif

  #if SORTING_DEBUG
    ofstream sfile("sorting_debug.txt");
  #endif

  //Reading the multiplicity generation method
  TObject* Multiplicity_Generation = (TObject*)infile.Get("Multiplicity_Generation");
  int N_mult = Multiplicity_Generation->GetUniqueID();

  vector<double> studied_multiplicities; //declared double to be used in the TGraph
  double std_studied_multiplicities[12] = {3,5,7,9,11,15,20,30,40,50,60,70};

//Given distribution
  if(N_mult == 0){
    for(int i = 0; i < 12; i++) studied_multiplicities.push_back(std_studied_multiplicities[i]);
  }
//Constant value
  else if(N_mult > 0){
    studied_multiplicities.push_back(N_mult);
  }
//Uniform distribution
  else{
    for(int i = 0; i < 12; i++){
      if(std_studied_multiplicities[i] < -1.*N_mult - 1.) studied_multiplicities.push_back(std_studied_multiplicities[i]);
    }
    studied_multiplicities.push_back(-1.*N_mult - 1.);
  }

  const int dim_mult = studied_multiplicities.size();

  //Reading TTree and branch
  TTree *tree = (TTree*)infile.Get("T");
  TBranch *b1 = tree->GetBranch("VertMult");
  TBranch *b2 = tree->GetBranch("HitsL1");
  TBranch *b3 = tree->GetBranch("HitsL2");

  //Defining addresses for data reading
  b1->SetAddress(&Vertex);
  b2->SetAddress(&HitsL1);
  b3->SetAddress(&HitsL2);



  //Loop ont the tree entries
  for(int ev = 0; ev < tree->GetEntries(); ev++){
    tree->GetEvent(ev);

    //Declaring vectors to store the reconstructed z values
    vector<double> reconstructed_z_values;

    if(ev%100000 == 0) cout << "Event #" << ev << endl;

    #if LOGGING
      lfile << "\nConsidering event #" << ev << endl;
      int vertex_counter = 0;
    #endif


    #if PRINT_EVENT
      cout << "Evento " << ev << "; Molteplicita = " << Vertex->GetMult() << endl;
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


    for(int i = 0; i < HitsL2->GetEntries(); i++){ //for cycle on 2nd layer

      #if LOGGING
        bool signal_vertex_check = 0; //bool to check if we find a vertex for each signal on layer 1
      #endif

      MySignal* InteractionOnLayer2 = (MySignal*)HitsL2->At(i);

      Tracklet->SetZ2(InteractionOnLayer2->GetZ()); //insertion of Z2 in the tracklet

      for(int j = 0; j < HitsL1->GetEntries(); j++){  //for cycle on 1st layer
        MySignal* InteractionOnLayer1 = (MySignal*)HitsL1->At(j);


        if(TMath::Abs(InteractionOnLayer2->GetPhi() - InteractionOnLayer1->GetPhi()) < 0.005/*6.*real_delta_phi*/){
          Tracklet->SetZ1(InteractionOnLayer1->GetZ()); //insertion of Z1 in the tracklet
          reconstructed_z = Tracklet->Intersection();
          reconstructed_z_values.push_back(reconstructed_z);

          #if LOGGING
            signal_vertex_check = 1;
            if(signal_vertex_check){
              lfile << "Reconstructed vertex Z = " << reconstructed_z << " thanks to hit #" << j << " on the 1st detector layer and to hit #" 
                                                                                            << i << " on the 2nd detector layer" << endl;
              if(InteractionOnLayer1->GetFlag() == InteractionOnLayer2->GetFlag())  lfile << "Vertex reconstructed using the same particles" << endl;
              else lfile << "Vertex reconstructed using particle " << InteractionOnLayer1->GetFlag() << " on the 1st layer and particle "
                                                                   << InteractionOnLayer2->GetFlag() << " on the second layer" << endl;
            }
            vertex_counter++;
          #endif

        }

      }//closing for cycle on 1st layer

      #if LOGGING
        if(!signal_vertex_check) lfile << "!! NOT FOUND ANY VERTEX FOR HIT ON LAYER 2 #" << i << endl;
        if(signal_vertex_check){
          lfile << "Real vertex z cohordinate = " << Vertex->GetZ() << endl;
          lfile << "Residual = " << Vertex->GetZ() - reconstructed_z << endl;
        }
      #endif

    }//closing for cycle on 2nd layer

    #if LOGGING
      lfile << "Number of hits on the 1st layer = " << HitsL1->GetEntries() << endl; 
      lfile << "Number of hits on the 2nd layer = " << HitsL2->GetEntries() << endl;
      lfile << "Found " << vertex_counter << " intersections" << endl;
    #endif

    //running window reconstruction
    //sorting the vector of reconstructed z values for the running window

    #if SORTING_DEBUG
      sfile << "Event #" << ev << endl;
      sfile << "Multiplicity = " << Vertex->GetMult() << "; #reconstructed z = " << reconstructed_z_values.size() << endl;
      sfile << "Before sorting" << endl;
      for(int i = 0; i < (int)(reconstructed_z_values.size()); i++){
        if(i == 0) sfile << "reconstructed_z_values = {" << reconstructed_z_values[i] << ", ";
        else if(i < (int)(reconstructed_z_values.size() - 1)) sfile << reconstructed_z_values[i] << ", ";
        else sfile << reconstructed_z_values[i] << "}" << endl;
      }
    #endif

    sort(reconstructed_z_values.begin(),reconstructed_z_values.end());  

    #if SORTING_DEBUG
      sfile << "After sorting" << endl;
      for(int i = 0; i < (int)(reconstructed_z_values.size()); i++){
        if(i == 0) sfile << "reconstructed_z_values = {" << reconstructed_z_values[i] << ", ";
        else if(i < (int)(reconstructed_z_values.size() - 1)) sfile << reconstructed_z_values[i] << ", ";
        else sfile << reconstructed_z_values[i] << "}" << endl;
      }
    #endif

//DA IMPLEMENTARE BENE
    bool reconstructable = 1;
    double reconstructed_vertex = RunningWindow->running_window(reconstructed_z_values,reconstructable);

    #if DEBUG
      cout << "Reconstructed flag set to " << reconstructable << " in the macro" << endl;
    #endif

    if(!reconstructable){
      #if DEBUG
        cout << "Running window reconstruction not possible (reconstructed_flag set to << " << reconstructable << ") in the macro, doubling the step and size values" << endl;
      #endif
      RunningWindow->SetSize(2.*window_size);
      RunningWindow->SetStep(2.*window_step);
      reconstructed_vertex = RunningWindow->running_window(reconstructed_z_values,reconstructable);
      if(!reconstructable){
        #if DEBUG
          cout << "!Running window reconstruction still not possible (reconstructed_flag set to << " << reconstructable << ") in the macro!" << endl;
        #endif
      }      
    }

    if(reconstructable){
      residual_z = Vertex->GetZ() - reconstructed_vertex;
      hResidual->Fill(residual_z*1.e4); //filling the histogram with the residual in micrometers
    }

//    cout << "Reconstructed vertex = " << reconstructed_vertex << " cm; Real vertex = " << Vertex->GetZ() << " cm; residual = " << 1.e4*(reconstructed_vertex - Vertex->GetZ()) << " um" << endl;

    RunningWindow->SetSize(window_size);  //resetting the values of the running window
    RunningWindow->SetStep(window_step);  //resetting the values of the running window

  }

//TCanvas declaration
  TCanvas *cResidual = new TCanvas("Residual","Residual",800,600);
  hResidual->Draw();


  #if LOGGING
    //closing log file
    lfile.close();
  #endif

  #if SORTING_DEBUG
    sfile.close();
  #endif

//delete pointers
  delete Vertex;
  delete Tracklet;


//Clock stop and time print
  Clock.Stop();
  Clock.Print();

}