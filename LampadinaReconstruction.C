#include "Riostream.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TEfficiency.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TMath.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TTree.h"
#include "TROOT.h"

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

void Reconstruction(double window_size = 0.35, double window_step = 0.175, const char* input_file = "simulation.root", const char* output_file = "analysis.root"){

  gROOT->SetBatch(kTRUE);

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
 
  //Declaring vectors to store the reconstructed z values
  vector<double> reconstructed_z_values;

  //Opening input file
  TFile infile(input_file);
  if(infile.IsZombie()){
		cout << "There is no file named " << input_file << " in the chosen directory" << endl;
		return;
	}

  #if LOGGING
    //Opening log file
    ofstream lfile("reconstruction_log.txt");
  #endif

  #if SORTING_DEBUG
    ofstream sfile("sorting_debug.txt");
  #endif

  //Reading the multiplicity/z generation method
  TObject* Multiplicity_Generation = (TObject*)infile.Get("Multiplicity_Generation");
  TObject* Z_Generation = (TObject*)infile.Get("Z_Generation");
  int N_mult = Multiplicity_Generation->GetUniqueID();
  int z_gen = Z_Generation->GetUniqueID();

  vector<double> studied_multiplicities; //declared double to be used in the TGraph
  double std_studied_multiplicities[12] = {3,5,7,9,11,15,20,30,40,50,60,69};

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

  //Vectors used to study the resolution
  //VS MULTIPLICITY
  vector<double> res_mult(dim_mult,0.);
  vector<double> res_mult_error(dim_mult,0.);
  vector<double> multiplicity_error(dim_mult,0.);
  //VS VERTEX Z
  vector<double> z_values;
  vector<double> standard_z_values = {-14., -12., -10., -8., -6., -4., -2., 0., 2., 4. , 6., 8., 10., 12., 14.}; 
  for(int i = 0; i < 15; i++){
    z_values.push_back(standard_z_values[i]);
  }
  const int dim_z = z_values.size();
  vector<double> z_error(dim_z,0.);
  vector<double> res_z(dim_z,0.);
  vector<double> res_z_error(dim_z,0.);

  //Objects to study the efficiency
  TEfficiency* effMult = new TEfficiency("effMult","Efficiency vs Multiplicity;Multiplicity;#epsilon",dim_mult,0.,60.);
  TEfficiency* effZ = new TEfficiency("effZ","Efficiency vs Vertex Z;Vertex Z;#epsilon",dim_z,-16.,16.);

//Declaring histograms
  //residual histogram
  TH1D* hResidual = new TH1D("Residual","Residual distribution",200,-1000.,1000.);
  hResidual->GetXaxis()->SetTitle("Residual [#mum]");
  hResidual->GetYaxis()->SetTitle("Counts");
  hResidual->SetLineColor(kBlack);
  //residual histograms for given multiplicities, an array for compactness
  TH1D* hResidualMult[dim_mult];
  char name[100];
  char title[100];
  //array filling
  for(int i = 0; i < dim_mult; i++){
    sprintf(name,"ResidualMult%d",i);
    sprintf(title,"%f < Multiplicity < %f", studied_multiplicities[i] - .5, studied_multiplicities[i] + .5);
    hResidualMult[i] = new TH1D(name, title, 400, -1000., 1000.);
    hResidualMult[i]->GetXaxis()->SetTitle("Residual [#mum]");
    hResidualMult[i]->GetYaxis()->SetTitle("Counts");
    hResidualMult[i]->SetLineColor(kBlack);
  }
  //residual histograms for given z values
  TH1D* hResidualZ[dim_z];
  //array filling
  for(int i = 0; i < dim_z; i++){
    sprintf(name,"ResidualZ%d",i);
    sprintf(title,"%f < Vertex Z < %f", z_values[i] - 1., z_values[i] + 1.);
    hResidualZ[i] = new TH1D(name, title, 400, -1000., 1000.);
    hResidualZ[i]->GetXaxis()->SetTitle("Residual [#mum]");
    hResidualZ[i]->GetYaxis()->SetTitle("Counts");
    hResidualZ[i]->SetLineColor(kBlack);
  }

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
      #if DEBUG
        if(!reconstructable){
          cout << "!Running window reconstruction still not possible (reconstructed_flag set to << " << reconstructable << ") in the macro!" << endl;
        }
      #endif        
    }

    //Filling the efficiency histograms
    effMult->Fill(reconstructable,Vertex->GetMult());
    effZ->Fill(reconstructable,Vertex->GetZ());

    if(reconstructable){
      residual_z = Vertex->GetZ() - reconstructed_vertex;
      hResidual->Fill(residual_z*1.e4); //filling the histogram with the residual in micrometers
      //looping on multiplicity histograms
      for(int i = 0; i < dim_mult; i++){
        if(Vertex->GetMult() > studied_multiplicities[i] - .5 && Vertex->GetMult() < studied_multiplicities[i] + .5){
          hResidualMult[i]->Fill(residual_z*1.e4);
          break;
        }
      }
      //looping on z histograms
      for(int i = 0; i < dim_z; i++){
        if(Vertex->GetZ() > z_values[i] - 1. && Vertex->GetZ() < z_values[i] + 1.){
          hResidualZ[i]->Fill(residual_z*1.e4);
          break;
        }
      }
    }

    reconstructed_z_values.clear();   //clearing the vector of reconstructed z values
    RunningWindow->SetSize(window_size);  //resetting the values of the running window
    RunningWindow->SetStep(window_step);  //resetting the values of the running window

  }//closing for cycle on tree entries

//Fitting the residual histograms for given multiplicities
  TF1* fResidualMult[dim_mult];
  for(int i = 0; i < dim_mult; i++){
    if(hResidualMult[i]->GetEntries() != 0){
      fResidualMult[i] = new TF1("fResidualMult","gaus",-1000.,1000.);
      hResidualMult[i]->Fit(fResidualMult[i],"R");
      //computing the resolution
      multiplicity_error[i] = .5;
      res_mult[i] = fResidualMult[i]->GetParameter(2);
      res_mult_error[i] = fResidualMult[i]->GetParError(2);
    }
  }
  TF1* fResidualZ[dim_z];
  for(int i = 0; i < dim_z; i++){
    if(hResidualZ[i]->GetEntries() != 0){
      fResidualZ[i] = new TF1("fResidualZ","gaus",-1000.,1000.);
      hResidualZ[i]->Fit(fResidualZ[i],"R");
      //computing the resolution
      z_error[i] = 1.;
      res_z[i] = fResidualZ[i]->GetParameter(2);
      res_z_error[i] = fResidualZ[i]->GetParError(2);
    }
  }

//TCanvas and TGraphErrors declaration
  char canvas_name[100];
  char canvas_title[100];
  switch(z_gen){
    case 1:
      sprintf(canvas_name,"cResMult - Gaussian");
      sprintf(canvas_title,"Resolution vs Multiplicity - Gaussian");
      break;
    case 2:
      sprintf(canvas_name,"cResMult - Uniform");
      sprintf(canvas_title,"Resolution vs Multiplicity - Uniform");
      break;
    default:
      sprintf(canvas_name,"cResMult - Gaussian");
      sprintf(canvas_title,"Resolution vs Multiplicity - Gaussian");
      break;
  }

  TCanvas* cResMult = new TCanvas(canvas_title,canvas_name,800,600);
  cResMult->cd();
  TGraphErrors* gResMult = new TGraphErrors(dim_mult,&(studied_multiplicities[0]),&(res_mult[0]),&(multiplicity_error[0]),&(res_mult_error[0]));
  gResMult->SetTitle("Resolution;Multiplicity;Resolution [#mum]");
  gResMult->SetMarkerStyle(20);
  gResMult->SetMarkerSize(1.5);
  gResMult->SetMarkerColor(kBlack);
  gResMult->Draw("AP");
  cResMult->SaveAs("cResMult.pdf");

  switch(z_gen){
    case 1:
      sprintf(canvas_name,"cResZ - Gaussian");
      sprintf(canvas_title,"Resolution vs Vertex Z - Gaussian");
      break;
    case 2:
      sprintf(canvas_name,"cResZ - Uniform");
      sprintf(canvas_title,"Resolution vs Vertex Z - Uniform");
      break;
    default:
      sprintf(canvas_name,"cResZ - Gaussian");
      sprintf(canvas_title,"Resolution vs Vertex Z - Gaussian");
      break;
  }

  TCanvas* cResZ = new TCanvas(canvas_title,canvas_name,800,600);
  cResZ->cd();
  TGraphErrors* gResZ = new TGraphErrors(dim_z,&(z_values[0]),&(res_z[0]),&(z_error[0]),&(res_z_error[0]));
  gResZ->SetTitle("Resolution;Vertex Z;Resolution [#mum]");
  gResZ->SetMarkerStyle(20);
  gResZ->SetMarkerSize(1.5);
  gResZ->SetMarkerColor(kBlack);
  gResZ->Draw("AP");
  cResZ->SaveAs("cResZ.pdf");

  //Efficiency plots
  TCanvas* cEffMult = new TCanvas("cEffMult","Efficiency vs Multiplicity",80,80,775,500);
  cEffMult->cd();
  effMult->SetMarkerStyle(33);
  effMult->SetMarkerColor(77);
  effMult->Draw("AP");
  cEffMult->SaveAs("cEffMult.pdf");

  TCanvas* cEffZ = new TCanvas("cEffZ","Efficiency vs Vertex Z",800,600);
  cEffZ->cd();
  effZ->SetMarkerStyle(33);
  effZ->SetMarkerColor(77);
  effZ->Draw("AP");
  cEffZ->Update();
  effZ->GetPaintedGraph()->SetMinimum(0.3);
  effZ->GetPaintedGraph()->SetMaximum(1.1);
  effZ->Draw("AP");
  cEffZ->SaveAs("cEffZ.pdf");

  //Opening output file
  TFile outfile(output_file,"RECREATE");

//Writing and saving the output file
  outfile.cd();
  outfile.Write();
  gResMult->Write("Resolution vs Multiplicity");
  gResZ->Write("Resolution vs Vertex Z");
  effMult->Write("Efficiency vs Multiplicity");
  effZ->Write("Efficiency vs Vertex Z");
  outfile.Close();

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
  delete RunningWindow;

//Clock stop and time print
  Clock.Stop();
  Clock.Print();

}
