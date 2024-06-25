#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMath.h"
#include "MyPoint.h"
#include "TClonesArray.h"
#include "Riostream.h"
#include "TRandom3.h"

void GeneraTree(){
 // Apertura del file di output
  TFile hfile("htree.root","RECREATE");
  // dichiarazione del TTree
  TTree *tree = new TTree("T","TTree con 2 branches");
  // se invertissi l'ordine dovrei scrivere
  // tree->SetDirectory(&hfile);

  // Dichiarazione di un TClonesArray
  TClonesArray *ptrhits = new TClonesArray("MyPoint",100);
  TClonesArray &hits = *ptrhits;
   
  // Definizione di una struct
  typedef struct{
    double X,Y,Z;
    int mult;} VTX;
    
    
  static VTX point;

  // Dichiarazione dei 2 branch del TTree
  tree->Branch("VertMult",&point.X,"X/D:Y:Z:mult/I");
  tree->Branch("Hits",&ptrhits);
      
  for(int i=0; i<100;i++){ // loop sugli eventi

    // Genero una molteplicita e un vertice
    int numpart=0;
    while(numpart<=0){
      numpart=(int)(0.5+gRandom->Gaus(50.,20.));
    }
    point.mult=numpart;
    point.X=gRandom->Gaus(0.,0.01);
    point.Y=gRandom->Gaus(0.,0.01);
    point.Z=gRandom->Gaus(0.,5.3);

    for (int j=0; j<numpart; j++){

      // Genero un hit in modo del tutto random (dummy)
      new(hits[j])MyPoint(-5.+gRandom->Rndm()*10.,5.+gRandom->Rndm()*10,15.+gRandom->Rndm()*10.);
    }

    // Debug
    printf("Evento %d - moltepl: %d\n",i,numpart);
    printf("x= %f ; y= %f; z= %f \n",point.X,point.Y,point.Z);
    printf("Entries nel TClonesArray: %d\n",ptrhits->GetEntries());
    for (int j=0; j<hits.GetEntries(); j++){
      MyPoint *tst=(MyPoint*)ptrhits->At(j);
      std::cout<<"MyPoint "<<j<<") x, y, z = "<<tst->GetX()<<"; "<<tst->GetY()<<"; "<<tst->GetZ()<<std::endl;
    }
    // fine del debug


    tree->Fill();
    ptrhits->Clear();

  }

  // Save all objects in this file
  hfile.Write();

  // Close the file. 
  hfile.Close();


}