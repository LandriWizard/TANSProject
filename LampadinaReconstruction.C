#include <Riostream.h>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "MyPoint.h"
#include "TMath.h"
#include "TClonesArray.h"
//using namespace std;

void Reconstruction(){
  // definizione struct
  typedef struct {
    double X,Y,Z;
    int mult;} VTX;
  static VTX point;
  // Dichiarazione TClonesArray
  TClonesArray *hits = new TClonesArray("MyPoint",100);
  //Apertura file di input
  TFile hfile("htree.root");
  //Lettura TTree  e branch
  TTree *tree = (TTree*)hfile.Get("T");
  TBranch *b1=tree->GetBranch("VertMult");
  TBranch *b2=tree->GetBranch("Hits");

  // Definizione degli indirizzi per la lettura dei dati su ttree
  b1->SetAddress(&point.X);
//  b2->SetAddress(*s);

  // loop sugli ingressi nel TTree
  for(int ev=0;ev<tree->GetEntries();ev++){
    tree->GetEvent(ev);
    std::cout<<"Evento "<<ev<<"; Molteplicita= "<<point.mult<<std::endl;
    std::cout<<"X,Y,Z = "<<point.X<<"; "<<point.Y<<"; "<<point.Z<<std::endl;
    int num=hits->GetEntries();
    std::cout<<"Numero di elementi nel TClonesArray "<<num<<std::endl;
    for (int j=0; j<num; j++){
      MyPoint *tst=(MyPoint*)hits->At(j);
      std::cout<<"MyPoint "<<j<<") x, y, z = "<<tst->GetX()<<"; "<<tst->GetY()<<"; "<<tst->GetZ()<<std::endl;
    }
  }
 
}