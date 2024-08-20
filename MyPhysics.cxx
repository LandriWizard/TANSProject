#include "MyPhysics.h"

#include "TMath.h"
#include "Riostream.h"

ClassImp(MyPhysics)

////////////////////////////////////////////////
// Class containing all the physics phenomena //
// that we need to simulate                   //
////////////////////////////////////////////////

//___________________________________________________________________________
MyPhysics::MyPhysics(): TObject(),
  dmR(0.),
  dmH(0.),
  dmScatteringTheta(0.),
  dmSmearingZ(0.),
  dmSmearingRPhi(0.){
  //Default constructor
}

//___________________________________________________________________________
MyPhysics::MyPhysics(double R, double H, double scattering_theta, double smearing_z, double smearing_rphi): TObject(),
  dmR(R),
  dmH(H),
  dmScatteringTheta(scattering_theta),
  dmSmearingZ(smearing_z),
  dmSmearingRPhi(smearing_rphi){
  //Standard constructor
}

//___________________________________________________________________________
MyPhysics::MyPhysics(const MyPhysics& source): TObject(source)
{
  dmR = source.dmR;
  dmH = source.dmH;
  dmScatteringTheta = source.dmScatteringTheta;
  dmSmearingZ = source.dmSmearingZ;
  dmSmearingRPhi = source.dmSmearingRPhi;
  //copy constructor  
}

//___________________________________________________________________________
MyPhysics::~MyPhysics()
{
  //Default destructor
}

//___________________________________________________________________________
MyPhysics& MyPhysics::operator=(const MyPhysics& source){
    if(this == &source) return *this;
    this->~MyPhysics();
    new(this) MyPhysics(source);
    return *this;
    //Copy operator
}

//Straight line transport implementation
MyPoint MyPhysics::Transport(MyPoint* Point, MyParticle* Particle){
  double x0 = Point->GetX();
  double y0 = Point->GetY();
  double z0 = Point->GetZ();

  double c1 = TMath::Sin(Particle->GetTheta())*TMath::Cos(Particle->GetPhi());
  double c2 = TMath::Sin(Particle->GetTheta())*TMath::Sin(Particle->GetPhi());
  double c3 = TMath::Cos(Particle->GetTheta());

  double delta = (x0*c1 + y0*c2)*(x0*c1 + y0*c2) - (c1*c1 + c2*c2)*(x0*x0 + y0*y0 - dmR*dmR);
  double tminus = (-(x0*c1 + y0*c2) - TMath::Sqrt(delta))/(c1*c1 + c2*c2);
  double tplus = (-(x0*c1 + y0*c2) + TMath::Sqrt(delta))/(c1*c1 + c2*c2);

  double t;
  if(tplus > 0.) t = tplus;
  else t = tminus;

  double x = x0 + c1*t;
  double y = y0 + c2*t;
  double z = z0 + c3*t;

  MyPoint hit(x,y,z);

  return hit;

}

//Multiple Coulomb scattering implementation
MyParticle MyPhysics::MultipleScattering(MyParticle* Particle){


    //deviazione dalla direzione di entrata

    double ThetaP = 0.;
    do {ThetaP = gRandom -> Gaus(0.,dmScatteringTheta);} while(ThetaP<0.); //vogliamo theta >=0 (il Do serve per farlo eseguire almeno una volta)
    double PhiP = gRandom -> Rndm()*2.*TMath::Pi();

    double theta = Particle->GetTheta();
    double phi = Particle->GetPhi();

    double mr[3][3];
    mr[0][0] = - TMath::Sin(phi); 
    mr[1][0] = TMath::Cos(phi); 
    mr[2][0] = 0.;
    mr[0][1] = - TMath::Cos(phi)*TMath::Cos(theta); 
    mr[1][1] = - TMath::Sin(phi)*TMath::Cos(theta);
    mr[2][1] = TMath::Sin(theta); 
    mr[0][2] = TMath::Cos(phi)*TMath::Sin(theta);
    mr[1][2] = TMath::Sin(phi)*TMath::Sin(theta);
    mr[2][2] = TMath::Cos(theta);
//    double rotation_matrix[3][3];
//    rotation_matrix[0][0] = -TMath::Sin(phi);
//    rotation_matrix[1][0] = TMath::Cos(phi);
//    rotation_matrix[2][0] = 0.;
//    rotation_matrix[0][1] = -TMath::Cos(phi)*TMath::Cos(theta);
//    rotation_matrix[1][1] = -TMath::Sin(phi)*TMath::Cos(theta);
//    rotation_matrix[2][1] = TMath::Sin(theta);
//    rotation_matrix[0][2] = TMath::Cos(phi)*TMath::Sin(theta);
//    rotation_matrix[1][2] = TMath::Sin(phi)*TMath::Sin(theta);
//    rotation_matrix[2][2] = TMath::Cos(theta);

    double scat[3]; //coordinate cartesiane della deviazione (cioè della particella dopo il multiscattering, rispetto alla direzione iniziale)
    scat[0] = TMath::Cos(PhiP)*TMath::Sin(ThetaP);
    scat[1] = TMath::Sin(PhiP)*TMath::Sin(ThetaP);
    scat[2] = TMath::Cos(ThetaP);
//    double cdp[3]; //director cosines in the primed reference system
//    cdp[0] = TMath::Cos(phiP)*TMath::Sin(thetaP);
//    cdp[1] = TMath::Sin(phiP)*TMath::Sin(thetaP);
//    cdp[2] = TMath::Cos(thetaP);

    double final_dir[3];

    for (int i = 0; i < 3; i++){
        final_dir[i] = 0.;
        for (int j = 0; j < 3; j++){
            final_dir[i]+=mr[i][j]*scat[j];
        }
    }
//    double temp_cd[3];
//    for(int i = 0; i < 3; i++){
//      temp_cd[i] = 0.;
//      for(int j = 0; j < 3; j++){
//        temp_cd[i] += rotation_matrix[i][j]*cdp[j];
//      }
//    }
//

    //Coordinate della particella dopo il multiscattering nel sistema di riferimento di partenza 
    double final_theta = TMath::ACos(final_dir[2]);
    double final_phi;
    
    //final_phi e' sempre ben definito poiche' se sin(theta) fosse 0 dovremmo avere un punto del tipo (0,0,z) ossia sull'asse del fascio -> non possibile 
    if(final_dir[1]>=0.) final_phi=TMath::ACos(final_dir[0]/(TMath::Sin(final_theta)));
    else final_phi=2.*TMath::Pi()-TMath::ACos(final_dir[0]/(TMath::Sin(final_theta)));
//    if(temp_cd[1] >= 0.) scattered_phi = TMath::ACos(temp_cd[0])/TMath::Sin(scattered_theta); 
//    else scattered_phi = 2.*TMath::Pi() - TMath::ACos(temp_cd[0])/TMath::Sin(scattered_theta);
    return MyParticle(final_theta,final_phi);


//    double theta = Particle->GetTheta();
//    double phi = Particle->GetPhi();
//
//    double thetaP = TMath::Abs(gRandom->Gaus(0.,dmScatteringTheta));
//    double phiP = gRandom->Uniform(0.,2.*TMath::Pi());
//
//    double rotation_matrix[3][3];
//    rotation_matrix[0][0] = -TMath::Sin(phi);
//    rotation_matrix[1][0] = TMath::Cos(phi);
//    rotation_matrix[2][0] = 0.;
//    rotation_matrix[0][1] = -TMath::Cos(phi)*TMath::Cos(theta);
//    rotation_matrix[1][1] = -TMath::Sin(phi)*TMath::Cos(theta);
//    rotation_matrix[2][1] = TMath::Sin(theta);
//    rotation_matrix[0][2] = TMath::Cos(phi)*TMath::Sin(theta);
//    rotation_matrix[1][2] = TMath::Sin(phi)*TMath::Sin(theta);
//    rotation_matrix[2][2] = TMath::Cos(theta);
//
//    double cdp[3]; //director cosines in the primed reference system
//    cdp[0] = TMath::Cos(phiP)*TMath::Sin(thetaP);
//    cdp[1] = TMath::Sin(phiP)*TMath::Sin(thetaP);
//    cdp[2] = TMath::Cos(thetaP);
//
//    double temp_cd[3];
//    for(int i = 0; i < 3; i++){
//      temp_cd[i] = 0.;
//      for(int j = 0; j < 3; j++){
//        temp_cd[i] += rotation_matrix[i][j]*cdp[j];
//      }
//    }
//
//    double scattered_theta = TMath::ACos(temp_cd[2]);
//    double scattered_phi;
//
//    if(temp_cd[1] >= 0.) scattered_phi = TMath::ACos(temp_cd[0])/TMath::Sin(scattered_theta); 
//    else scattered_phi = 2.*TMath::Pi() - TMath::ACos(temp_cd[0])/TMath::Sin(scattered_theta);
//
//    return MyParticle(scattered_theta,scattered_phi);
}

//Smearing implementation
MySignal MyPhysics::SmearingOn(MySignal* Signal){
  int flag = Signal->GetFlag();
  double r = Signal->GetR();
  double z_rec = Signal->GetZ() + gRandom->Gaus(0.,dmSmearingZ);
  double phi_rec = Signal->GetPhi() + gRandom->Gaus(0.,dmSmearingRPhi)/r;
  if(phi_rec < 0.) phi_rec += 2.*TMath::Pi();
  else if(phi_rec > 2.*TMath::Pi()) phi_rec -= 2.*TMath::Pi();

  return MySignal(r,z_rec,phi_rec,flag);
}
