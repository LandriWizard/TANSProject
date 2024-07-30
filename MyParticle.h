#ifndef MYPARTICLE_H
#define MYPARTICLE_H


#include "TObject.h"

class MyParticle : public TObject
{

  public:

    MyParticle();
    MyParticle(double Theta, double Phi);
    MyParticle(const MyParticle& source);
    virtual ~MyParticle();
    MyParticle& operator=(const MyParticle& source);		

    double GetTheta() const {return dmTheta;} 
    double GetPhi() const {return dmPhi;}

    void SetTheta(double Theta) {dmTheta = Theta;}
    void SetPhi(double Phi) {dmPhi = Phi;}


  private:

    double dmTheta;
    double dmPhi;

  ClassDef(MyParticle,1)
};



#endif