#ifndef MYTRACKLET_H
#define MYTRACKLET_H

#include "TObject.h"
#include "MySignal.h"

class MyTracklet : public TObject{

  public:
    MyTracklet();
    MyTracklet(double R1, double R2, double Z1, double Z2);
    MyTracklet(MySignal* InnerSignal, MySignal* OuterSignal);
    MyTracklet(const MyTracklet& source);
    virtual ~MyTracklet();

    MyTracklet& operator=(const MyTracklet& source);

    double Intersection();

    double GetR1() const {return dmR1;}
    double GetR2() const {return dmR2;}
    double GetZ1() const {return dmZ1;}
    double GetZ2() const {return dmZ2;}

    void SetR1(double R1) {dmR1 = R1;}
    void SetR2(double R2) {dmR2 = R2;}
    void SetZ1(double Z1) {dmZ1 = Z1;}
    void SetZ2(double Z2) {dmZ2 = Z2;}



  private:
    double dmR1;
    double dmR2;
    double dmZ1;
    double dmZ2;


  ClassDef(MyTracklet,1)

};
#endif