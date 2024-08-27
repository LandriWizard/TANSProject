#ifndef MYRUNNINGWINDOW_H
#define MYRUNNINGWINDOW_H

#include "TObject.h"

using std::vector;

class MyRunningWindow : public TObject{

  public:

    MyRunningWindow();
    MyRunningWindow(double size, double step);
		MyRunningWindow(const MyRunningWindow& source);
		virtual ~MyRunningWindow();
		MyRunningWindow& operator=(const MyRunningWindow& source);

    double GetSize() const {return dmSize;}
    double GetStep() const {return dmStep;}

    void SetSize(double size) {dmSize = size;}
    void SetStep(double step) {dmStep = step;}

    double running_window(vector<double> vector, bool &reconstructed_flag);

    private:
      double dmSize;
      double dmStep;
      double Average(vector<double> vector, double upper_bound, double lower_bound);


  ClassDef(MyRunningWindow,1)
};

#endif