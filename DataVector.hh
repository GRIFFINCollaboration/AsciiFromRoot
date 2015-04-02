#ifndef DATA_VECTOR_HH
#define DATA_VECTOR_HH

#include <iostream>
#include <string>

#include "TMath.h"
#include "TVectorT.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

using namespace std;

//fRebin:
//0 = sum the values  and use square root of sum of squares of errors as errors
//1 = average the values and use min/max as errors
//2 = average the values and use square root of sum of squares of errors as errors

class DataVector : public TObject
{
public:
  DataVector();
  DataVector(const char*, DataVector, int index = -999, int rebin = 0);
  DataVector(int size, const char* name = NULL, const char* title = NULL, int index = -999, int rebin = 0);
  DataVector(TVectorD);
  DataVector(TVectorD, TVectorD);
  DataVector(TVectorD, TVectorD, TVectorD);
  ~DataVector();

  void Clear();

  //set functions
  void Set(TVectorD);
  void Set(TVectorD, TVectorD);
  void Set(TVectorD, TVectorD, TVectorD);
  void Set(const char*, DataVector, int index = -999, int rebin = 0);

  //calculate the errors from the values
  void SetStatisticalError(double factor = -1.);
  void SetSystematicalError(double factor = -1.);

  void Values(TVectorD values)
  {
    fValues.ResizeTo(values);
    fValues = values;
    fStatisticalError.ResizeTo(fValues);
    fSystematicalErrorLow.ResizeTo(fValues);
    fSystematicalErrorHigh.ResizeTo(fValues);
  };

  void Value(int i, double value)
  {
    fValues[i] = value;
  }

  //the size of the errors has alreay been set => has to match

  //symmetric errors
  void StatisticalError(TVectorD error)
  {
    if(error.GetNoElements() != fStatisticalError.GetNoElements())
      {
	throw Form("%s: Error, sizes of given and existing error vector are different: %d != %d",__PRETTY_FUNCTION__,error.GetNoElements(),fStatisticalError.GetNoElements());
      }
    fStatisticalError = error;
  };

  void StatisticalError(int i, double error)
  {
    fStatisticalError[i] = error;
  }

  void SystematicalErrorLow(TVectorD errorLow)
  {
    if(errorLow.GetNoElements() != fSystematicalErrorLow.GetNoElements())
      {
	throw Form("%s: Error, sizes of given and existing error vector are different: %d != %d",__PRETTY_FUNCTION__,errorLow.GetNoElements(),fSystematicalErrorLow.GetNoElements());
      }
    fSystematicalErrorLow = errorLow;
  };

  void SystematicalErrorLow(int i, double errorLow)
  {
    fSystematicalErrorLow[i] = errorLow;
  }

  void SystematicalErrorHigh(TVectorD errorHigh)
  {
    if(errorHigh.GetNoElements() != fSystematicalErrorHigh.GetNoElements())
      {
	throw Form("%s: Error, sizes of given and existing error vector are different: %d != %d",__PRETTY_FUNCTION__,errorHigh.GetNoElements(),fSystematicalErrorHigh.GetNoElements());
      }
    fSystematicalErrorHigh = errorHigh;
  };

  void SystematicalErrorHigh(int i, double errorHigh)
  {
    fSystematicalErrorHigh[i] = errorHigh;
  }

  //symmetric errors
  void SystematicalError(TVectorD error)
  {
    if(error.GetNoElements() != fSystematicalErrorLow.GetNoElements())
      {
	throw Form("%s: Error, sizes of given and existing error vector are different: %d != %d",__PRETTY_FUNCTION__,error.GetNoElements(),fSystematicalErrorLow.GetNoElements());
      }
    fSystematicalErrorLow = error;
    fSystematicalErrorHigh = error;
  };

  void SystematicalError(int i, double error)
  {
    fSystematicalErrorLow[i] = error;
    fSystematicalErrorHigh[i] = error;
  }

  void ScaleStatisticalError(int, double);

  void Name(string name)
  {
    fName = name;
  };

  void Title(string title)
  {
    fTitle = title;
  };

  void Name(const char* name)
  {
    fName = name;
  };

  void Title(const char* title)
  {
    fTitle = title;
  };

  void Index(int index)
  {
    fIndex = index;
  };

  void Rebin(int rebin)
  {
    fRebin = rebin;
  };

  void ResizeTo(int i, const char* name = NULL, const char* title = NULL, int index = -999, int rebin = 0)
  {
    fValues.ResizeTo(i);
    fStatisticalError.ResizeTo(i);
    fSystematicalErrorLow.ResizeTo(i);
    fSystematicalErrorHigh.ResizeTo(i);
    if(name != NULL)
      {
	fName = name;
      }
    if(title != NULL)
      {
	fTitle = title;
      }
    if(index != -999)
      {
	fIndex = index;
      }
    fRebin = rebin;
  };

  //get functions
  TVectorD* Values()
  {
    return &fValues;
  };

  double& Value(int i)
  {
    return fValues[i];
  };

  TVectorD* StatisticalError()
  {
    return &fStatisticalError;
  };

  double& StatisticalError(int i)
  {
    return fStatisticalError[i];
  };

  TVectorD* SystematicalErrorLow()
  {
    return &fSystematicalErrorLow;
  };

  double& SystematicalErrorLow(int i)
  {
    return fSystematicalErrorLow[i];
  };

  TVectorD* SystematicalErrorHigh()
  {
    return &fSystematicalErrorHigh;
  };

  double& SystematicalErrorHigh(int i)
  {
    return fSystematicalErrorHigh[i];
  };

  double SystematicalError(int i)
  {
    return (fSystematicalErrorLow[i]+fSystematicalErrorHigh[i])/2.;
  };

  const char* Name() const
  {
    return fName.c_str();
  };

  const char* Title() const
  {
    return fTitle.c_str();
  };

  int Index() const
  {
    return fIndex;
  };

  int Rebin() const
  {
    return fRebin;
  };

  int Size() const
  {
    return fValues.GetNoElements();
  };

  //swap errors
  void SwapSystematicalErrors(int);

  //rebin the vector
  void Rebin(size_t, vector<int>, int VerboseLevel = 0);

  //remove all entries that are zero (including errors)
  void PurgeZeros(vector<int>);

  //index operator
  double& operator[](int index)
  {
    return fValues[index];
  };

  //basic arithmetic operators
  //DataVector
  DataVector& operator+= (const DataVector&);
  DataVector& operator-= (const DataVector&);
  DataVector& operator*= (const DataVector&);
  DataVector& operator/= (const DataVector&);

  const DataVector operator+ (const DataVector rhs)
  {
    return DataVector(*this) += rhs;
  };
  const DataVector operator- (const DataVector rhs)
  {
    return DataVector(*this) -= rhs;
  };
  const DataVector operator* (const DataVector rhs)
  {
    return DataVector(*this) *= rhs;
  };
  const DataVector operator/ (const DataVector rhs)
  {
    return DataVector(*this) /= rhs;
  };
  
  //TVectorD
  DataVector& operator+= (const TVectorD&);
  DataVector& operator-= (const TVectorD&);
  DataVector& operator*= (const TVectorD&);
  DataVector& operator/= (const TVectorD&);

  const DataVector operator+ (const TVectorD rhs)
  {
    return DataVector(*this) += rhs;
  };
  const DataVector operator- (const TVectorD rhs)
  {
    return DataVector(*this) -= rhs;
  };
  const DataVector operator* (const TVectorD rhs)
  {
    return DataVector(*this) *= rhs;
  };
  const DataVector operator/ (const TVectorD rhs)
  {
    return DataVector(*this) /= rhs;
  };

  //double  
  DataVector& operator+= (const double&);
  DataVector& operator-= (const double&);
  DataVector& operator*= (const double&);
  DataVector& operator/= (const double&);

  const DataVector operator+ (const double rhs)
  {
    return DataVector(*this) += rhs;
  };
  const DataVector operator- (const double rhs)
  {
    return DataVector(*this) -= rhs;
  };
  const DataVector operator* (const double rhs)
  {
    return DataVector(*this) *= rhs;
  };
  const DataVector operator/ (const double rhs)
  {
    return DataVector(*this) /= rhs;
  };
  
protected:
  TVectorD fValues;
  TVectorD fStatisticalError;
  TVectorD fSystematicalErrorLow;
  TVectorD fSystematicalErrorHigh;

  int fIndex;
  string fName;
  string fTitle;
  int fRebin;

  ClassDef(DataVector, 2);
};

class DataGraph : public TObject
{
public:
  DataGraph();
  DataGraph(DataVector&,DataVector&);
  DataGraph(DataVector&,TVectorD&,TVectorD&, const char* = NULL, const char* = NULL);
  ~DataGraph();

  //overloaded member functions of TGraph
  void Draw(Option_t* = NULL, const char* = NULL);

  void SetMarkerColor(Color_t, const char* = NULL);
  void SetLineColor(Color_t, const char* = NULL);
  void SetLineStyle(Style_t, const char* = NULL);
  void SetLineWidth(Width_t, const char* = NULL);

  double GetMaximum()
  {
    return fStatisticalGraph.GetMaximum();
  };

  int Fit(const char* formula, Option_t* option = "", Option_t* goption = "", Axis_t  xmin = 0, Axis_t  xmax = 0)
  {
    return fStatisticalGraph.Fit(formula, option, goption, xmin, xmax);
  };

  int Fit(TF1* f1, Option_t* option = "", Option_t* goption = "", Axis_t xmin = 0, Axis_t xmax = 0)
  {
    return fStatisticalGraph.Fit(f1, option, goption, xmin, xmax);
  };

  //these are used in the WriteGraph routine
  //they return the values of fStatisticalGraph unless other specified
  int GetN(bool systematical = false);
  int GetPoint(int, double&, double&, bool systematical = false);
  double GetErrorX(int, bool systematical = false);
  double GetErrorXlow(int, bool systematical = false);
  double GetErrorXhigh(int, bool systematical = false);
  double GetErrorY(int, bool systematical = false);
  double GetErrorYlow(int, bool systematical = false);
  double GetErrorYhigh(int, bool systematical = false);

  //set some plotting properties
  TAxis* GetXaxis(bool systematical = false);
  TAxis* GetYaxis(bool systematical = false);

  //get functions
  const char* Name();
  const char* Title();
  TGraphErrors GetStatisticalGraph()
  {
    return fStatisticalGraph;
  }
  TGraphAsymmErrors GetSystematicalGraph()
  {
    return fSystematicalGraph;
  }

  //set functions
  void Set(DataVector, DataVector);
  void SetName(const char*);
  void SetTitle(const char*);

  void Scale(double);

  void SetPoint(int, double, double, bool systematical = false);
  void SetPointError(int, double, double);
  void SetPointError(int, double, double, double, double);

  //access to y (should be the same for stat. and syst.)
  double GetY(int index)
  {
    return fStatisticalGraph.GetY()[index];
  };

  //index operator to access x
  double& operator[](int index)
  {
    return fStatisticalGraph.GetX()[index];
  };

protected:
  TGraphErrors fStatisticalGraph;
  TGraphAsymmErrors fSystematicalGraph;

  ClassDef(DataGraph, 2);
};

#endif
