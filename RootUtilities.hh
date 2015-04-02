////////////////////////////////////////////////////////////////////////////////
// provide basic functions to write special root objects to ascii-file and
// other functions
// this implementation needs c++11 (i.e. gcc 4.8 or newer)
////////////////////////////////////////////////////////////////////////////////

#ifndef __UTIL_HH
#define __UTIL_HH

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cctype>
#include <stdexcept>

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TClassTable.h"
#include "TSystem.h"
#include "TCollection.h"
#include "TKey.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TF2.h"
#include "TText.h"
#include "TMath.h"
#include "TVector3.h"
#include "TRotation.h"
#include "TSpline.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

#include "DataVector.hh"

using namespace std;
using namespace TMath;

template <class T> int WriteHistogram(T* histogram, const char* filename=NULL, const char* format=NULL, double xLow=-1e200, double xHigh=1e200, double yLow=-1e200, double yHigh=1e200, double zLow=-1e200, double zHigh=1e200);
extern template int WriteHistogram(TH1* histogram, const char* filename=NULL, const char* format=NULL, double xLow=-1e200, double xHigh=1e200, double yLow=-1e200, double yHigh=1e200, double zLow=-1e200, double zHigh=1e200);
extern template int WriteHistogram(TH2* histogram, const char* filename=NULL, const char* format=NULL, double xLow=-1e200, double xHigh=1e200, double yLow=-1e200, double yHigh=1e200, double zLow=-1e200, double zHigh=1e200);
extern template int WriteHistogram(TH3* histogram, const char* filename=NULL, const char* format=NULL, double xLow=-1e200, double xHigh=1e200, double yLow=-1e200, double yHigh=1e200, double zLow=-1e200, double zHigh=1e200);
//extern template int WriteHistogram<TH1C*>(TH1C*, const char* filename=NULL, const char* format=NULL, double xLow=-1e200, double xHigh=1e200, double yLow=-1e200, double yHigh=1e200, double zLow=-1e200, double zHigh=1e200);
//extern template int WriteHistogram<TH1S*>(TH1S*, const char* filename=NULL, const char* format=NULL, double xLow=-1e200, double xHigh=1e200, double yLow=-1e200, double yHigh=1e200, double zLow=-1e200, double zHigh=1e200);
//extern template int WriteHistogram<TH1I*>(TH1I*, const char* filename=NULL, const char* format=NULL, double xLow=-1e200, double xHigh=1e200, double yLow=-1e200, double yHigh=1e200, double zLow=-1e200, double zHigh=1e200);
//extern template int WriteHistogram<TH1F*>(TH1F*, const char* filename=NULL, const char* format=NULL, double xLow=-1e200, double xHigh=1e200, double yLow=-1e200, double yHigh=1e200, double zLow=-1e200, double zHigh=1e200);
//extern template int WriteHistogram<TH1D*>(TH1D*, const char* filename=NULL, const char* format=NULL, double xLow=-1e200, double xHigh=1e200, double yLow=-1e200, double yHigh=1e200, double zLow=-1e200, double zHigh=1e200);

int WriteSpline(TSpline3* s, const char* filename, double xlow, double xhigh, double xstep);
int WriteFunction(TF1* function, const char* filename, int numberOfSteps, double xmin = 1e9, double xmax = -1e9);
int WriteFunction(TF2* function, const char* filename, int numberOfSteps);
int WriteDataGraph(DataGraph* graph, const char* filename, const char* format, const char* header = NULL);

int get_peak_area(double sigma, double dsigma, double height, double dheight, double keVperbin);
int get_normed_intensity(double area, double darea, double in, double din);
int get_counts_above_threshold(TH1F* h, double xlow=-1e200, double xhigh=1e200, double thresh=0);
int tranform_theta_phi_from_hardware_to_MB_frame(double theta, double phi);
int tranform_theta_phi_alpha_from_hardware_to_MB_frame(double theta, double phi, double alpha);
TGraph* plot_spline(char* splinedata, int datapoints);

bool CreateGraph(TGraphErrors*&, vector<double>, vector<double>, vector<double> ex = vector<double>(0), vector<double> ey = vector<double>(0));

template <class T>
int WriteGraph(T* graph, const char* filename, const char* format)
{
  // writes spectrum to file or stdout if no filename is given
  ofstream output;
  streambuf* coutBuffer = cout.rdbuf();
  if(filename != NULL)
    {
      output.open(filename);
      if(output.rdstate() != 0)
	{
	  cerr<<"Error, couldn't open file "<<filename<<endl;
	  return 1;
	}

      cout.rdbuf(output.rdbuf());
    }

  //possible formats are (combinations possible):
  //gnuplot: x,y,xlow,xhigh,ylow,yhigh, and lows and highs are the absolute values (if TGraphAsymmErrors)
  //      or x,y,delta x,delta y (if TGraphErrors)
  //fresco: write x,y,errror
  //error: write errors
  //nozero: zero suppression

  //format flags
  bool WriteErrors = false;
  bool WriteZeros = true;
  bool Gnuplot = false;
  bool Fresco = false;

  //check the format
  string formatString;
  if(format != NULL)
    {
      formatString = format;
      transform(formatString.begin(), formatString.end(), formatString.begin(), (int (*)(int))tolower);

      if(formatString.find("error") != string::npos)
	{
	  WriteErrors = true;
	}
      if(formatString.find("nozero") != string::npos)
	{
	  WriteZeros = false;
	}
      if(formatString.find("gnuplot") != string::npos)
	{
	  Gnuplot = true;
	}
      if(formatString.find("fresco") != string::npos)
	{
	  Fresco = true;
	}

      if(Gnuplot && Fresco)
	{
	  cerr<<"Error, can't write gnuplot and fresco format at once!"<<endl;

	  return 1;
	}
    }

  double x, y;
  
  for(int point = 0; point < graph->GetN(); point++)
    {
      graph->GetPoint(point, x, y);
      if(y != 0 || WriteZeros)
	{
	  if(!WriteErrors && !Gnuplot && !Fresco)
	    {
	      cout<<setw(8)<<x<<" "<<setw(8)<<y<<endl;
	    }
	  else if(Gnuplot)
	    {
	      if(strcmp(graph->IsA()->GetName(),"TGraphAsymmErrors") == 0)
		{
		  //x,y,xlow,xhigh,ylow,yhigh, and lows and highs are the absolute values
		  cout<<setw(8)<<x<<" "<<setw(8)<<y<<" "<<setw(8)<<x-graph->GetErrorXlow(point)<<" "<<setw(8)<<x+graph->GetErrorXhigh(point)<<" "<<setw(8)<<y-graph->GetErrorYlow(point)<<" "<<setw(8)<<y+graph->GetErrorYhigh(point)<<endl;
		}
	      else
		{
		  //x,y,delta x,delta y
		  cout<<setw(8)<<x<<" "<<setw(8)<<y<<" "<<setw(8)<<graph->GetErrorX(point)<<" "<<setw(8)<<graph->GetErrorY(point)<<endl;
		}
	    }
	  else if(Fresco)
	    {
	      //x,y,yerror
	      cout<<setw(8)<<x<<" "<<setw(8)<<y<<" "<<setw(8)<<graph->GetErrorY(point)<<endl;
	    }
	  else
	    {
	      cout<<setw(8)<<x<<" "<<setw(8)<<graph->GetErrorXlow(point)<<" "<<setw(8)<<graph->GetErrorXhigh(point)<<" "<<setw(8)<<y<<" "<<setw(8)<<graph->GetErrorYlow(point)<<" "<<setw(8)<<graph->GetErrorYhigh(point)<<endl;
	    }
	}
    }

  cout.rdbuf(coutBuffer);

  output.close();

  return 0;
}

template <class T>
int ReadHistogram(T* histogram, const char* filename)
{
  if(histogram == NULL || filename == NULL)
    {
      cerr<<"Either histogram or filename were NULL: "<<histogram<<", "<<filename<<endl;

      return 1;
    }

  ifstream input(filename);

  string line;
  istringstream Stream;

  int Dimension = histogram->GetDimension();
  vector<double> x(Dimension);

  while(getline(input, line).good())
    {
      if(line[0] == '#')
	{
	  continue;
	}
      Stream.clear();
      Stream.str(line);

      for(int i = 0; i < Dimension; i++)
	{		  
	  Stream>>x[i];
	}

      switch(Dimension)
	{
	case 1:
	  histogram->Fill(x[0]);
	  break;
	case 2:
	  histogram->Fill(x[0],x[1]);
	  break;
//	case 3:
//	  histogram->Fill(x[0],x[1],x[2]);
//	  break;
	default:
	  break;
	}
    }

  return 0;
}


#endif
