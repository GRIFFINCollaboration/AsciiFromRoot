#include <iostream>

#include "TFile.h"
#include "TFolder.h"
#include "TH1F.h"
#include "TH2F.h"

#include "CommandLineInterface.hh"
#include "RootUtilities.hh"
#include "DataVector.hh"

int main(int argc, char** argv) {
  string RootFileName;
  vector<string> HistogramNames;
  string Format;
  bool Symmetrize = false;
  bool Truncate = false;
  bool Debug = false;

  CommandLineInterface* interface = new CommandLineInterface();

  interface->Add("This program creates ascii files from the provided histograms");
  interface->Add("(each in a individual file with a .dat appended to the name)");
  interface->Add("formats are:");
  interface->Add("  histograms: e = write errors (default), z = write zeros, v = only values, p = pm3d format, and r = radware format (binary, not ASCII)");
  interface->Add("      graphs: error = write errors , nozeros = don't write zeros, gnuplot = xyerrorbars format and fresco = fresco input format");
  interface->Add("The radware format for matrices gets automatically chosen based on the number of channels and resolution of the TH2 histogram:");
  interface->Add("  histograms with more than 4096 bins per axis create a (truncated) 4byte 8192x8192 file (unless the -4k flag is given)");
  interface->Add("  histograms with up to 4096 bins per axis will be written as 2byte (TH2C or TH2S) or 4byte (TH2I, TH2F, or TH2F) file");
  interface->Add("-4k", "truncate matrix to 4096x4096 entries",&Truncate);
  interface->Add("-sm", "symmetrize matrix before writing",&Symmetrize);
  interface->Add("-rf", "<name of root file>",&RootFileName);
  interface->Add("-hn", "<histogram names>",&HistogramNames);
  interface->Add("-of", "<output format (optional)>",&Format);
  interface->Add("-d",  "<turn on debugging messages>",&Debug);

  interface->CheckFlags(argc, argv);

  if(RootFileName.empty() || HistogramNames.size() == 0) {
      cerr<<"Sorry, I need the name of the root-file to read from as well as at least one histogram name"<<endl;

      return 1;
  }

  TFile* rootFile = new TFile(RootFileName.c_str());
  TObject* object;

  for(size_t i = 0; i < HistogramNames.size(); i++) {
    object = rootFile->Get(HistogramNames.at(i).c_str());

    if(Debug) {
      cout<<"got "<<i+1<<". object "<<HistogramNames.at(i)<<" from file: "<<object<<endl;
    }

    //if we couldn't find the object, loop through all object in the file, and if they're folders search in those
    if(object == NULL) {
      TIter nextkey(rootFile->GetListOfKeys());
      TKey* key;
      while((key = (TKey*)nextkey()) != NULL) {
        object = key->ReadObj();
        if(strcmp(object->ClassName(),"TFolder") == 0) {
          TFolder* folder = (TFolder*) object;
          object = folder->FindObject(folder->FindFullPathName(HistogramNames.at(i).c_str()));
          break;
        } else if(Debug) {
          cout<<"object "<<object->GetName()<<" is a "<<object->ClassName()<<", not a TFolder"<<endl;
        }
      }
      if(key == NULL || object == NULL) {
	cerr<<"couldn't find histogram "<<HistogramNames.at(i)<<" in file "<<rootFile->GetName()<<endl;

	continue;
      }
    }

    //check what type this object is: 1D- or 2D-histogram, graph w/ or w/o errors, or a data graph (custom class) 
    if(strncmp(object->ClassName(),"TH1",3) == 0)	{
      if(Debug) {
	cout<<"calling TH1 with format "<<Format<<endl;
      }

      if(Format.empty()) {
	WriteHistogram((TH1*) object,Form("%s.dat",HistogramNames.at(i).c_str()),"e");
      } else {
	if(Debug) {
	  cout<<"calling WriteHistogram((TH1*) object,Form(\"%s.dat\",HistogramNames.at(i).c_str()),Format.c_str()) = WriteHistogram((TH1*) "<<object<<", "<<Form("%s.dat",HistogramNames.at(i).c_str())<<","<<Format.c_str()<<")"<<endl;
	}
	WriteHistogram((TH1*) object,Form("%s.dat",HistogramNames.at(i).c_str()),Format.c_str());
      }

      if(Debug) {
	cout<<"done"<<endl;
      }
    } else if(strncmp(object->ClassName(),"TH2",3) == 0) {
      if(Debug) {
	cout<<"calling TH2 with format "<<Format<<endl;
      }

      if(Truncate) {
	((TH2*)object)->SetBins(4096,0.,4096.,4096,0.,4096.);
      }

      if(Symmetrize) {
	TObject* tmpObject = object->Clone("tmpHist");
	for(int i = 1; i <= ((TH2*)object)->GetNbinsX(); ++i) {
	  for(int j = 1; j <= ((TH2*)object)->GetNbinsX(); ++j) {
	    ((TH2*)object)->SetBinContent(i,j,((TH2*)object)->GetBinContent(i,j)+((TH2*)tmpObject)->GetBinContent(j,i));
	  }
	}
	delete tmpObject;
      }

      WriteHistogram((TH2*) object,Form("%s.dat",HistogramNames.at(i).c_str()),Format.c_str());

      if(Debug) {
	cout<<"done"<<endl;
      }
    } else if(strncmp(object->ClassName(),"TSpline",7) == 0) {
      if(Debug) {
	cout<<"calling TSpline"<<endl;
      }

      WriteSpline((TSpline3*) object,Form("%s.dat",HistogramNames.at(i).c_str()), ((TSpline3*) object)->GetXmin(), ((TSpline3*) object)->GetXmax(), (((TSpline3*) object)->GetXmax()-((TSpline3*) object)->GetXmin())/((TSpline3*) object)->GetNp());

      if(Debug) {
	cout<<"done"<<endl;
      }
    } else if(strcmp(object->ClassName(),"TGraph") == 0) {
      if(Debug) {
	cout<<"calling TGraph"<<endl;
      }

      if(Format.empty()) {
	WriteGraph((TGraph*) object,Form("%s.dat",HistogramNames.at(i).c_str()),"error");
      } else {
	WriteGraph((TGraph*) object,Form("%s.dat",HistogramNames.at(i).c_str()),Format.c_str());
      }

      if(Debug) {
	cout<<"done"<<endl;
      }
    } else if(strcmp(object->ClassName(),"TGraphErrors") == 0) {
      if(Debug) {
	cout<<"calling TGraphErrors"<<endl;
      }

      if(Format.empty()) {
	WriteGraph((TGraphErrors*) object,Form("%s.dat",HistogramNames.at(i).c_str()),"gnuploterror");
      } else {
	WriteGraph((TGraphErrors*) object,Form("%s.dat",HistogramNames.at(i).c_str()),Format.c_str());
      }

      if(Debug) {
	cout<<"done"<<endl;
      }
    } else if(strcmp(object->ClassName(),"TGraphAsymmErrors") == 0) {
      if(Debug) {
	cout<<"calling TGraphAsymmErrors"<<endl;
      }

      if(Format.empty()) {
	WriteGraph((TGraphAsymmErrors*) object,Form("%s.dat",HistogramNames.at(i).c_str()),"gnuploterror");
      } else {
	WriteGraph((TGraphAsymmErrors*) object,Form("%s.dat",HistogramNames.at(i).c_str()),Format.c_str());
      }

      if(Debug) {
	cout<<"done"<<endl;
      }
    } /*else if(strcmp(object->ClassName(),"DataGraph") == 0) {
      if(Debug) {
	cout<<"calling DataGraph"<<endl;
      }

      if(Format.empty()) {
	WriteDataGraph((DataGraph*) object,Form("%s.dat",HistogramNames.at(i).c_str()),"error");
      } else {
	WriteDataGraph((DataGraph*) object,Form("%s.dat",HistogramNames.at(i).c_str()),Format.c_str());
      }

      if(Debug) {
	cout<<"done"<<endl;
      }
    } */else {
      cout<<"Sorry, class "<<object->ClassName()<<" not implemented yet!"<<endl;
    }
  }//for(size_t i = 0; i < HistogramNames.size(); i++)

  rootFile->Close();

  return 0;
}
