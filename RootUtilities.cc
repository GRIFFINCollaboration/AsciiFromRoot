#include "RootUtilities.hh"

using namespace std;
using namespace TMath;

template int WriteHistogram(TH1* histogram, const char* filename=NULL, const char* format=NULL, double xLow=-1e200, double xHigh=1e200, double yLow=-1e200, double yHigh=1e200, double zLow=-1e200, double zHigh=1e200);
template int WriteHistogram(TH2* histogram, const char* filename=NULL, const char* format=NULL, double xLow=-1e200, double xHigh=1e200, double yLow=-1e200, double yHigh=1e200, double zLow=-1e200, double zHigh=1e200);
template int WriteHistogram(TH3* histogram, const char* filename=NULL, const char* format=NULL, double xLow=-1e200, double xHigh=1e200, double yLow=-1e200, double yHigh=1e200, double zLow=-1e200, double zHigh=1e200);

template <class T>
int WriteHistogram(T* histogram, const char* filename=NULL, const char* format=NULL, 
		   double xLow=-1e200, double xHigh=1e200, 
		   double yLow=-1e200, double yHigh=1e200, 
		   double zLow=-1e200, double zHigh=1e200) {
  //the number of bins is 1 for dimensions that aren't covered by histogram type
  int xBinLow, xBinHigh, xNumberOfBins;
  int yBinLow, yBinHigh, yNumberOfBins;
  int zBinLow, zBinHigh, zNumberOfBins;

  xNumberOfBins = histogram->GetNbinsX();
  yNumberOfBins = histogram->GetNbinsY();
  zNumberOfBins = histogram->GetNbinsZ();

  xBinLow = histogram->GetXaxis()->FindBin(xLow);
  yBinLow = histogram->GetYaxis()->FindBin(yLow);
  zBinLow = histogram->GetZaxis()->FindBin(zLow);

  if(xBinLow < 1)
    xBinLow = 1;
  if(yBinLow < 1)
    yBinLow = 1;
  if(zBinLow < 1)
    zBinLow = 1;
  
  xBinHigh = histogram->GetXaxis()->FindBin(xHigh);
  yBinHigh = histogram->GetYaxis()->FindBin(yHigh);
  zBinHigh = histogram->GetZaxis()->FindBin(zHigh);

  if(xBinHigh > xNumberOfBins)
    xBinHigh = xNumberOfBins;
  if(yBinHigh > yNumberOfBins)
    yBinHigh = yNumberOfBins;
  if(zBinHigh > zNumberOfBins)
    zBinHigh = zNumberOfBins;

  // writes spectrum to file or stdout if no filename is given
  ofstream output;
  streambuf* coutBuffer = cout.rdbuf();
  if(filename != NULL) {
    output.open(filename);
    if(output.rdstate() != 0) {
      cerr<<"Error, couldn't open file "<<filename<<endl;
      return 1;
    }
    
    cout.rdbuf(output.rdbuf());
  }

  //possible formats are (combinations possible):
  //p: write 2d histogram to be used as input for gnuplots pm3d mode
  //e: write errors
  //z: zero suppression
  //v: values only (no x,y or z)

  //format flags
  bool WriteErrors = false;
  bool WriteZeros = true;
  bool WriteValuesOnly = false;
  bool WritePm3d = false;
  bool WriteRadware = false;
  
  //check the format
  string formatString;
  if(format != NULL) {
      formatString = format;
      transform(formatString.begin(), formatString.end(), formatString.begin(), (int (*)(int))tolower);

      if(formatString.find("e") != string::npos) {
	WriteErrors = true;
      }
      if(formatString.find("z") != string::npos) {
	WriteZeros = false;
      }
      if(formatString.find("v") != string::npos) {
	WriteValuesOnly = true;
      }
      if(formatString.find("p") != string::npos) {
	WritePm3d = true;
      }
      if(formatString.find("r") != string::npos) {
	WriteRadware = true;
      }

      if(WritePm3d && (WriteErrors || WriteValuesOnly)) {
	cerr<<"Error, can't write pm3d format and errors or only values!"<<endl;
	return 1;
      }
      if(WritePm3d && (yBinHigh <= yBinLow || zBinHigh > zBinLow)) {
	cerr<<"Error, pm3d format requested but histogram doesn't seem to be a TH2: yBinLow = "<<yBinLow<<", yBinHigh = "<<yBinHigh<<", zBinLow = "<<zBinLow<<", zBinHigh = "<<zBinHigh<<", class name = "<<histogram->IsA()->GetName()<<endl;
	return 2;
      }
      if(filename == NULL && WriteRadware) {
	cerr<<"Error, won't write binary radware format to stdout, need to have a file name!"<<endl;
	return 3;
      }
  }

  if(!WritePm3d && !WriteRadware) {
    //default: write x(,y,z) and counts
    //if no y (or z) axis exist the respective loop will execute only once (BinLow == BinHigh == 1)
    //only write the y (or z) if BinHigh > BinLow
    for(int i = xBinLow; i <= xBinHigh; i++) {
      for(int j = yBinLow; j <= yBinHigh; j++) {
	for(int k = zBinLow; k <= zBinHigh; k++) {
	  if(WriteZeros || histogram->GetBinContent(i,j,k) > 0) {
	    if(!WriteValuesOnly) {
	      cout<<setw(8)<<histogram->GetXaxis()->GetBinCenter(i);
	      if(yBinHigh > yBinLow)
		cout<<" "<<setw(8)<<histogram->GetYaxis()->GetBinCenter(j);
	      if(zBinHigh > zBinLow)
		cout<<" "<<setw(8)<<histogram->GetZaxis()->GetBinCenter(k);
	    }
	    cout<<" "<<setw(8)<<histogram->GetBinContent(i,j,k);
	    if(WriteErrors)
	      cout<<" "<<setw(8)<<histogram->GetBinError(i,j,k);
	    cout<<endl;
	  }
	}
      }
    }
  } else if(WritePm3d) {
    // special formatting for gnuplots pm3d (set pm3d map ; splot 'filename')
    // writes spectrum to file or stdout if no filename is given
    // all counts are written 
    // blank lines between different x
    // if all values are zero for multiple consecutive x-values only the first and last are written
    // if values of consecutive y are all the same only the first is printed
    // after each x-set the same set is printed again at the next x-value

    //we need to get the lower edge of the next bins as well
    xBinHigh++;
    yBinHigh++;

    vector<vector<double> > values(xBinHigh-xBinLow+1);
    vector<double> x(xBinHigh-xBinLow+1);
    vector<vector<double> > y(xBinHigh-xBinLow+1);
      
    size_t erasedX = 0;
    vector<size_t> erasedY(xBinHigh-xBinLow+1,0);
      
    //get the values
    for(int i = xBinLow; i <= xBinHigh; i++) {
      x.at(i-xBinLow-erasedX) = histogram->GetXaxis()->GetBinLowEdge(i);

      y.at(i-xBinLow-erasedX).resize(yBinHigh-yBinLow+1);
      values.at(i-xBinLow-erasedX).resize(yBinHigh-yBinLow+1);
	  
      //read the first y-value in any case
      y.at(i-xBinLow-erasedX).at(0) = histogram->GetYaxis()->GetBinLowEdge(yBinLow);
      values.at(i-xBinLow-erasedX).at(0) = histogram->GetBinContent(i,yBinLow);

      for(int j = yBinLow+1; j <= yBinHigh; j++) {
	//check whether this value is different than the last (or the last value)
	if(histogram->GetBinContent(i,j) != values.at(i-xBinLow-erasedX).at(j-yBinLow-erasedY.at(i-xBinLow-erasedX)-1) || j == yBinHigh) {
	  y.at(i-xBinLow-erasedX).at(j-yBinLow-erasedY.at(i-xBinLow-erasedX)) = histogram->GetYaxis()->GetBinLowEdge(j);
	  values.at(i-xBinLow-erasedX).at(j-yBinLow-erasedY.at(i-xBinLow-erasedX)) = histogram->GetBinContent(i,j);
	} else {
	  //erase this y and increase yBinLow
	  values.at(i-xBinLow-erasedX).erase(values.at(i-xBinLow-erasedX).begin()+j-yBinLow-erasedY.at(i-xBinLow-erasedX));
	  y.at(i-xBinLow-erasedX).erase(y.at(i-xBinLow-erasedX).begin()+j-yBinLow-erasedY.at(i-xBinLow-erasedX));
	  erasedY.at(i-xBinLow-erasedX)++;
	}
      }
      //check if this x and the last had the same values (and are not the first x)
      //=> if so erase this x (i) and increase erasedX
      if(xBinLow < i && i < xBinHigh) {
	if(values.at(i-xBinLow-erasedX) == values.at(i-xBinLow-erasedX-1)) {
	  values.erase(values.begin()+i-xBinLow-erasedX);
	  y.erase(y.begin()+i-xBinLow-erasedX);
	  x.erase(x.begin()+i-xBinLow-erasedX);
	  erasedY.erase(erasedY.begin()+i-xBinLow-erasedX);
	  erasedX++;
	}
      } else {//if(i > xBinLow)
	cout<<"not erasing values since i = "<<i<<" not between "<<xBinLow<<" and "<<xBinHigh-erasedX<<endl;
      }
    }//for(int i = xBinLow; i <= xBinHigh; i++)
      
    //printing of values (and instructions)
    cout<<"#set pm3d map corners2color c1"<<endl
	<<"#set palette rgbformulae 22,13,-31"<<endl
	<<"#splot 'filename'"<<endl;
    for(size_t i = 0; i < values.size(); i++) {
      for(size_t j = 0; j < values.at(i).size(); j++) {
	cout<<setw(8)<<x.at(i)<<setw(8)<<y.at(i).at(j)<<setw(8)<<values.at(i).at(j)<<endl;
      }
      cout<<endl;
      if(i < values.size()-1) {
	for(size_t j = 0; j < values.at(i).size(); j++) {
	  cout<<setw(8)<<x.at(i+1)<<setw(8)<<y.at(i).at(j)<<setw(8)<<values.at(i).at(j)<<endl;
	}
	cout<<endl;
      } else {
	for(size_t j = 0; j < values.at(i).size(); j++) {
	  cout<<setw(8)<<x.at(i)+histogram->GetXaxis()->GetBinWidth(1)<<setw(8)<<y.at(i).at(j)<<setw(8)<<values.at(i).at(j)<<endl;
	}
	cout<<endl;
      }
    }
  } else if(WriteRadware) {//else if(WritePm3d)
    //we're going to write binary data => need to write to file (and re-open it as binary)
    cout.rdbuf(coutBuffer);
    if(output.is_open()) {
      output.close();
    }
    output.open(filename, ofstream::out | ofstream::binary);
    if(output.rdstate() != 0) {
      cerr<<"Error, couldn't open file "<<filename<<endl;
      return 1;
    }

    //if the high bin of the y-axis is less than it's low bin, it's a 1D-histogram
    if(yBinHigh <= yBinLow) {
      //spe header format (value - #bytes)
      //24 - 4, name - 8, #elements - 4, 3 x (1 - 4), 24 - 4, 4*#elements - 4
      uint32_t header[9];
      header[0] = 24;
      char name[8];
      strncpy(name,histogram->GetName(),sizeof(name));
      header[1] = ((uint32_t)name[0])<<24 | ((uint32_t)name[1])<<16 | ((uint32_t)name[2])<<8 | name[3];
      header[2] = ((uint32_t)name[4])<<24 | ((uint32_t)name[5])<<16 | ((uint32_t)name[6])<<8 | name[7];
      header[3] = histogram->GetNbinsX();
      header[4] = 1;
      header[5] = 1;
      header[6] = 1;
      header[7] = 24;
      header[8] = 4*header[3];
      output.write((const char*) header, sizeof(header));

      float* values = new float[header[3]];
      for(uint32_t i = 1; i <= header[3]; ++i) {
	values[i-1] = histogram->GetBinContent(i);
      }
      output.write((const char*) values, header[3]*sizeof(values[0]));
      output.write((const char*) header+8, sizeof(header[8]));//footer
      //cerr<<"nof bins = "<<histogram->GetNbinsX()<<", header[3] = "<<header[3]<<", sizeof(header) = "<<sizeof(header)<<", sizeof(values[0]) = "<<sizeof(values[0])<<endl;
      delete[] values;
    } else if(zBinHigh <= zBinLow) {
      //.mat = matrix of 2-byte words 4096*4096
      //.spn or .m4b = matrix of 4-byte words 4096*4096
      //.sec = large matrix of 4-byte words 8192*8192
      //first we check the dimensions, if any are larger that 4096, we use the large 4-byte matrix
      if(histogram->GetNbinsX() > 4096 || histogram->GetNbinsY() > 4096) {
	uint32_t* values = new uint32_t[(0x1)<<26];
	//cerr<<"writing large matrix, sizeof(values) = "<<sizeof(values)<<" = 4*"<<((0x1)<<26)<<endl;
	for(int i = 0; i < ((0x1)<<26); ++i) {
	  values[i] = histogram->GetBinContent((i%8192)+1,(i/8192)+1);
	}
	output.write((const char*) values, (0x1)<<28);//4 bytes per value
      } else {
	if(strcmp(histogram->IsA()->GetName(),"TH2C") == 0 || strcmp(histogram->IsA()->GetName(),"TH2S") == 0) {
	  //these two classes use characters and shorts, i.e. we can fit them in a 2-byte matrix
	  uint16_t* values = new uint16_t[(0x1)<<24];
	  for(int i = 0; i < ((0x1)<<24); ++i) {
	    values[i] = histogram->GetBinContent((i%4096)+1,(i/4096)+1);
	  }
	  output.write((const char*) values, (0x1)<<25);//2 bytes per value
	} else {
	  //these other classes use integers, floats, or doubles, i.e. we need to use a 4-byte matrix
	  uint32_t* values = new uint32_t[(0x1)<<24];
	  for(int i = 0; i < ((0x1)<<24); ++i) {
	    values[i] = histogram->GetBinContent((i%4096)+1,(i/4096)+1);
	  }
	  output.write((const char*) values, (0x1)<<26);//4 bytes per value
	}
      }
    } else {
      //.gcub = cube of 4-byte words 1024*1024*1024
      uint32_t* values = new uint32_t[(0x1)<<30];
      //cerr<<"writing cube, sizeof(values) = "<<sizeof(values)<<" = 4*"<<((0x1)<<30)<<endl;
      for(uint32_t i = 0; i < ((0x1)<<30); ++i) {
	values[i] = histogram->GetBinContent((i%(1024*1024))+1,((i/1024)%1024)+1,(i/1024/1024)+1);
      }
      //1<<32 is larger than the (internal) 32bit, so we use streamsize
      streamsize tmp = 0x1;
      output.write((const char*) values, tmp<<32);//4 bytes per value
    }
    output.close();
  }

  cout.rdbuf(coutBuffer);

  if(output.is_open()) {
    output.close();
  }

  return 0;
}

int WriteSpline(TSpline3* s, const char* filename, double xlow, double xhigh, double xstep) 
{
  // writes spectrum to file or stdout if no filename is given
  ofstream output;
  streambuf* coutBuffer = cout.rdbuf();
  if(filename != NULL)
    {
      output.open(filename);
      if(output.rdstate() != 0) {
	  cerr<<"Error, couldn't open file "<<filename<<endl;
	  return 1;
	}

      cout.rdbuf(output.rdbuf());
    }

  for(double x = xlow; x <= xhigh; x += xstep) 
    {
      cout<<x<<" "<<s->Eval(x)<<endl;;
    }

  cout.rdbuf(coutBuffer);

  output.close();

  return 0;
}

int WriteFunction(TF1* function, const char* filename, int numberOfSteps, double xmin, double xmax)
{
  // writes spectrum to file or stdout if no filename is given
  ofstream output;
  streambuf* coutBuffer = cout.rdbuf();
  if(filename != NULL)
    {
      output.open(filename);
      if(output.rdstate() != 0) {
	  cerr<<"Error, couldn't open file "<<filename<<endl;
	  return 1;
	}

      cout.rdbuf(output.rdbuf());
    }

  if(xmin > xmax)
    {
      function->GetRange(xmin, xmax);
    }

  for(double x = xmin; x <= xmax; x += (xmax - xmin)/numberOfSteps)
    {
      cout<<setw(8)<<x<<" "<<setw(8)<<function->Eval(x)<<endl;
    }

  cout.rdbuf(coutBuffer);

  output.close();

  return 0;
}

int WriteFunction(TF2* function, const char* filename, int numberOfSteps)
{
  // writes spectrum to file or stdout if no filename is given
  ofstream output;
  streambuf* coutBuffer = cout.rdbuf();
  if(filename != NULL)
    {
      output.open(filename);
      if(output.rdstate() != 0) {
	  cerr<<"Error, couldn't open file "<<filename<<endl;
	  return 1;
	}

      cout.rdbuf(output.rdbuf());
    }

  double xmin,xmax;
  double ymin,ymax;

  function->GetRange(xmin, ymin, xmax, ymax);

  for(double x = xmin; x <= xmax; x += (xmax - xmin)/numberOfSteps)
    {
      for(double y = ymin; y <= ymax; y += (ymax - ymin)/numberOfSteps) {
	  cout<<setw(8)<<x<<" "<<setw(8)<<y<<" "<<setw(8)<<function->Eval(x,y)<<endl;
	}
    }

  cout.rdbuf(coutBuffer);

  output.close();

  return 0;
}

/*int WriteDataGraph(DataGraph* graph, const char* filename, const char* format, const char* header)
{
  // writes spectrum to file or stdout if no filename is given
  ofstream output;
  streambuf* coutBuffer = cout.rdbuf();
  if(filename != NULL)
    {
      output.open(filename);
      if(output.rdstate() != 0) {
	  cerr<<"Error, couldn't open file "<<filename<<endl;
	  return 1;
	}

      cout.rdbuf(output.rdbuf());
    }

  //possible formats are (combinations possible):
  //gnuplot: x,y,xlow,xhigh,ylow,yhigh, and lows and highs are the absolute values
  //fresco: write x,y,errror
  //error: write errors
  //nozero: zero suppression
  //comb: combine the stat. and syst. errors

  //format flags
  bool WriteErrors = false;
  bool WriteZeros = true;
  bool Gnuplot = false;
  bool Fresco = false;
  bool CombineErrors = false;

  //check the format
  string formatString;
  if(format != NULL)
    {
      formatString = format;
      transform(formatString.begin(), formatString.end(), formatString.begin(), (int (*)(int))tolower);

      if(formatString.find("error") != string::npos) {
	  WriteErrors = true;
	}
      if(formatString.find("nozero") != string::npos) {
	  WriteZeros = false;
	}
      if(formatString.find("gnuplot") != string::npos) {
	  Gnuplot = true;
	}
      if(formatString.find("fresco") != string::npos) {
	  Fresco = true;
	}
      if(formatString.find("comb") != string::npos) {
	  CombineErrors = true;
	}

      if(Gnuplot && Fresco) {
	  cerr<<"Error, can't write gnuplot and fresco format at once!"<<endl;

	  return 1;
	}
    }

  double x, y;
  
  if(header != NULL)
    {
      cout<<header<<endl;
    }
  
  for(int point = 0; point < graph->GetN(); point++)
    {
      graph->GetPoint(point, x, y);
      if(y != 0 || WriteZeros) {
	  if(!WriteErrors && !Gnuplot && !Fresco) {
	      cout<<setw(8)<<x<<" "<<setw(8)<<y<<endl;
	    }
	  else if(Gnuplot) {
	      //x,y,xlow_stat,xhigh_stat,ylow_stat,yhigh_stat,xlow_sys,xhigh_sys,ylow_sys,yhigh_sys, and lows and highs are the absolute values
	      if(CombineErrors) {
		  cout<<setw(8)<<x<<" "<<setw(8)<<y<<" "
		      <<setw(8)<<x-Sqrt(Power(graph->GetErrorXlow(point),2)+Power(graph->GetErrorXlow(point, true),2))<<" "
		      <<setw(8)<<x+Sqrt(Power(graph->GetErrorXhigh(point),2)+Power(graph->GetErrorXhigh(point, true),2))<<" "
		      <<setw(8)<<y-Sqrt(Power(graph->GetErrorYlow(point),2)+Power(graph->GetErrorYlow(point, true),2))<<" "
		      <<setw(8)<<y+Sqrt(Power(graph->GetErrorYhigh(point),2)+Power(graph->GetErrorYhigh(point, true),2))<<" "
		      <<endl;
		}
	      else
		{
		  cout<<setw(8)<<x<<" "<<setw(8)<<y<<" "
		      <<setw(8)<<x-graph->GetErrorXlow(point)<<" "<<setw(8)<<x+graph->GetErrorXhigh(point)<<" "
		      <<setw(8)<<y-graph->GetErrorYlow(point)<<" "<<setw(8)<<y+graph->GetErrorYhigh(point)<<" "
		      <<setw(8)<<x-graph->GetErrorXlow(point, true)<<" "<<setw(8)<<x+graph->GetErrorXhigh(point, true)<<" "
		      <<setw(8)<<y-graph->GetErrorYlow(point, true)<<" "<<setw(8)<<y+graph->GetErrorYhigh(point, true)<<endl;
		}
	    }
	  else if(Fresco) {
	      //x,y,yerror
	      if(CombineErrors) {
		  cout<<setw(8)<<x<<" "<<setw(8)<<y<<" "<<setw(8)<<Sqrt(Power(Max(graph->GetErrorYhigh(point),      graph->GetErrorYlow(point)),2)+
									Power(Max(graph->GetErrorYhigh(point, true),graph->GetErrorYlow(point, true)),2))<<endl;
		}
	      else
		{
		  cout<<setw(8)<<x<<" "<<setw(8)<<y<<" "<<setw(8)<<graph->GetErrorY(point)<<endl;
		}
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
}*/

bool CreateGraph(TGraphErrors*& graph, vector<double> x, vector<double> y, vector<double> ex, vector<double> ey)
{
  //don't check the length of the vectors, just create a graph based on the x-vector (but check this one)
  if(x.size() == 0)
    {
      cerr<<"Error, can't create a graph with zero x-entries"<<endl;
      
      return false;
    }
  double* xArray = new double[x.size()];
  double* exArray = new double[x.size()];
  double* yArray = new double[x.size()];
  double* eyArray = new double[x.size()];

  for(size_t i = 0; i < x.size(); i++)
    {
      try
	{
	  xArray[i] = x.at(i);
	}
      catch(std::out_of_range& error) {
	  cout<<"Caught std::out_of_range for x.at("<<i<<"/"<<x.size()<<"): "<<error.what()<<endl;
	  return false;
	}//catch
      try
	{
	  if(ex.size() != 0) {
	      exArray[i] = ex.at(i);
	    }
	  else
	    {
	      exArray[i] = 0.;
	    }
	}
      catch(std::out_of_range& error) 
	{
	  cout<<"Caught std::out_of_range for ex.at("<<i<<"/"<<ex.size()<<"): "<<error.what()<<endl;
	  return false;
	}//catch
      try
	{
	  yArray[i] = y.at(i);
	}
      catch(std::out_of_range& error) 
	{
	  cout<<"Caught std::out_of_range for y.at("<<i<<"/"<<y.size()<<"): "<<error.what()<<endl;
	  return false;
	}//catch
      try
	{
	  if(ey.size() != 0) {
	      eyArray[i] = ey.at(i);
	    }
	  else
	    {
	      eyArray[i] = 0.;
	    }
	}
      catch(std::out_of_range& error) 
	{
	  cout<<"Caught std::out_of_range for ey.at("<<i<<"/"<<ey.size()<<"): "<<error.what()<<endl;
	  return false;
	}//catch
    }

  graph = new TGraphErrors(x.size(),xArray,yArray,exArray,eyArray);

  if(graph == NULL)
    {
      return false;
    }

  return true;
}

int get_peak_area(double sigma, double dsigma, double height, double dheight, double keVperbin) {

  /* calculate intensity=area of peak */
  double area = (sqrt(2*(TMath::Pi())) * sigma * height ) / (double)keVperbin;

  /* calculate error */ 
  double darea = ( sqrt(2*(TMath::Pi())) / (double)keVperbin) * sqrt( pow(height*dsigma,2) + pow(sigma*dheight,2) ); 

  /* output: area, darea */
  printf("\narea: %f; darea: %f\n\n", area, darea);

  return 0;
}


int get_normed_intensity(double area, double darea, double in, double din) {

  /* calculate normed intensity */
  double normin = area/in;
  
  /* calculate error of normed intensity */
  double dnormin = sqrt( pow(darea/in,2) + pow((area*din)/(in*in),2) );

  /* output: normin, dnormin */
  printf("\nnormed int.: %f; d(normed int.): %f\n\n", normin, dnormin);

  return 0;
}


int get_counts_above_threshold(TH1F* h, double xlow, double xhigh, double thresh) {

  // determines counts above certain threshold in a given range 
  /* get bin_low and bin_high */
  int blow, bhigh, len;
  double counts=0;
  len=h->GetNbinsX();
  blow=h->FindBin(xlow);
  if (blow<1) blow=1;
  
  bhigh=h->FindBin(xhigh);
  if (bhigh>len) bhigh=len;


  /* loop over bins */
  for (int i=blow; i<=bhigh; i++) {
    //fprintf(fp,"%e %e %e\n",h->GetBinCenter(i),h->GetBinContent(i),h->GetBinError(i));
    /* check if value above threshold */
    if(h->GetBinContent(i) >= thresh) {
      counts += h->GetBinContent(i) - thresh;
      printf("BinContent(%d): %f, thresh: %f -> diff: %f -> counts: %f\n", i, h->GetBinContent(i), thresh, h->GetBinContent(i) - thresh, counts);
    }
  }
  
  /* print result to stdout */
  printf("\nHistogram: %s\nCounts above %f in range: %f <-> %f: %f\n\n", h->GetName(), thresh, xlow, xhigh, counts);
  
  return 0;
}


int tranform_theta_phi_from_hardware_to_MB_frame(double theta, double phi) {

  /* copy&paste from '/home/cb/oniederm/offl_ana/is_410_oct_03/geometry/geometry.cc:get_positions_segments()' */

  /* angles in rad */
  phi   = phi/180.0*(TMath::Pi());
  theta = theta/180.0*(TMath::Pi());
  
  /* ********* */
  /* IMPORTANT */
  /* ********* */
  /* update theta and phi to be angles in the standard MINIBALL frame of reference */
  /* theta and phi passed to this routine are values read from the MINIBALL frame (HARDWARE), which
   * correspond to value in a system with the z-axis pointing upward, i.e. the x-direction of the
   * MINIBALL frame of reference (FoR).  The following lines transform the theta and phi values. */

  /* check if phi > PI -> rotation angles change */
  if(phi>(TMath::Pi())) {
    phi = 2*(TMath::Pi()) - phi;
    theta = -theta;
  }

  TVector3 trans;
  trans.SetXYZ(1,0,0);
  trans.RotateY(-phi);
  trans.RotateX(-theta);

  theta=trans.Theta();
  phi=trans.Phi();

  /* done transforming theta and phi */

  /* END of copy&paste */
 
  /* angles in deg */
  theta = theta*180.0/(TMath::Pi());
  if(phi < 0) {
    phi += 2*(TMath::Pi()); 
  }
  phi = phi*180.0/(TMath::Pi());

  /* output of result */
  printf("\nMINIBALL frame: Theta = %f, Phi = %f\n\n", theta, phi);

  return 0;
}


int tranform_theta_phi_alpha_from_hardware_to_MB_frame(double theta, double phi, double alpha) {

  /* copy& paste from 'geometry.cc' */
  /* angles in rad */
  alpha = alpha/180.0*(TMath::Pi());  
  phi   = phi/180.0*(TMath::Pi());
  theta = theta/180.0*(TMath::Pi());
  
  /* ********* */
  /* IMPORTANT */
  /* ********* */
  /* update theta and phi to be angles in the standard MINIBALL frame of reference */
  /* theta and phi passed to this routine are values read from the MINIBALL frame (HARDWARE), which
   * correspond to value in a system with the z-axis pointing upward, i.e. the x-direction of the
   * MINIBALL frame of reference (FoR).  The following lines transform the theta and phi values. */

  /* check if phi > PI -> rotation angles change */
  if(phi>(TMath::Pi())) {
    phi = 2*(TMath::Pi()) - phi;
    theta = -theta;
  }

  TVector3 trans;
  trans.SetXYZ(1,0,0);
  trans.RotateY(-phi);
  trans.RotateX(-theta);

  theta=trans.Theta();
  phi=trans.Phi();


  /* done transforming theta and phi */

  /* find alpha offset */
  
  TVector3 vd,va,vn;
  const TVector3 ex(1,0,0);
  
  vd.SetXYZ(0,0,1);
  va.SetXYZ(1,0,0);
  
  vd.RotateY(theta);
  vd.RotateZ(phi);
  va.RotateY(theta);
  va.RotateZ(phi);
  
  vn=ex-vd*(ex*vd);
  vn.SetMag(1.);

  double spat=vn*(va.Cross(vd));
  double alpha_offset=va.Angle(vn);  //*spat;
  if (spat<0) alpha_offset*=-1;

  /* END of copy&paste */
 
  /* angles in deg */
  theta = theta*180.0/(TMath::Pi());
  if(phi < 0) {
    phi += 2*(TMath::Pi()); 
  }
  phi = phi*180.0/(TMath::Pi());
  alpha = alpha*180.0/(TMath::Pi());

  /* output */
  printf("\nMINIBALL frame: Theta = %f, Phi = %f, Alpha-Rotation-Angle = %f\n\n", theta, phi, -alpha-alpha_offset);

#if 0         
  /* rotation of all segments */
  for(i=0;i<24;i++) {

    v2[i].SetXYZ(x[i],y[i],z[i]);
    v2[i].RotateZ(-alpha-alpha_offset);
    v2[i].RotateY(theta);
    v2[i].RotateZ(phi); 
  }
#endif

  return 0;
}

TGraph* plot_spline(char *splinedata, int datapoints)
{
  /* variables needed */
  double* x = new double[datapoints];
  double* y = new double[datapoints];
  ifstream varfile;
  //int i;
  int datactr = 0;

  /* open input file */
  varfile.open(splinedata, ios::in);
  if (!varfile.is_open())
    {
      printf("%s: ERROR: couldn't open file '%s'\n", __func__, splinedata);

      delete x;
      delete y;

      return(NULL);
    }
  
  /* read table */
  while (!varfile.eof() && datactr<datapoints) 
    {
      varfile>>x[datactr]>>y[datactr];
      datactr++;
    }
  datactr--;

  cout<<"read "<<datactr<<" data points"<<endl;
  //for(i=0;i<datactr;i++)
  //cout<<x[i]<<" "<<y[i]<<endl;

  /* create TSpline3 object */
  TGraph*   graph  = new TGraph(datactr, x, y);

  graph->SetLineWidth(4);

  /* close input file */
  varfile.close();

  delete x;
  delete y;

  return(graph);
}
