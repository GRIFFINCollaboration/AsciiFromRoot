#include "DataVector.hh"

using namespace std;
using namespace TMath;

ClassImp(DataVector)
ClassImp(DataGraph)

DataVector::DataVector()
{
  fValues.Clear();
  fStatisticalError.Clear();
  fSystematicalErrorLow.Clear();
  fSystematicalErrorHigh.Clear();

  fIndex = -999;
  fName.clear();
  fTitle.clear();
  fRebin = 0;
}

DataVector::DataVector(const char* name, DataVector vector, int index, int rebin)
{
  fValues.ResizeTo(vector.Size());
  fStatisticalError.ResizeTo(vector.Size());
  fSystematicalErrorLow.ResizeTo(vector.Size());
  fSystematicalErrorHigh.ResizeTo(vector.Size());
  fValues = *(vector.Values());
  fStatisticalError  = *(vector.StatisticalError());
  fSystematicalErrorLow  = *(vector.SystematicalErrorLow());
  fSystematicalErrorHigh = *(vector.SystematicalErrorHigh());

  fName = name;
  fTitle = vector.Title();
  if(index != -999)
    {
      fIndex = index;
    }
  else
    {
      fIndex = vector.Index();
    }

  if(rebin != 0)
    {
      fRebin = rebin;
    }
}

DataVector::DataVector(int size, const char* name, const char* title, int index, int rebin)
{
  fValues.ResizeTo(size);
  fStatisticalError.ResizeTo(size);
  fSystematicalErrorLow.ResizeTo(size);
  fSystematicalErrorHigh.ResizeTo(size);

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
}

DataVector::DataVector(TVectorD values)
{
  //this constructor uses the Sqrt of the values as statistical errors
  fValues = values;
  fStatisticalError.ResizeTo(fValues);
  fSystematicalErrorLow.ResizeTo(fValues);
  fSystematicalErrorHigh.ResizeTo(fValues);
  for(int i = 0; i < fValues.GetNoElements(); i++)
    {
      fStatisticalError[i] = Sqrt(fValues[i]);
      fSystematicalErrorLow[i] = 0.;
      fSystematicalErrorHigh[i] = 0.;
    }
}

DataVector::DataVector(TVectorD values, TVectorD error)
{
  //this constructor copies the error to statistical errors
  if(values.GetNoElements() != error.GetNoElements())
    {
      throw Form("%s: Error, sizes of value and error vector are different: %d != %d",__PRETTY_FUNCTION__,values.GetNoElements(),error.GetNoElements());
    }

  fValues = values;
  fStatisticalError.ResizeTo(fValues);
  fSystematicalErrorLow.ResizeTo(fValues);
  fSystematicalErrorHigh.ResizeTo(fValues);
  for(int i = 0; i < fValues.GetNoElements(); i++)
    {
      fStatisticalError[i] = error[i];
      fSystematicalErrorLow[i] = 0.;
      fSystematicalErrorHigh[i] = 0.;
    }
}

DataVector::DataVector(TVectorD values, TVectorD statisticalError, TVectorD systematicalError)
{
  if(values.GetNoElements() != statisticalError.GetNoElements() || values.GetNoElements() != systematicalError.GetNoElements())
    {
      throw Form("%s: Error, sizes of value and error vectors are different: %d != %d != %d",__PRETTY_FUNCTION__,values.GetNoElements(),statisticalError.GetNoElements(),systematicalError.GetNoElements());
    }

  fValues.ResizeTo(values);
  fStatisticalError.ResizeTo(values);
  fSystematicalErrorLow.ResizeTo(values);
  fSystematicalErrorHigh.ResizeTo(values);

  fValues = values;
  fStatisticalError = statisticalError;
  fSystematicalErrorLow = systematicalError;
  fSystematicalErrorHigh = systematicalError;
}

DataVector::~DataVector()
{
}

void DataVector::Clear()
{
  fValues.Clear();
  fStatisticalError.Clear();
  fSystematicalErrorLow.Clear();
  fSystematicalErrorHigh.Clear();

  fIndex = -999;
  fName.clear();
  fTitle.clear();
  fRebin = 0;
}

//calculate the errors from the values, if factor is negative use the sqrt
void DataVector::SetStatisticalError(double factor)
{
  if(factor < 0)
    {
      for(int i = 0; i < fValues.GetNoElements(); i++)
	{
	  fStatisticalError[i] = Sqrt(fValues[i]);
	}
    }
  else
    {
      for(int i = 0; i < fValues.GetNoElements(); i++)
	{
	  fStatisticalError[i] = factor*fValues[i];
	}
    }
}

void DataVector::SetSystematicalError(double factor)
{
  if(factor < 0)
    {
      for(int i = 0; i < fValues.GetNoElements(); i++)
	{
	  fSystematicalErrorLow[i] = Sqrt(fValues[i]);
	  fSystematicalErrorHigh[i] = Sqrt(fValues[i]);
	}
    }
  else
    {
      for(int i = 0; i < fValues.GetNoElements(); i++)
	{
	  fSystematicalErrorLow[i] = factor*fValues[i];
	  fSystematicalErrorHigh[i] = factor*fValues[i];
	}
    }
}

void DataVector::Set(TVectorD values)
{
  //use the Sqrt of the values as statistical errors
  fValues.ResizeTo(values);
  fStatisticalError.ResizeTo(values);
  fSystematicalErrorLow.ResizeTo(values);
  fSystematicalErrorHigh.ResizeTo(values);

  fValues = values;
  for(int i = 0; i < fValues.GetNoElements(); i++)
    {
      fStatisticalError[i] = Sqrt(fValues[i]);
      fSystematicalErrorLow[i] = 0.;
      fSystematicalErrorHigh[i] = 0.;
    }
}

void DataVector::Set(TVectorD values, TVectorD error)
{
  //copy the (symmetric) error to errorLow and errorHigh
  if(values.GetNoElements() != error.GetNoElements())
    {
      throw Form("%s: Error, sizes of value and error vector are different: %d != %d",__PRETTY_FUNCTION__,values.GetNoElements(),error.GetNoElements());
    }

  fValues = values;
  fStatisticalError.ResizeTo(fValues);
  fSystematicalErrorLow.ResizeTo(fValues);
  fSystematicalErrorHigh.ResizeTo(fValues);
  for(int i = 0; i < fValues.GetNoElements(); i++)
    {
      fStatisticalError[i] = error[i];
      fSystematicalErrorLow[i] = 0.;
      fSystematicalErrorHigh[i] = 0.;
    }
}

void DataVector::Set(TVectorD values, TVectorD statisticalError, TVectorD systematicalError)
{
  if(values.GetNoElements() != statisticalError.GetNoElements() || values.GetNoElements() != systematicalError.GetNoElements())
    {
      throw Form("%s: Error, sizes of value and error vectors are different: %d != %d != %d",__PRETTY_FUNCTION__,values.GetNoElements(),statisticalError.GetNoElements(),systematicalError.GetNoElements());
    }

  fValues.ResizeTo(values);
  fStatisticalError.ResizeTo(values);
  fSystematicalErrorLow.ResizeTo(values);
  fSystematicalErrorHigh.ResizeTo(values);

  fValues = values;
  fStatisticalError = statisticalError;
  fSystematicalErrorLow = systematicalError;
  fSystematicalErrorHigh = systematicalError;
}

void DataVector::Set(const char* name, DataVector vector, int index, int rebin)
{
  fValues.ResizeTo(vector.Size());
  fStatisticalError.ResizeTo(vector.Size());
  fSystematicalErrorLow.ResizeTo(vector.Size());
  fSystematicalErrorHigh.ResizeTo(vector.Size());

  fValues = *(vector.Values());
  fStatisticalError = *(vector.StatisticalError());
  fSystematicalErrorLow  = *(vector.SystematicalErrorLow());
  fSystematicalErrorHigh = *(vector.SystematicalErrorHigh());

  fName = name;
  fTitle = vector.Title();
  if(index != -999)
    {
      fIndex = index;
    }
  else
    {
      fIndex = vector.Index();
    }
  if(rebin != 0)
    {
      fRebin = rebin;
    }
}

void DataVector::SwapSystematicalErrors(int index)
{
  if(index < 0 || index > fSystematicalErrorLow.GetNoElements())
    {
      throw Form("%s: Error, index out of bounds 0 - %d",__PRETTY_FUNCTION__,fSystematicalErrorLow.GetNoElements());
    }

  swap(fSystematicalErrorLow[index],fSystematicalErrorHigh[index]);
}

void DataVector::ScaleStatisticalError(int index, double factor)
{
  if(index < 0 || index > fStatisticalError.GetNoElements())
    {
      throw Form("%s: Error, index out of bounds 0 - %d",__PRETTY_FUNCTION__,fStatisticalError.GetNoElements());
    }

  fStatisticalError[index] *= factor;
}

//fRebin:
//0 = sum the values  and use square root of sum of squares of errors as errors
//1 = average the values and use min/max as errors
//2 = average the values and use square root of sum of squares of errors as errors
void DataVector::Rebin(size_t RebinSize, vector<int> NumberOfPoints, int VerboseLevel)
{
  int i,j;
  double tmpValue;
  double tmpError;
  double tmpLow;
  double tmpHigh;
  double tmpIndex;
  int SumOfPoints;
  TVectorD values_Rebin(RebinSize);
  TVectorD statisticalError_Rebin(RebinSize);
  TVectorD systematicalErrorLow_Rebin(RebinSize);
  TVectorD systematicalErrorHigh_Rebin(RebinSize);

  //-------------------- rebin (according to NumberOfPoints) --------------------
  if(fRebin == 0)
    {
      //add all values together, error is sqrt of sum of squares of errors
      tmpValue = 0.;
      tmpError = 0.;
      tmpLow = 0.;
      tmpHigh = 0.;
      
      SumOfPoints = NumberOfPoints.at(0);
      for(i = 0, j = 0; i < fValues.GetNoElements(); i++)
	{
	  if(i == SumOfPoints)
	    {
	      if(VerboseLevel > 0)
	        {
	          cout<<"(i == "<<SumOfPoints<<") => "<<fName<<"_Rebin["<<j<<"] = "<<tmpValue<<endl;
	        }
	      values_Rebin[j] = tmpValue;
	      statisticalError_Rebin[j] = Sqrt(tmpError);
	      systematicalErrorLow_Rebin[j] = Sqrt(tmpLow);
	      systematicalErrorHigh_Rebin[j] = Sqrt(tmpHigh);
	      //reset all values
	      tmpValue = 0.;
	      tmpError = 0.;
	      tmpLow = 0.;
	      tmpHigh = 0.;
	      j++;
	      SumOfPoints += NumberOfPoints.at(j);
	    }//if(i == SumOfPoints)
	  if(VerboseLevel > 0)
	    {
	      cout<<fName<<"["<<i<<"] = "<<fValues[i]<<" +- "<<fStatisticalError[i]<<" + "<<fSystematicalErrorHigh[i]<<" - "<<fSystematicalErrorLow[i]<<", ";
	    }
	  tmpValue += fValues[i];
	  tmpError += Power(fStatisticalError[i],2);
	  tmpLow +=  Power(fSystematicalErrorLow[i],2);
	  tmpHigh += Power(fSystematicalErrorHigh[i],2);
	}//for(i = 0, j = 0; i < fValues.GetNoElements(); i++)
	  
      //store the last point
      if(VerboseLevel > 0)
        {
          cout<<"=> "<<fName<<"_Rebin["<<j<<"] = "<<tmpValue<<endl;
        }
      values_Rebin[j] = tmpValue;
      statisticalError_Rebin[j] = Sqrt(tmpError);
      systematicalErrorLow_Rebin[j] = Sqrt(tmpLow);
      systematicalErrorHigh_Rebin[j] = Sqrt(tmpHigh);
    }//if(fRebin == 0)
  else
    {
      //create an average of the values, rebin = 1 ignores any statistical errors
      //rebin > 10: ignore entries < 0
      tmpValue = 0.;
      tmpError = 0.;
      tmpIndex = 0;
      //rebin = 1 is only for values with systematical error only
      if(fRebin%10 == 1)
	{
	  tmpLow =  fValues[0]-fSystematicalErrorLow[0];
	  tmpHigh = fValues[0]+fSystematicalErrorHigh[0];
	}
      else
	{
	  tmpError = Power(fStatisticalError[0],2);
	  tmpLow =  Power(fSystematicalErrorLow[0],2);
	  tmpHigh = Power(fSystematicalErrorHigh[0],2);
	}
      SumOfPoints = NumberOfPoints.at(0);
      for(i = 0, j = 0; i < fValues.GetNoElements(); i++)
	{
	  //if this point is == SumOfPoints it is the first point of the next group => save the results and start anew
	  if(i == SumOfPoints)
	    {
	      if(VerboseLevel > 0)
	        {
	          cout<<"(i == "<<SumOfPoints<<") => "<<fName<<"_Rebin["<<j<<"] = "<<tmpValue/tmpIndex<<" = "<<tmpValue<<"/"<<tmpIndex<<endl;
	        }
	      values_Rebin[j] = tmpValue/tmpIndex;
	      if(fRebin == 1)
		{
		  statisticalError_Rebin[j] = 0.;
		  systematicalErrorLow_Rebin[j] = tmpValue/tmpIndex-tmpLow;
		  systematicalErrorHigh_Rebin[j] = tmpHigh-tmpValue/tmpIndex;
		}
	      else
		{
		  statisticalError_Rebin[j] = Sqrt(tmpError)/tmpIndex;
		  systematicalErrorLow_Rebin[j] = Sqrt(tmpLow)/tmpIndex;
		  systematicalErrorHigh_Rebin[j] = Sqrt(tmpHigh)/tmpIndex;
		}
	      //reset all values
	      tmpValue = fValues[i];
	      //set low and high to the actual values of upper and lower borders, not just their difference
	      if(fRebin == 1)
		{
		  tmpLow =  fValues[i]-fSystematicalErrorLow[i];
		  tmpHigh = fValues[i]+fSystematicalErrorHigh[i];
		}
	      else
		{
		  tmpError = Power(fStatisticalError[i],2);
		  tmpLow =  Power(fSystematicalErrorLow[i],2);
		  tmpHigh = Power(fSystematicalErrorHigh[i],2);
		}
	      tmpIndex = 1;
	      j++;
	      SumOfPoints += NumberOfPoints.at(j);
	    }//if(i == SumOfPoints)
	  else if(fValues[i] > 0. || fRebin < 10)//if the value is negative and rebin greate/equal 10 we do nothing
	    {
	      tmpValue += fValues[i];
	      tmpIndex++;
	      if(fRebin == 1)
		{
		  if(fValues[i] - fSystematicalErrorLow[i] < tmpLow)
		    {
		      tmpLow = fValues[i] - fSystematicalErrorLow[i];
		    }
		  if(fValues[i] + fSystematicalErrorHigh[i] > tmpHigh)
		    {
		      tmpHigh = fValues[i] + fSystematicalErrorHigh[i];
		    }
		}
	      else
		{
		  tmpError += Power(fStatisticalError[i],2);
		  tmpLow +=  Power(fSystematicalErrorLow[i],2);
		  tmpHigh += Power(fSystematicalErrorHigh[i],2);
		}
	    }
	      
	  if(VerboseLevel > 0)
	    {
	      cout<<fName<<"["<<i<<"] = "<<fValues[i]<<" + "<<fSystematicalErrorHigh[i]<<" - "<<fSystematicalErrorLow[i]<<", ";
	    }
	}//for(i = 0, j = 0; i < fValues.GetNoElements(); i++)
	  
      //store the last point
      if(VerboseLevel > 0)
        {
          cout<<"=> "<<fName<<"_Rebin.at("<<i<<")["<<j<<"] = "<<tmpValue/tmpIndex<<" = "<<tmpValue<<"/"<<tmpIndex<<endl;
        }
      values_Rebin[j] = tmpValue/tmpIndex;
      if(fRebin == 1)
	{
	  statisticalError_Rebin[j] = 0.;
	  systematicalErrorLow_Rebin[j] = tmpValue/tmpIndex-tmpLow;
	  systematicalErrorHigh_Rebin[j] = tmpHigh-tmpValue/tmpIndex;
	}
      else
	{
	  statisticalError_Rebin[j] = Sqrt(tmpError)/tmpIndex;
	  systematicalErrorLow_Rebin[j] = Sqrt(tmpLow)/tmpIndex;
	  systematicalErrorHigh_Rebin[j] = Sqrt(tmpHigh)/tmpIndex;
	}
    }//else of if(fRebin == 0)

  fValues.ResizeTo(values_Rebin.GetNoElements());
  fStatisticalError.ResizeTo(statisticalError_Rebin.GetNoElements());
  fSystematicalErrorLow.ResizeTo(systematicalErrorLow_Rebin.GetNoElements());
  fSystematicalErrorHigh.ResizeTo(systematicalErrorHigh_Rebin.GetNoElements());  
  
  fValues = values_Rebin;
  fStatisticalError = statisticalError_Rebin;
  fSystematicalErrorLow  = systematicalErrorLow_Rebin;
  fSystematicalErrorHigh = systematicalErrorHigh_Rebin;  
}

//remove entries that are zero and have zero error
void DataVector::PurgeZeros(vector<int> RemovedIndices)
{
  RemovedIndices.clear();
  for(int i = 0; i < fValues.GetNoElements(); i++)
    {
      //if all vectors are zero at this index remove it
      if(fValues[i] == 0. && fStatisticalError[i] == 0. && fSystematicalErrorLow[i] == 0. && fSystematicalErrorHigh[i] == 0.)
	{
	  cout<<"Error, this function is not yet fully implemented, i.e. it's doing nothing!"<<endl;
	}
    }
}

//-------------------- operators --------------------
//-------------------- DataVector --------------------
DataVector& DataVector::operator+= (const DataVector& rhs)
{
  //check for self-assignment
  if(this != &rhs)
    {
      if(this->Size() != rhs.Size())
	{
	  cout<<"Error, sizes of operands don't match: "<<this->Size()<<" != "<<rhs.Size()<<endl;

	  return *this;
	}

      fValues += rhs.fValues;
      for(int i = 0; i < fValues.GetNoElements(); i++)
	{
	  fStatisticalError[i] = Sqrt(Power(rhs.fStatisticalError[i],2)+Power(fStatisticalError[i],2));
	  fSystematicalErrorLow[i] = Sqrt(Power(rhs.fSystematicalErrorLow[i],2)+Power(fSystematicalErrorLow[i],2));
	  fSystematicalErrorHigh[i] = Sqrt(Power(rhs.fSystematicalErrorHigh[i],2)+Power(fSystematicalErrorHigh[i],2));
	}
    }

  return *this;
}

DataVector& DataVector::operator-= (const DataVector& rhs)
{
  //check for self-assignment
  if(this != &rhs)
    {
      if(this->Size() != rhs.Size())
	{
	  cout<<"Error, sizes of operands don't match: "<<this->Size()<<" != "<<rhs.Size()<<endl;

	  return *this;
	}

      fValues -= rhs.fValues;
      for(int i = 0; i < fValues.GetNoElements(); i++)
	{
	  fStatisticalError[i] = Sqrt(Power(rhs.fStatisticalError[i],2)+Power(fStatisticalError[i],2));
	  fSystematicalErrorLow[i] = Sqrt(Power(rhs.fSystematicalErrorLow[i],2)+Power(fSystematicalErrorLow[i],2));
	  fSystematicalErrorHigh[i] = Sqrt(Power(rhs.fSystematicalErrorHigh[i],2)+Power(fSystematicalErrorHigh[i],2));
	}
    }

  return *this;
}

DataVector& DataVector::operator*= (const DataVector& rhs)
{
  //check for self-assignment
  if(this != &rhs)
    {
      if(this->Size() != rhs.Size())
	{
	  cout<<"Error, sizes of operands don't match: "<<this->Size()<<" != "<<rhs.Size()<<endl;

	  return *this;
	}

      for(int i = 0; i < fValues.GetNoElements(); i++)
	{
	  fValues[i] *= rhs.fValues[i];
	  fStatisticalError[i] = Sqrt(Power(rhs.fStatisticalError[i]*fValues[i],2)+Power(fStatisticalError[i]*rhs.fValues[i],2));
	  fSystematicalErrorLow[i] = Sqrt(Power(rhs.fSystematicalErrorLow[i]*fValues[i],2)+Power(fSystematicalErrorLow[i]*rhs.fValues[i],2));
	  fSystematicalErrorHigh[i] = Sqrt(Power(rhs.fSystematicalErrorHigh[i]*fValues[i],2)+Power(fSystematicalErrorHigh[i]*rhs.fValues[i],2));
	}
    }

  return *this;
}

DataVector& DataVector::operator/= (const DataVector& rhs)
{
  //check for self-assignment
  if(this != &rhs)
    {
      if(this->Size() != rhs.Size())
	{
	  cout<<"Error, sizes of operands don't match: "<<this->Size()<<" != "<<rhs.Size()<<endl;

	  return *this;
	}

      for(int i = 0; i < fValues.GetNoElements(); i++)
	{
	  if(rhs.fValues[i] != 0)
	    {
	      fValues[i] /= rhs.fValues[i];
	      fStatisticalError[i] = Sqrt(Power(rhs.fStatisticalError[i]*fValues[i]/rhs.fValues[i]/rhs.fValues[i],2)+Power(fStatisticalError[i]/rhs.fValues[i],2));
	      fSystematicalErrorLow[i] = Sqrt(Power(rhs.fSystematicalErrorLow[i]*fValues[i]/rhs.fValues[i]/rhs.fValues[i],2)+Power(fSystematicalErrorLow[i]/rhs.fValues[i],2));
	      fSystematicalErrorHigh[i] = Sqrt(Power(rhs.fSystematicalErrorHigh[i]*fValues[i]/rhs.fValues[i]/rhs.fValues[i],2)+Power(fSystematicalErrorHigh[i]/rhs.fValues[i],2));
	    }
	  else
	    {
	      fValues[i] = 0.;
	      fStatisticalError[i] = 0.;
	      fSystematicalErrorLow[i] = 0.;
	      fSystematicalErrorHigh[i] = 0.;
	    }
	}
    }

  return *this;
}

//-------------------- TVectorD --------------------
//adding or subtracting a double (w/o error) leads to no change in the error
DataVector& DataVector::operator+= (const TVectorD& rhs)
{
  if(this->Size() != rhs.GetNoElements())
    {
      cout<<"Error, sizes of operands don't match: "<<this->Size()<<" != "<<rhs.GetNoElements()<<endl;
      
      return *this;
    }

  fValues += rhs;

  return *this;
}

DataVector& DataVector::operator-= (const TVectorD& rhs)
{
  if(this->Size() != rhs.GetNoElements())
    {
      cout<<"Error, sizes of operands don't match: "<<this->Size()<<" != "<<rhs.GetNoElements()<<endl;
      
      return *this;
    }

  fValues -= rhs;

  return *this;
}

DataVector& DataVector::operator*= (const TVectorD& rhs)
{
  if(this->Size() != rhs.GetNoElements())
    {
      cout<<"Error, sizes of operands don't match: "<<this->Size()<<" != "<<rhs.GetNoElements()<<endl;
      
      return *this;
    }

  for(int i = 0; i < fValues.GetNoElements(); i++)
    {
      fValues[i] *= rhs[i];
      fStatisticalError[i] *= rhs[i];
      fSystematicalErrorLow[i] *= rhs[i];
      fSystematicalErrorHigh[i] *= rhs[i];
    }

  return *this;
}

DataVector& DataVector::operator/= (const TVectorD& rhs)
{
  if(this->Size() != rhs.GetNoElements())
    {
      cout<<"Error, sizes of operands don't match: "<<this->Size()<<" != "<<rhs.GetNoElements()<<endl;
      
      return *this;
    }

  for(int i = 0; i < fValues.GetNoElements(); i++)
    {
      if(rhs[i] != 0.)
	{
	  fValues[i] /= rhs[i];
	  fStatisticalError[i] /= rhs[i];
	  fSystematicalErrorLow[i] /= rhs[i];
	  fSystematicalErrorHigh[i] /= rhs[i];
	}
      else
	{
	  fValues[i] = 0.;
	  fStatisticalError[i] = 0.;
	  fSystematicalErrorLow[i] = 0.;
	  fSystematicalErrorHigh[i] = 0.;
	}
    }

  return *this;
}

//-------------------- double --------------------
//adding or subtracting a double (w/o error) leads to no change in the error
DataVector& DataVector::operator+= (const double& rhs)
{
  fValues += rhs;

  return *this;
}

DataVector& DataVector::operator-= (const double& rhs)
{
  fValues -= rhs;

  return *this;
}

DataVector& DataVector::operator*= (const double& rhs)
{
  for(int i = 0; i < fValues.GetNoElements(); i++)
    {
      fValues[i] *= rhs;
      fStatisticalError[i] *= rhs;
      fSystematicalErrorLow[i] *= rhs;
      fSystematicalErrorHigh[i] *= rhs;
    }

  return *this;
}

DataVector& DataVector::operator/= (const double& rhs)
{
  if(rhs == 0)
    {
      cerr<<"Error, trying to divide DataVector by zero!"<<endl;

      return *this;
    }

  for(int i = 0; i < fValues.GetNoElements(); i++)
    {
      fValues[i] /= rhs;
      fStatisticalError[i] /= rhs;
      fSystematicalErrorLow[i] /= rhs;
      fSystematicalErrorHigh[i] /= rhs;
    }

  return *this;
}


//--------------------------------------------------------------------------------
//data graph class
DataGraph::DataGraph()
{
}

DataGraph::DataGraph(DataVector& x, DataVector& y)
{
  if(x.Size() != y.Size())
    {
      throw Form("%s: Error, sizes don't match: %d != %d",__PRETTY_FUNCTION__, x.Size(), y.Size());
    }

  if(strncmp(x.Name(),"RecoilThetaCm_",14) == 0)//compare the first 14 characters of the strings
    {
      if(strstr(x.Name(),"_Rebin") != NULL)
	{
	  fStatisticalGraph.SetName(Form("%s_vs_ThetaCm_Rebin", y.Name()));
	  fSystematicalGraph.SetName(Form("sys_%s_vs_ThetaCm_Rebin", y.Name()));
	}
      else
	{
	  fStatisticalGraph.SetName(Form("%s_vs_ThetaCm", y.Name()));
	  fSystematicalGraph.SetName(Form("sys_%s_vs_ThetaCm", y.Name()));
	}
    }
  else
    {
      fStatisticalGraph.SetName(Form("%s_vs_%s", y.Name(), x.Name()));
      fSystematicalGraph.SetName(Form("sys_%s_vs_%s", y.Name(), x.Name()));
    }
  fStatisticalGraph.SetTitle(Form("statistical errors of %s vs. %s", y.Title(), x.Title()));
  fSystematicalGraph.SetTitle(Form("systematical errors of %s vs. %s", y.Title(), x.Title()));

  fStatisticalGraph.Set(x.Size());
  fSystematicalGraph.Set(x.Size());

  for(int i = 0; i < x.Size(); i++)
    {
      fStatisticalGraph.SetPoint(i,x.Value(i),y.Value(i));
      fStatisticalGraph.SetPointError(i,x.StatisticalError(i),y.StatisticalError(i));
      fSystematicalGraph.SetPoint(i,x.Value(i),y.Value(i));
      fSystematicalGraph.SetPointEXlow(i,x.SystematicalErrorLow(i));
      fSystematicalGraph.SetPointEXhigh(i,x.SystematicalErrorHigh(i));
      fSystematicalGraph.SetPointEYlow(i,y.SystematicalErrorLow(i));
      fSystematicalGraph.SetPointEYhigh(i,y.SystematicalErrorHigh(i));
    }
}

DataGraph::DataGraph(DataVector& x, TVectorD& yValue, TVectorD& yError, const char* Name, const char* Title)
{
  if(x.Size() != yValue.GetNoElements())
    {
      throw Form("%s: Error, sizes don't match: %d != %d",__PRETTY_FUNCTION__, x.Size(), yValue.GetNoElements());
    }

  fStatisticalGraph.SetName(Name);
  fStatisticalGraph.SetTitle(Title);
  fSystematicalGraph.SetName(Form("sys_%s", Name));
  fSystematicalGraph.SetTitle(Form("systematical errors of %s", Title));

  fStatisticalGraph.Set(x.Size());
  fSystematicalGraph.Set(x.Size());

  for(int i = 0; i < x.Size(); i++)
    {
      fStatisticalGraph.SetPoint(i,x.Value(i),yValue[i]);
      fStatisticalGraph.SetPointError(i,x.StatisticalError(i),yError[i]);
      fSystematicalGraph.SetPoint(i,x.Value(i),yValue[i]);
      fSystematicalGraph.SetPointEXlow(i,x.SystematicalErrorLow(i));
      fSystematicalGraph.SetPointEXhigh(i,x.SystematicalErrorHigh(i));
      fSystematicalGraph.SetPointEYlow(i,0.);
      fSystematicalGraph.SetPointEYhigh(i,0.);
    }
}

DataGraph::~DataGraph()
{
}

void DataGraph::Draw(Option_t* option, const char* selection)
{
  if(selection != NULL)
    {
      if(strcmp(selection,"stat") == 0)
	{
	  fStatisticalGraph.Draw(option);
	}
      else if(strcmp(selection,"sys") == 0)
	{
	  fSystematicalGraph.Draw(option);
	}
      else
	{
	  cerr<<"selection "<<selection<<" is not implemented, choose either 'stat' or 'sys'"<<endl;
	}
    }
  else
    {
      //Option_t is a const char
      string Option;
      if(option != NULL)
	{
	  Option = option;
	  if(Option.find("[]") != string::npos)
	    {
	      Option.erase(Option.find("[]"),2);
	    }
	}
      else
	{
	  Option = "";
	}
      
      Option.insert(0,"a");
      
      fStatisticalGraph.Draw(Option.c_str());
      Option.erase(Option.begin());
      Option.append("[]");
      fSystematicalGraph.Draw(Option.c_str());
    }
}

void DataGraph::SetMarkerColor(Color_t color, const char* selection)
{
  if(selection != NULL)
    {
      if(strcmp(selection,"stat") == 0)
	{
	  fStatisticalGraph.SetMarkerColor(color);
	  return;
	}
      else if(strcmp(selection,"sys") == 0)
	{
	  fSystematicalGraph.SetMarkerColor(color);
	  return;
	}
      fStatisticalGraph.SetMarkerColor(color);
      fSystematicalGraph.SetMarkerColor(color);
      return;
    }
  else
    {
      fStatisticalGraph.SetMarkerColor(color);
      fSystematicalGraph.SetMarkerColor(color);
    }
}

void DataGraph::SetLineColor(Color_t color, const char* selection)
{
  if(selection != NULL)
    {
      if(strcmp(selection,"stat") == 0)
	{
	  fStatisticalGraph.SetLineColor(color);
	  return;
	}
      else if(strcmp(selection,"sys") == 0)
	{
	  fSystematicalGraph.SetLineColor(color);
	  return;
	}
      fStatisticalGraph.SetLineColor(color);
      fSystematicalGraph.SetLineColor(color);
      return;
    }
  else
    {
      fStatisticalGraph.SetLineColor(color);
      fSystematicalGraph.SetLineColor(color);
    }
}

void DataGraph::SetLineStyle(Style_t style, const char* selection)
{
  if(selection != NULL)
    {
      if(strcmp(selection,"stat") == 0)
	{
	  fStatisticalGraph.SetLineStyle(style);
	  return;
	}
      else if(strcmp(selection,"sys") == 0)
	{
	  fSystematicalGraph.SetLineStyle(style);
	  return;
	}
      fStatisticalGraph.SetLineStyle(style);
      fSystematicalGraph.SetLineStyle(style);
      return;
    }
  else
    {
      fStatisticalGraph.SetLineStyle(style);
      fSystematicalGraph.SetLineStyle(style);
    }
}

void DataGraph::SetLineWidth(Width_t width, const char* selection)
{
  if(selection != NULL)
    {
      if(strcmp(selection,"stat") == 0)
	{
	  fStatisticalGraph.SetLineWidth(width);
	  return;
	}
      else if(strcmp(selection,"sys") == 0)
	{
	  fSystematicalGraph.SetLineWidth(width);
	  return;
	}
      fStatisticalGraph.SetLineWidth(width);
      fSystematicalGraph.SetLineWidth(width);
      return;
    }
  else
    {
      fStatisticalGraph.SetLineWidth(width);
      fSystematicalGraph.SetLineWidth(width);
    }
}

int DataGraph::GetN(bool systematical)
{
  if(systematical)
    {
      return fSystematicalGraph.GetN();
    }
  return fStatisticalGraph.GetN();
}

int DataGraph::GetPoint(int bin, double& x, double& y, bool systematical)
{
  if(systematical)
    {
      return fSystematicalGraph.GetPoint(bin, x, y);
    }
  return fStatisticalGraph.GetPoint(bin, x, y);
}

double DataGraph::GetErrorX(int bin, bool systematical)
{
  if(systematical)
    {
      return fSystematicalGraph.GetErrorX(bin);
    }
  return fStatisticalGraph.GetErrorX(bin);
}

double DataGraph::GetErrorXlow(int bin, bool systematical)
{
  if(systematical)
    {
      return fSystematicalGraph.GetErrorXlow(bin);
    }
  return fStatisticalGraph.GetErrorXlow(bin);
}

double DataGraph::GetErrorXhigh(int bin, bool systematical)
{
  if(systematical)
    {
      return fSystematicalGraph.GetErrorXhigh(bin);
    }
  return fStatisticalGraph.GetErrorXhigh(bin);
}

double DataGraph::GetErrorY(int bin, bool systematical)
{
  if(systematical)
    {
      return fSystematicalGraph.GetErrorY(bin);
    }
  return fStatisticalGraph.GetErrorY(bin);
}

double DataGraph::GetErrorYlow(int bin, bool systematical)
{
  if(systematical)
    {
      return fSystematicalGraph.GetErrorYlow(bin);
    }
  return fStatisticalGraph.GetErrorYlow(bin);
}

double DataGraph::GetErrorYhigh(int bin, bool systematical)
{
  if(systematical)
    {
      return fSystematicalGraph.GetErrorYhigh(bin);
    }
  return fStatisticalGraph.GetErrorYhigh(bin);
}

TAxis* DataGraph::GetXaxis(bool systematical)
{
  if(systematical)
    {
      return fSystematicalGraph.GetXaxis();
    }
  return fStatisticalGraph.GetXaxis();
}

TAxis* DataGraph::GetYaxis(bool systematical)
{
  if(systematical)
    {
      return fSystematicalGraph.GetYaxis();
    }
  return fStatisticalGraph.GetYaxis();
}

void DataGraph::Set(DataVector x, DataVector y)
{
  if(x.Size() != y.Size())
    {
      throw Form("%s: Error, sizes don't match: %d != %d",__PRETTY_FUNCTION__, x.Size(), y.Size());
    }

  fStatisticalGraph.SetName(Form("%s_vs_%s", y.Name(), x.Name()));
  fStatisticalGraph.SetTitle(Form("%s vs. %s", y.Title(), x.Title()));
  fSystematicalGraph.SetName(Form("sys_%s_vs_%s", y.Name(), x.Name()));
  fSystematicalGraph.SetTitle(Form("systematical errors of %s vs. %s", y.Title(), x.Title()));

  fStatisticalGraph.Set(x.Size());
  fSystematicalGraph.Set(x.Size());

  for(int i = 0; i < x.Size(); i++)
    {
      fStatisticalGraph.SetPoint(i,x.Value(i),y.Value(i));
      fStatisticalGraph.SetPointError(i,x.StatisticalError(i),y.StatisticalError(i));
      fSystematicalGraph.SetPoint(i,x.Value(i),y.Value(i));
      fSystematicalGraph.SetPointEXlow(i,x.SystematicalErrorLow(i));
      fSystematicalGraph.SetPointEXhigh(i,x.SystematicalErrorHigh(i));
      fSystematicalGraph.SetPointEYlow(i,y.SystematicalErrorLow(i));
      fSystematicalGraph.SetPointEYhigh(i,y.SystematicalErrorHigh(i));
    }
}

void DataGraph::SetName(const char* name)
{
  fStatisticalGraph.SetName(name);
  fSystematicalGraph.SetName(Form("sys_%s", name));
}

void DataGraph::SetTitle(const char* title)
{
  fStatisticalGraph.SetTitle(title);
  fSystematicalGraph.SetTitle(Form("systematical errors of %s", title));
}

const char* DataGraph::Name()
{
  return fStatisticalGraph.GetName();
}

const char* DataGraph::Title()
{
  return fStatisticalGraph.GetTitle();
}

void DataGraph::Scale(double scale)
{
  double x,y;
  for(int i = 0; i < fStatisticalGraph.GetN(); i++)
    {
      fStatisticalGraph.GetPoint(i,x,y);
      fStatisticalGraph.SetPoint(i,x,y*scale);
      fStatisticalGraph.SetPointError(i,fStatisticalGraph.GetErrorX(i)*scale,fStatisticalGraph.GetErrorY(i)*scale);
      fSystematicalGraph.GetPoint(i,x,y);
      fSystematicalGraph.SetPoint(i,x,y*scale);
      fSystematicalGraph.SetPointError(i,fSystematicalGraph.GetErrorXlow(i),fSystematicalGraph.GetErrorXhigh(i),
				       fSystematicalGraph.GetErrorYlow(i)*scale,fSystematicalGraph.GetErrorYhigh(i)*scale);
    }
}

void DataGraph::SetPoint(int bin, double x, double y, bool systematical)
{
  if(systematical)
    {
      fSystematicalGraph.SetPoint(bin, x, y);
    }
  fStatisticalGraph.SetPoint(bin, x, y);
}

void DataGraph::SetPointError(int bin, double ex, double ey)
{
  fStatisticalGraph.SetPointError(bin, ex, ey);
}

void DataGraph::SetPointError(int bin, double exLow, double exHigh, double eyLow, double eyHigh)
{
  fSystematicalGraph.SetPointError(bin, exLow, exHigh, eyLow, eyHigh);
}

