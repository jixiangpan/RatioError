#include<iostream>
#include<fstream>
#include<cmath>
#include "stdlib.h"
using namespace std;

#include<map>
#include<vector>
#include<set>

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THStack.h"
#include "TF1.h"
#include "TLine.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"

#include "TRandom3.h"
#include "TGaxis.h"
#include "TStyle.h"

#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TROOT.h"

#include "TF1Convolution.h"

void func_canv_margin(TCanvas *canv, double left, double right, double top, double bot)
{
  canv->SetLeftMargin(left);
  canv->SetRightMargin(right);
  canv->SetTopMargin(top);
  canv->SetBottomMargin(bot);
}

/////////////////////////////////////////////////////////////////////////////////// class

class TAnalysis
{
 public:
  TAnalysis() {
    num_component = 0;
    
    data_range_low = 0;
    data_range_hgh = 0;

    ratio_range_low = 0;
    ratio_range_hgh = 0;
    
    val_meas = 0;
    val_pred = 0;
    
    ratio_1sigma_low = 0;
    ratio_1sigma_hgh = 0;
     
  }

  static double pdf_func_sterling_poisson(double *x, double *par);
  static double func_Y2X(double *x, double *par);
  static double pdf_func_Y2X(double *x, double *par);
    
  void Add_component(double num, double weight);
  TF1 *Get_component(int i) {return  map_pdf_func_component[i];}
  int Get_num_component() { return num_component; }
  
  TF1 *Summation_func(int i, int j);// summation of two components
  void Summation_total_func();
  TF1 *Get_func_summation_pred() { return func_summation_pred; };

  void Func_ratio_meas2pred();
  TF1 *Get_Func_ratio_meas2pred() { return func_ratio_meas2pred; };

  void Get_ratio_lower_upper();
  
  void Clear();

  double data_range_low;
  double data_range_hgh;

  double ratio_range_low;
  double ratio_range_hgh;
  
  double val_meas;
  double val_pred;

  double ratio_1sigma_low;
  double ratio_1sigma_hgh;
  
 private:
  int num_component;
  map<int, TF1*>map_pdf_func_component;
  double map_map_pars[100][4];

  TF1 *func_summation_pred;
  TF1 *func_ratio_meas2pred;

  static TF1 *fX_func;
  static TF1 *fY_func;
};

TF1* TAnalysis::fX_func = NULL;
TF1* TAnalysis::fY_func = NULL;

////
////

void TAnalysis::Clear()
{
  num_component = 0;
  
  map_pdf_func_component.clear();
    
  for(int i=0; i<100; i++)
    for(int j=0; j<4; j++)
      map_map_pars[i][j] = 0;

  func_summation_pred = NULL;

  data_range_low = 0;
  data_range_hgh = 0;

  ratio_range_low = 0;
  ratio_range_hgh = 0;
  
  val_meas = 0;
  val_pred = 0;

  ratio_1sigma_low = 0;
  ratio_1sigma_hgh = 0;
}

double TAnalysis::pdf_func_sterling_poisson(double *x, double *par)
{  
  double val_meas   = par[0];
  double val_weight = par[1];
  double xrange_low = par[2];
  double xrange_hgh = par[3];

  double xcur = x[0];
  double val_pred = xcur/val_weight;

  double result = 0;

  if( xcur>=xrange_low && xcur<=xrange_hgh ) {
    if( val_meas>0 ) {
      // https://dlmf.nist.gov/5.11#i
      // https://en.wikipedia.org/wiki/Stirling%27s_approximation
      const int cn = 6;
      double ci[cn] = {1., 1./12, 1./288, -139./51840, -571./2488320, 163879./209018880};
      for(int i=0; i<cn; i++) {
        result += ci[i]/pow( val_meas, i );
      }
      result = exp( val_meas - val_pred + val_meas*log(val_pred/val_meas) ) / ( sqrt(2*TMath::Pi()*val_meas) * result  ) /val_weight;
    }
    else {
      result = exp(-val_pred);
    }
  }
  
  return result;
}


double TAnalysis::func_Y2X(double *x, double *par)
{
  double result = 0;
  
  double xcur = x[0];
  double Y2X = par[0];
  double range_low = par[1];
  double range_hgh = par[2];
  
  if( xcur<range_low || xcur>range_hgh ) {
    result = 0;
  }
  else {
    double eval_fX = fX_func->Eval( xcur );
    double eval_fY = fY_func->Eval( Y2X * xcur );
    double abs_X = fabs( xcur );
    result = abs_X * eval_fX * eval_fY;
  }
  
  return result;
}

double TAnalysis::pdf_func_Y2X(double *x, double *par)
{
  double result = 0;

  double z = x[0];
  double range_low = par[0];
  double range_hgh = par[1];

  TF1 *roofunc_Y2X = new TF1("roofunc_Y2X", func_Y2X, range_low, range_hgh, 3);
  roofunc_Y2X->SetParameter(0, z);
  roofunc_Y2X->SetParameter(1, range_low);
  roofunc_Y2X->SetParameter(2, range_hgh);
  roofunc_Y2X->SetNpx(6000);
  
  /// Exact solution (cannot reach tolerance because of roundoff error):
  double val_integration = roofunc_Y2X->Integral( range_low, range_hgh );    
  result = val_integration;

  delete roofunc_Y2X;  
  return result;  
}


/////////////////////////////

void TAnalysis::Add_component(double num, double weight)
{
  num_component++;
  map_map_pars[num_component][0] = num;
  map_map_pars[num_component][1] = weight;
  map_map_pars[num_component][2] = data_range_low;
  map_map_pars[num_component][3] = data_range_hgh;  
  
  TString roostr = TString::Format("map_pdf_func_component_%02d", num_component);
  map_pdf_func_component[num_component] = new TF1( roostr, pdf_func_sterling_poisson, data_range_low, data_range_hgh, 4 );
  map_pdf_func_component[num_component]->SetParameters( num, weight, data_range_low, data_range_hgh );
  map_pdf_func_component[num_component]->SetNpx(60000);

  if(num_component==1) {
    val_meas = num * weight;
  }
  else {
    val_pred += num * weight;
  }
  
}

TF1 *TAnalysis::Summation_func(int i, int j)
{
  TF1Convolution *fconv = new TF1Convolution(map_pdf_func_component[i], map_pdf_func_component[j], data_range_low, data_range_hgh, true);
  fconv->SetNofPointsFFT(10000);

  TF1 *func_temp = new TF1("func_temp", fconv, data_range_low, data_range_hgh, 0);
  func_temp->SetNpx(60000);

  return func_temp;
}

void TAnalysis::Summation_total_func()
{
  map<int, TF1Convolution*>fconv;
  map<int, TF1*>func_conv;
  TString roostr = "";
  
  if( num_component<=2 ) {
    if( num_component==2 ) {
      cerr<<" Warning: num_component of prediction is 1"<<endl;
      func_summation_pred = map_pdf_func_component[2];
    }
    else {
      cerr<<" Error: no predction"<<endl;
    }
  }
  else {
    int line = 3;
    fconv[line] = new TF1Convolution( map_pdf_func_component[line-1], map_pdf_func_component[line], data_range_low, data_range_hgh, true );
    fconv[line]->SetNofPointsFFT(10000);
    roostr = TString::Format("func_conv_%02d", line);
    func_conv[line] = new TF1(roostr, fconv[line], data_range_low, data_range_hgh, 0);
    func_conv[line]->SetNpx(60000);
  
    if( num_component>3 ) {
      for( int idx=4; idx<=num_component; idx++ ) {
	line++;
	fconv[line] = new TF1Convolution( func_conv[line-1], map_pdf_func_component[line], data_range_low, data_range_hgh, true );
	fconv[line]->SetNofPointsFFT(10000);
	roostr = TString::Format("func_conv_%02d", line);
	func_conv[line] = new TF1(roostr, fconv[line], data_range_low, data_range_hgh, 0);
	func_conv[line]->SetNpx(60000);
      }      
    }

    func_summation_pred = func_conv[line];    
  }    
}

void TAnalysis::Func_ratio_meas2pred()
{
  TAnalysis::fY_func = map_pdf_func_component[1];
  
  if( num_component==2 ) TAnalysis::fX_func = map_pdf_func_component[2];
  else TAnalysis::fX_func = func_summation_pred;

  TF1 *func_ratio_Y2X = new TF1("pdf_func_Y2X", TAnalysis::pdf_func_Y2X, ratio_range_low, ratio_range_hgh, 2);
  func_ratio_Y2X->SetParameters(data_range_low, data_range_hgh);
  func_ratio_Y2X->SetNpx(60000);
  func_ratio_meas2pred = func_ratio_Y2X;
}

void TAnalysis::Get_ratio_lower_upper()
{
  gErrorIgnoreLevel = 6001;
  
  double val_ratio = val_meas/val_pred;
  double value_1sigma = 1-TMath::Prob(1, 1);

  double central_ratio = val_ratio;
  
  int line = 0;
  double val_typeA = 0;
  double val_typeB = 0;

  ////////////////////////

  line = 0;
  val_typeA = ratio_range_low;
  val_typeB = central_ratio;
  while( fabs(val_typeA-val_typeB)>1e-4 ) {
    line++;
    double val_mid = (val_typeA+val_typeB)/2;
    double y_mid = func_ratio_meas2pred->Integral(val_mid, central_ratio);
    if( y_mid > value_1sigma/2 ) {
      val_typeA = val_mid;
    }
    else {
      val_typeB = val_mid;
    }
    //cout<<TString::Format(" -------> %3d, %10.6f %10.6f", line, val_typeA, val_typeB)<<endl;
    if( line>30 ) {
      cerr<<" Error: cannot find ratio_1sigma_low"<<endl;
      break;
    }
  }
  
  ratio_1sigma_low = (val_typeA+val_typeB)/2;
     
  ////////////////////////
  
  line = 0;
  val_typeA = central_ratio;
  val_typeB = ratio_range_hgh;
  while( fabs(val_typeA-val_typeB)>1e-4 ) {
    line++;
    double val_mid = (val_typeA+val_typeB)/2;
    double y_mid = func_ratio_meas2pred->Integral(central_ratio, val_mid);
    if( y_mid > value_1sigma/2 ) {
      val_typeB = val_mid;
    }
    else {
      val_typeA = val_mid;
    }
    //cout<<TString::Format(" -------> %3d, %10.6f %10.6f", line, val_typeA, val_typeB)<<endl;
    if( line>30 ) {
      cerr<<" Error: cannot find ratio_1sigma_hgh"<<endl;
      break;
    }
  }
  
  ratio_1sigma_hgh = (val_typeA+val_typeB)/2;
     
  ////////////////////////  
  
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// MAIN ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

void read_cal()
{
  //ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1.E-5); 
  //gErrorIgnoreLevel = 6001;
  
  TString roostr = "";

  cout<<endl<<" Hello World"<<endl<<endl;

  TAnalysis *exampleA = new TAnalysis();
  exampleA->data_range_low  = 0;
  exampleA->data_range_hgh  = 600;
  exampleA->ratio_range_low = 0;
  exampleA->ratio_range_hgh = 2;
  
  exampleA->Add_component(460, 1);// the first is for measurment
  exampleA->Add_component(400, 0.25);// the following is for prediciton component
  exampleA->Add_component(200, 1);
  exampleA->Add_component(300, 0.5);

  exampleA->Summation_total_func();
  
  /////////////////////////////////
  
  roostr = "canv_testA";
  TCanvas *canv_testA = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_testA, 0.15, 0.2,0.1,0.15);
    
  (exampleA->Get_component(2))->Draw();
  (exampleA->Get_component(2))->SetLineColor(kBlue);
  
  (exampleA->Get_component(3))->Draw("same");
  (exampleA->Get_component(3))->SetLineColor(kRed);
   
  (exampleA->Get_component(4))->Draw("same");
  (exampleA->Get_component(4))->SetLineColor(kGreen+1);
   
  (exampleA->Get_component(1))->Draw("same");
  (exampleA->Get_component(1))->SetLineColor(kGray+1);

  (exampleA->Get_func_summation_pred())->Draw("same");
  (exampleA->Get_func_summation_pred())->SetLineColor(kBlack);
  
  ///////////////////////////////////
  
  exampleA->Func_ratio_meas2pred();
  
  exampleA->Get_ratio_lower_upper();

  cout<<endl<<TString::Format(" Ratio: %5.3f, (%5.3f, %5.3f)",
			      exampleA->val_meas/exampleA->val_pred,
			      exampleA->ratio_1sigma_low,
			      exampleA->ratio_1sigma_hgh
			      )<<endl<<endl;
  
  roostr = "canv_testB";
  TCanvas *canv_testB = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_testB, 0.15, 0.2,0.1,0.15);  
  (exampleA->Get_Func_ratio_meas2pred())->Draw();
  cout<<endl<<" Integral "<<(exampleA->Get_Func_ratio_meas2pred())->Integral(0,2)<<endl<<endl;
  
  /////////
/*  
  exampleA->Clear();
  
  exampleA->data_range_low  = 0;
  exampleA->data_range_hgh  = 600;
  exampleA->ratio_range_low = 0;
  exampleA->ratio_range_hgh = 2;
  
  exampleA->Add_component(400, 1);// the first is for measurment
  exampleA->Add_component(200, 1);// the following is for prediciton component
  exampleA->Add_component(100, 2);

  exampleA->Summation_total_func();
  
  exampleA->Func_ratio_meas2pred();
  
  exampleA->Get_ratio_lower_upper();

  cout<<endl<<TString::Format(" Ratio: %5.3f, (%5.3f, %5.3f)",
			      exampleA->val_meas/exampleA->val_pred,
			      exampleA->ratio_1sigma_low,
			      exampleA->ratio_1sigma_hgh
			      )<<endl<<endl;
  
  cout<<endl<<" Integral "<<(exampleA->Get_Func_ratio_meas2pred())->Integral(0,2)<<endl<<endl;
*/
}
