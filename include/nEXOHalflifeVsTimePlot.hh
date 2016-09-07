#ifndef nEXOHalflifeVsTimePlot_hh
#define nEXOHalflifeVsTimePlot_hh

#include <iostream>
#include <map>

#include "TNamed.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TGaxis.h"
#include "TPaveText.h"
#include "TText.h"
#include "TLegend.h"
#include "TF1.h"

#include "nEXOSensPlot.hh"
#include "nEXONuclearMatrixElement.hh"
#include "nEXOUtils.hh"

class nEXOHalflifeVsTimePlot : public nEXOSensPlot
{
public:
  nEXOHalflifeVsTimePlot(const char* name = 0, const char* title = 0);
  virtual ~nEXOHalflifeVsTimePlot();

  bool SetGraphPoints(const char* name, const char* title, size_t n, double* yrs, double* hls, double unit = 1e26);
  
  void SetAxisLimits(Double_t xMin=0., Double_t xMax=10., Double_t yMin=1.e25, Double_t yMax=2.e29);
  void SetNME(nEXONuclearMatrixElement::NME_t nme);

  TObject* GetPlot();

  TCanvas* CreateEmptyCanvas();
  TGraph* GetNormalGraph(nEXONuclearMatrixElement& nme);
  TGraph* GetInvertedGraph(nEXONuclearMatrixElement& nme);

  std::map<TString, Int_t> fLineColors;
  std::map<TString, Int_t> fLineWidths;
  std::map<TString, Int_t> fLineStyles;
  std::map<TString, TF1*> fSmoothFunction;
  
protected:

  nEXONuclearMatrixElement* fNME;

  std::map<TString, TGraph> fGraphs;
    
  TString fXaxisTitle;
  TString fYaxisTitle;
  
  Double_t fXmin;
  Double_t fXmax;
  Double_t fYmin;
  Double_t fYmax;
  
  void CreateAvailableGraphs();
  void SetGraphProperties(const char* name);
  void AddGraph(TString name, TString title, TGraph& graph);
  void PlotGraph(TGraph& graph, TCanvas& canvas, TLegend& leg, TF1* smooth = 0);
  void PlotEXO200(TCanvas& canvas, TLegend& leg);
  
  ClassDef(nEXOHalflifeVsTimePlot,1) // Base class for nEXO sensitivity plots
};

#endif

