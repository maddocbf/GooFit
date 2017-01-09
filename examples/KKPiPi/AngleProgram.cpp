#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TTreeReader.h"
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <sstream>
#include "math.h"
#include "tgmath.h"
#include <complex>
#include <iostream>
#include <cmath>
#include "TGraph.h"
#include "TH2.h"
#include "TGraphErrors.h"
#include <vector>
#include <array>
#include "TChain.h"
#include "TLorentzVector.h"
#include "fit1MeV_Gaussian.C"
#include "TLegend.h"
#include "TPaveText.h"
#include "TColor.h"

double phiCalculation(TVector3 pAB, TVector3 pCD, TVector3 pA, TVector3 pB, TVector3 pC, TVector3 pD)
{
    TVector3 zhat = pAB.Unit();
    TVector3 yhat = (pA.Cross(pB)).Unit();
    TVector3 xhat = (yhat.Cross(zhat)).Unit();

    TVector3 yhatPrime = (pC.Cross(pD)).Unit();

    Double_t cosPhi = (yhat.Dot(yhatPrime));
    Double_t sinPhi = (xhat.Dot(yhatPrime));
     
    Double_t phi    = TMath::ATan2(sinPhi,cosPhi);
    Double_t altPhi = TMath::ATan2(sinPhi,cosPhi);

    if (phi < 0)
    {
      phi = 2.0*TMath::Pi()+phi;
    }
    
    return phi; 
}

double sin2ThetaAfn(TVector3 pAprime_3vec, TVector3 pAB, TVector3 pCD, TVector3 pA, TVector3 pB, TVector3 pC, TVector3 pD)
{
      TVector3 zhat = pAB.Unit();
      TVector3 yhat = (pA.Cross(pB)).Unit();
      TVector3 xhat = (yhat.Cross(zhat)).Unit();

      Double_t paPrimeX = pAprime_3vec.Dot(xhat);
      Double_t paPrimeY = pAprime_3vec.Dot(yhat);
      Double_t paPrimeZ = pAprime_3vec.Dot(zhat);
      
      Double_t paPrimeMag = sqrt(paPrimeX*paPrimeX+paPrimeY*paPrimeY+paPrimeZ*paPrimeZ);
      Double_t paPrimePT  = sqrt(paPrimeX*paPrimeX+paPrimeY*paPrimeY);
      Double_t cosThetaA  = paPrimeZ/paPrimeMag;
      Double_t sinThetaA  = paPrimePT/paPrimeMag;
      Double_t sin2ThetaA = 2.0*sinThetaA*cosThetaA;

      return sin2ThetaA; 
}

double sin2ThetaCfn(TVector3 pCprime_3vec, TVector3 pAB, TVector3 pCD, TVector3 pA, TVector3 pB, TVector3 pC, TVector3 pD)
{

    TVector3 zhat = pAB.Unit();
    TVector3 yhat = (pA.Cross(pB)).Unit();
    TVector3 xhat = (yhat.Cross(zhat)).Unit();

    Double_t pcPrimeX = pCprime_3vec.Dot(xhat);
    Double_t pcPrimeY = pCprime_3vec.Dot(yhat);
    Double_t pcPrimeZ = pCprime_3vec.Dot(zhat);
      
    Double_t pcPrimeMag = sqrt(pcPrimeX*pcPrimeX+pcPrimeY*pcPrimeY+pcPrimeZ*pcPrimeZ);
    Double_t pcPrimePT  = sqrt(pcPrimeX*pcPrimeX+pcPrimeY*pcPrimeY);
    Double_t cosThetaC  = pcPrimeZ/pcPrimeMag;
    Double_t sinThetaC  = pcPrimePT/pcPrimeMag;
    Double_t sin2ThetaC = 2.0*sinThetaC*cosThetaC;

    return sin2ThetaC; 
}

 
int main(int argc, char *argv[])
{ 

  TString DecayOne = "KK PiPi"; 
  TString DecayTwo = " KPi KPi"; 

  TCanvas* mAB2mCD2Canvas = new TCanvas("mAB2mCD2Canvas", "mAB2mCD2Canvas", 10,10,800,600); 
  TCanvas* mAD2mBC2Canvas = new TCanvas("mAD2mBC2Canvas", "mAD2mBC2Canvas", 10,10, 800,600); 
  TCanvas* mAB2mBC2Canvas = new TCanvas("mAB2mBC2Canvas", "mAB2mBC2Canvas", 10,10,800,600);
  TCanvas* mAB2mAD2Canvas = new TCanvas("mAB2mAD2Canvas", "mAB2mAD2Canvas", 10,10,800,600);
  TCanvas* mCD2mBC2Canvas = new TCanvas("mCD2mBC2Canvas", "mCD2mBC2Canvas", 10,10,800,600);
  TCanvas* mCD2mAD2Canvas = new TCanvas("mCD2mAD2Canvas", "mCD2mAD2Canvas", 10,10,800,600);
  TCanvas* mAB2mCD2ResCanvas = new TCanvas("mAB2mCD2ResCanvas", "mAB2mCD2ResCanvas", 10, 10, 800,600); 
  TCanvas* mAD2mBC2ResCanvas = new TCanvas("mAD2mBC2ResCanvas", "mAD2mBC2ResCanvas", 10, 10, 800,600); 
  TCanvas* mAB2mCD2LogCanvas = new TCanvas("mAB2mCD2LogCanvas", "mAB2mCD2LogCanvas", 10,10,800,600); 
  TCanvas* mAD2mBC2LogCanvas = new TCanvas("mAD2mBC2LogCanvas", "mAD2mBC2LogCanvas", 10,10,800,600);  
  TCanvas* mAB2mBC2LogCanvas = new TCanvas("mAB2mBC2LogCanvas", "mAB2mBC2LogCanvas", 10,10,800,600);
  TCanvas* mAB2mAD2LogCanvas = new TCanvas("mAB2mAD2LogCanvas", "mAB2mAD2LogCanvas", 10,10,800,600);
  TCanvas* mCD2mBC2LogCanvas = new TCanvas("mCD2mBC2LogCanvas", "mCD2mBC2LogCanvas", 10,10,800,600);
  TCanvas* mCD2mAD2LogCanvas = new TCanvas("mCD2mAD2LogCanvas", "mCD2mAD2LogCanvas", 10,10,800,600);  

  TH2* mAB2mCD2Hist = new TH2D("mAB2mCD2Hist", "(KK) (#pi#pi) 2D Histogram; m(KK)^{2} GeV^{2}; m(#pi #pi)^{2} GeV^{2}", 100,0.9,2.5, 100,0.0,1.0); 
  TH2* mAD2mBC2Hist = new TH2D("mAD2mBC2Hist", "(K#pi) (K#pi) 2D Histogram; m(K^{-}#pi^{+})^{2} GeV^{2}; m(K^{+}#pi^{-})^{2} GeV^{2}", 100, 0.2,1.8, 100,0.2,1.8);
  TH2* mAB2mBC2Hist = new TH2D("mAB2mBC2Hist", "(KK) (K^{+}#pi^{-}) 2D Histogram; m(KK)^{2} GeV^{2}; m(K^{+} #pi^{-})^{2} GeV^{2}", 100,0.8,2.5, 100,0.2,1.8); 
  TH2* mAB2mAD2Hist = new TH2D("mAB2mAD2Hist", "(KK) (K^{-}#pi^{+}) 2D Histogram; m(KK)^{2} GeV^{2}; m(K^{-} #pi^{+})^{2} GeV^{2}", 100,0.8,2.5, 100,0.2,1.8);
  TH2* mCD2mBC2Hist = new TH2D("mCD2mBC2Hist", "(#pi#pi) (K^{+} #pi^{-}) 2D Histogram; m(#pi #pi)^{2} GeV^{2}; m(K^{+} #pi^{-})^{2} GeV^{2}", 100,0.0,1.0, 100,0.2,1.8); 
  TH2* mCD2mAD2Hist = new TH2D("mCD2mAD2Hist", "(#pi#pi) (K^{-} #pi^{+}) 2D Histogram; m(#pi #pi)^{2} GeV^{2}; m(K^{-} #pi^{+})^{2} GeV^{2}", 100,0.0,1.0, 100,0.2,1.8);
  TH2* mAB2mCD2ResHist = new TH2D("mAB2mCD2ResHist", "(KK) (#pi#pi) 2D Histogram; m(KK)^{2} GeV^{2}; m(#pi #pi)^{2} GeV^{2}", 100,0.9,2.5, 100,0.0,1.0); 
  TH2* mAD2mBC2ResHist = new TH2D("mAD2mBC2ResHist", "(K#pi) (K#pi) 2D Histogram; m(K^{-}#pi^{+})^{2} GeV^{2}; m(K^{+}#pi^{-})^{2} GeV^{2}", 100, 0.2,1.8, 100,0.2,1.8);
  TH2* mAB2mCD2LogHist = new TH2D("mAB2mCD2LogHist", "(KK) (#pi#pi) 2D Histogram; m(KK)^{2} GeV^{2}; m(#pi #pi)^{2} GeV^{2}", 100,0.9,2.5, 100,0.0,1.); 
  TH2* mAD2mBC2LogHist = new TH2D("mAD2mBC2LogHist", "(K#pi) (K#pi) 2D Histogram; m(K^{-}#pi^{+})^{2} GeV^{2}; m(K^{+}#pi^{-})^{2} GeV^{2}", 100, 0.2,1.8, 100,0.2,1.8);
  TH2* mAB2mBC2LogHist = new TH2D("mAB2mBC2LogHist", "(KK) (K^{+}#pi^{-}) 2D Histogram; m(KK)^{2} GeV^{2}; m(K^{+} #pi^{-})^{2} GeV^{2}", 100,0.8,2.5, 100,0.2,1.8); 
  TH2* mAB2mAD2LogHist = new TH2D("mAB2mAD2LogHist", "(KK) (K^{-}#pi^{+}) 2D Histogram; m(KK)^{2} GeV^{2}; m(K^{-} #pi^{+})^{2} GeV^{2}", 100,0.8,2.5, 100,0.2,1.8);
  TH2* mCD2mBC2LogHist = new TH2D("mCD2mBC2LogHist", "(#pi#pi) (K^{+} #pi^{-}) 2D Histogram; m(#pi #pi)^{2} GeV^{2}; m(K^{+} #pi^{-})^{2} GeV^{2}", 100,0.0,1.0, 100,0.2,1.8); 
  TH2* mCD2mAD2LogHist = new TH2D("mCD2mAD2LogHist", "(#pi#pi) (K^{-} #pi^{+}) 2D Histogram; m(#pi #pi)^{2} GeV^{2}; m(K^{-} #pi^{+})^{2} GeV^{2}", 100,0.0,1.0, 100,0.2,1.8);

  

  TCanvas* mABCanvas = new TCanvas("mABCanvas", "mABCanvas", 10,10,800,600); 
  TCanvas* mCDCanvas = new TCanvas("mCDCanvas", "mCDCanvas", 10,10,800,600); 
  TCanvas* mADCanvas = new TCanvas("mADCanvas", "mADCanvas", 10,10,800,600); 
  TCanvas* mBCCanvas = new TCanvas("mBCCanvas", "mBCCanvas", 10,10,800,600); 
  TH1* mABHist = new TH1D("mABHist", "m(KK) Histogram; m(KK) MeV; events per 1 MeV", 650, 950, 1600); 
  TH1* mCDHist = new TH1D("mCDHist", "m(#pi#pi) Histogram; m(#pi#pi) MeV; events per 1 MeV", 700, 200, 900); 
  TH1* mADHist = new TH1D("mADHist", "m(K^{-}#pi^{+}) Histogram; m(K^{-}#pi^{+}) MeV; events per 1 MeV", 600, 600, 1200); 
  TH1* mBCHist = new TH1D("mBCHist", "m(K^{+}#pi^{-}) Histogram; m(K^{+}#pi^{-}) MeV; events per 1 MeV", 600, 600, 1200); 

 

  TCanvas* cosThetaACanvas = new TCanvas("CosThetaACanvas", "CosThetaACanvas", 10,10,800,600); 
  TCanvas* cosThetaCCanvas = new TCanvas("CosThetaCCanvas", "CosThetaCCanvas", 10,10,800,600); 
  TCanvas* cosThetaACCanvas = new TCanvas("CosThetaACCanvas", "CosThetaACCanvas", 10,10,800,600); 
  TH1* cosThetaAHist = new TH1D("cosThetaAHist", "cosThetaAHist; #phi(radians); events per", 100, -1,1); 
  TH1* cosThetaCHist = new TH1D("cosThetaCHist", "cosThetaCHist; #phi(radians); events per", 100, -1,1); 
  TH2* cosThetaAcosThetaCHist = new TH2D("cosThetaAcosThetaCHist", "cosThetaAcosThetaCHist",100,-1,1,100,-1,1); 

  

  TCanvas* quad1CanvasAB = new TCanvas("quad1CanvasAB", "quad1CanvasAB", 10,10,800,600); 
  TCanvas* quad2CanvasAB = new TCanvas("quad2CanvasAB", "quad2CanvasAB", 10,10,800,600); 
  TCanvas* quad3CanvasAB = new TCanvas("quad3CanvasAB", "quad3CanvasAB", 10,10,800,600); 
  TCanvas* quad4CanvasAB = new TCanvas("quad4CanvasAB", "quad4CanvasAB", 10,10,800,600); 
  TH1* quad1HistAB = new TH1D("quad1: D0" + DecayOne, "D0 Quadrant 1: Resonance Pairs: " + DecayOne + "; #phi(radians); events per 2#pi/100", 100,0,6.28); 
  TH1* quad2HistAB = new TH1D("quad2: D0" + DecayOne, "D0 Quadrant 2: Resonance Pairs: " + DecayOne + "; #phi(radians); events per 2#pi/100", 100,0,6.28); 
  TH1* quad3HistAB = new TH1D("quad3: D0" + DecayOne, "D0 Quadrant 3: Resonance Pairs: " + DecayOne + "; #phi(radians); events per 2#pi/100", 100,0,6.28); 
  TH1* quad4HistAB = new TH1D("quad4: D0" + DecayOne, "D0 Quadrant 4: Resonance Pairs: " + DecayOne + "; #phi(radians); events per 2#pi/100", 100,0,6.28); 


  TCanvas* quad1Phi1CanvasAB = new TCanvas("quad1Phi1CanvasAB", "quad1Phi1CanvasAB", 10,10,800,600); 
  TCanvas* quad2Phi1CanvasAB = new TCanvas("quad2Phi1CanvasAB", "quad2Phi1CanvasAB", 10,10,800,600); 
  TCanvas* quad3Phi1CanvasAB = new TCanvas("quad3Phi1CanvasAB", "quad3Phi1CanvasAB", 10,10,800,600); 
  TCanvas* quad4Phi1CanvasAB = new TCanvas("quad4Phi1CanvasAB", "quad4Phi1CanvasAB", 10,10,800,600); 
  TH1* quad1Phi1HistAB = new TH1D("quad1: D0" + DecayOne, "D0 Quadrant 1: Resonance Pairs: " + DecayOne + "; #phi(radians); events per 2#pi/100", 100,0,6.28); 
  TH1* quad2Phi1HistAB = new TH1D("quad2: D0" + DecayOne, "D0 Quadrant 2: Resonance Pairs: " + DecayOne + "; #phi(radians); events per 2#pi/100", 100,0,6.28); 
  TH1* quad3Phi1HistAB = new TH1D("quad3: D0" + DecayOne, "D0 Quadrant 3: Resonance Pairs: " + DecayOne + "; #phi(radians); events per 2#pi/100", 100,0,6.28); 
  TH1* quad4Phi1HistAB = new TH1D("quad4: D0" + DecayOne, "D0 Quadrant 4: Resonance Pairs: " + DecayOne + "; #phi(radians); events per 2#pi/100", 100,0,6.28); 
 
  TCanvas* quad1Phi2CanvasAB = new TCanvas("quad1PHi2CanvasAB", "quad1PHi2CanvasAB", 10,10,800,600); 
  TCanvas* quad2Phi2CanvasAB = new TCanvas("quad2PHi2CanvasAB", "quad2PHi2CanvasAB", 10,10,800,600); 
  TCanvas* quad3Phi2CanvasAB = new TCanvas("quad3PHi2CanvasAB", "quad3PHi2CanvasAB", 10,10,800,600); 
  TCanvas* quad4Phi2CanvasAB = new TCanvas("quad4PHi2CanvasAB", "quad4PHi2CanvasAB", 10,10,800,600); 
  TH1* quad1Phi2HistAB = new TH1D("quad1: D0" + DecayOne, "D0 Quadrant 1: Resonance Pairs: " + DecayOne + "; #phi(radians); events per 2#pi/100", 100,0,6.28); 
  TH1* quad2Phi2HistAB = new TH1D("quad2: D0" + DecayOne, "D0 Quadrant 2: Resonance Pairs: " + DecayOne + "; #phi(radians); events per 2#pi/100", 100,0,6.28); 
  TH1* quad3Phi2HistAB = new TH1D("quad3: D0" + DecayOne, "D0 Quadrant 3: Resonance Pairs: " + DecayOne + "; #phi(radians); events per 2#pi/100", 100,0,6.28); 
  TH1* quad4Phi2HistAB = new TH1D("quad4: D0" + DecayOne, "D0 Quadrant 4: Resonance Pairs: " + DecayOne + "; #phi(radians); events per 2#pi/100", 100,0,6.28); 

  TCanvas* quad1DSCanvasAB = new TCanvas("quad1DSCanvasAB", "quad1DSCanvasAB", 10,10,800,600); 
  TCanvas* quad2DSCanvasAB = new TCanvas("quad2DSCanvasAB", "quad2DSCanvasAB", 10,10,800,600); 
  TCanvas* quad3DSCanvasAB = new TCanvas("quad3DSCanvasAB", "quad3DSCanvasAB", 10,10,800,600); 
  TCanvas* quad4DSCanvasAB = new TCanvas("quad4DSCanvasAB", "quad4DSCanvasAB", 10,10,800,600);
  TH1* quad1DSHistAB = new TH1D("quad1: D0-bar" + DecayOne, "D0-bar Quadrant 1: Resonance Pairs: " + DecayOne + "; #phi(radians); events per 2#pi/100", 100,0,6.28); 
  TH1* quad2DSHistAB = new TH1D("quad2: D0-bar" + DecayOne, "D0-bar Quadrant 2: Resonance Pairs: " + DecayOne + "; #phi(radians); events per 2#pi/100", 100,0,6.28); 
  TH1* quad3DSHistAB = new TH1D("quad3: D0-bar" + DecayOne, "D0-bar Quadrant 3: Resonance Pairs: " + DecayOne + "; #phi(radians); events per 2#pi/100", 100,0,6.28); 
  TH1* quad4DSHistAB = new TH1D("quad4: D0-bar" + DecayOne, "D0-bar Quadrant 4: Resonance Pairs: " + DecayOne + "; #phi(radians); events per 2#pi/100", 100,0,6.28); 

  TCanvas* quad1Phi1DSCanvasAB = new TCanvas("quad1DSPhi1CanvasAB", "quad1DSPhi1CanvasAB", 10,10,800,600); 
  TCanvas* quad2Phi1DSCanvasAB = new TCanvas("quad2DSPhi1CanvasAB", "quad2DSPhi1CanvasAB", 10,10,800,600); 
  TCanvas* quad3Phi1DSCanvasAB = new TCanvas("quad3DSPhi1CanvasAB", "quad3DSPhi1CanvasAB", 10,10,800,600); 
  TCanvas* quad4Phi1DSCanvasAB = new TCanvas("quad4DSPhi1CanvasAB", "quad4DSPhi1CanvasAB", 10,10,800,600);
  TH1* quad1Phi1DSHistAB = new TH1D("quad1: D0-bar" + DecayOne, "D0-bar Quadrant 1: Resonance Pairs: " + DecayOne + "; #phi(radians); events per 2#pi/100", 100,0,6.28); 
  TH1* quad2Phi1DSHistAB = new TH1D("quad2: D0-bar" + DecayOne, "D0-bar Quadrant 2: Resonance Pairs: " + DecayOne + "; #phi(radians); events per 2#pi/100", 100,0,6.28); 
  TH1* quad3Phi1DSHistAB = new TH1D("quad3: D0-bar" + DecayOne, "D0-bar Quadrant 3: Resonance Pairs: " + DecayOne + "; #phi(radians); events per 2#pi/100", 100,0,6.28); 
  TH1* quad4Phi1DSHistAB = new TH1D("quad4: D0-bar" + DecayOne, "D0-bar Quadrant 4: Resonance Pairs: " + DecayOne + "; #phi(radians); events per 2#pi/100", 100,0,6.28);
  
  TCanvas* quad1Phi2DSCanvasAB = new TCanvas("quad1DSPhi2CanvasAB", "quad1DSPhi2CanvasAB", 10,10,800,600); 
  TCanvas* quad2Phi2DSCanvasAB = new TCanvas("quad2DSPhi2CanvasAB", "quad2DSPhi2CanvasAB", 10,10,800,600); 
  TCanvas* quad3Phi2DSCanvasAB = new TCanvas("quad3DSPhi2CanvasAB", "quad3DSPhi2CanvasAB", 10,10,800,600); 
  TCanvas* quad4Phi2DSCanvasAB = new TCanvas("quad4DSPhi2CanvasAB", "quad4DSPhi2CanvasAB", 10,10,800,600);
  TH1* quad1Phi2DSHistAB = new TH1D("quad1: D0-bar" + DecayOne, "D0-bar Quadrant 1: Resonance Pairs: " + DecayOne + "; #phi(radians); events per 2#pi/100", 100,0,6.28); 
  TH1* quad2Phi2DSHistAB = new TH1D("quad2: D0-bar" + DecayOne, "D0-bar Quadrant 2: Resonance Pairs: " + DecayOne + "; #phi(radians); events per 2#pi/100", 100,0,6.28); 
  TH1* quad3Phi2DSHistAB = new TH1D("quad3: D0-bar" + DecayOne, "D0-bar Quadrant 3: Resonance Pairs: " + DecayOne + "; #phi(radians); events per 2#pi/100", 100,0,6.28); 
  TH1* quad4Phi2DSHistAB = new TH1D("quad4: D0-bar" + DecayOne, "D0-bar Quadrant 4: Resonance Pairs: " + DecayOne + "; #phi(radians); events per 2#pi/100", 100,0,6.28);

  TCanvas* quad1CanvasAD = new TCanvas("quad1CanvasAD", "quad1CanvasAD", 10,10,800,600); 
  TCanvas* quad2CanvasAD = new TCanvas("quad2CanvasAD", "quad2CanvasAD", 10,10,800,600); 
  TCanvas* quad3CanvasAD = new TCanvas("quad3CanvasAD", "quad3CanvasAD", 10,10,800,600); 
  TCanvas* quad4CanvasAD = new TCanvas("quad4CanvasAD", "quad4CanvasAD", 10,10,800,600); 
  TH1* quad1HistAD = new TH1D("quad1: D0" + DecayTwo, "D0 Quadrant 1: Resonance Pairs: " + DecayTwo + "; #phi(radians); events per 2#pi/100", 100,0,6.28); 
  TH1* quad2HistAD = new TH1D("quad2: D0" + DecayTwo, "D0 Quadrant 2: Resonance Pairs: " + DecayTwo + "; #phi(radians); events per 2#pi/100", 100,0,6.28); 
  TH1* quad3HistAD = new TH1D("quad3: D0" + DecayTwo, "D0 Quadrant 3: Resonance Pairs: " + DecayTwo + "; #phi(radians); events per 2#pi/100", 100,0,6.28); 
  TH1* quad4HistAD = new TH1D("quad4: D0" + DecayTwo, "D0 Quadrant 4: Resonance Pairs: " + DecayTwo + "; #phi(radians); events per 2#pi/100", 100,0,6.28); 


  
  TCanvas* quad1DSCanvasAD = new TCanvas("quad1DSCanvasAD", "quad1DSCanvasAD", 10,10,800,600); 
  TCanvas* quad2DSCanvasAD = new TCanvas("quad2DSCanvasAD", "quad2DSCanvasAD", 10,10,800,600); 
  TCanvas* quad3DSCanvasAD = new TCanvas("quad3DSCanvasAD", "quad3DSCanvasAD", 10,10,800,600); 
  TCanvas* quad4DSCanvasAD = new TCanvas("quad4DSCanvasAD", "quad4DSCanvasAD", 10,10,800,600); 
  TH1* quad1DSHistAD = new TH1D("quad1: D0-bar" + DecayTwo, "D0-bar Quadrant 1: Resonance Pairs: " + DecayTwo + "; #phi(radians); events per 2#pi/100", 100,0,6.28); 
  TH1* quad2DSHistAD = new TH1D("quad2: D0-bar" + DecayTwo, "D0-bar Quadrant 2: Resonance Pairs: " + DecayTwo + "; #phi(radians); events per 2#pi/100", 100,0,6.28); 
  TH1* quad3DSHistAD = new TH1D("quad3: D0-bar" + DecayTwo, "D0-bar Quadrant 3: Resonance Pairs: " + DecayTwo + ";  #phi(radians); events per 2#pi/100", 100,0,6.28); 
  TH1* quad4DSHistAD = new TH1D("quad4: D0-bar" + DecayTwo, "D0-bar Quadrant 4: Resonance Pairs: " + DecayTwo + "; #phi(radians); events per 2#pi/100", 100,0,6.28); 

  

  TCanvas phiCanvas("phiCanvas", "phiCanvas", 800, 600); 
  TCanvas phiHist2Canvas("phi2Canvas", "phi2Canvas", 800, 600); 
  TH1D phiHist("phiHist", "phiHist", 100, 0, 6.28); 
  TH1D phiHist2Hist("phi2Hist", "phi2Hist", 100, 0,6.28);



  TCanvas* quad1CompareABCanvas = new TCanvas("quad1CompareABCanvas", "quad1CompareABCanvas", 10,10,800,600); 
  TCanvas* quad2CompareABCanvas = new TCanvas("quad2CompareABCanvas", "quad2CompareABCanvas", 10,10,800,600); 
  TCanvas* quad3CompareABCanvas = new TCanvas("quad3CompareABCanvas", "quad3CompareABCanvas", 10,10,800,600); 
  TCanvas* quad4CompareABCanvas = new TCanvas("quad4CompareABCanvas", "quad4CompareABCanvas", 10,10,800,600); 

  TCanvas* quad1Phi1CompareABCanvas = new TCanvas("quad1CompareABCanvas", "quad1CompareABCanvas", 10,10,800,600); 
  TCanvas* quad2Phi1CompareABCanvas = new TCanvas("quad2CompareABCanvas", "quad2CompareABCanvas", 10,10,800,600); 
  TCanvas* quad3Phi1CompareABCanvas = new TCanvas("quad3CompareABCanvas", "quad3CompareABCanvas", 10,10,800,600); 
  TCanvas* quad4Phi1CompareABCanvas = new TCanvas("quad4CompareABCanvas", "quad4CompareABCanvas", 10,10,800,600); 

  TCanvas* quad1Phi2CompareABCanvas = new TCanvas("quad1CompareABCanvas", "quad1CompareABCanvas", 10,10,800,600); 
  TCanvas* quad2Phi2CompareABCanvas = new TCanvas("quad2CompareABCanvas", "quad2CompareABCanvas", 10,10,800,600); 
  TCanvas* quad3Phi2CompareABCanvas = new TCanvas("quad3CompareABCanvas", "quad3CompareABCanvas", 10,10,800,600); 
  TCanvas* quad4Phi2CompareABCanvas = new TCanvas("quad4CompareABCanvas", "quad4CompareABCanvas", 10,10,800,600); 

  TCanvas* quad1CompareADCanvas = new TCanvas("quad1CompareADCanvas", "quad1CompareADCanvas", 10,10,800,600); 
  TCanvas* quad2CompareADCanvas = new TCanvas("quad2CompareADCanvas", "quad2CompareADCanvas", 10,10,800,600); 
  TCanvas* quad3CompareADCanvas = new TCanvas("quad3CompareADCanvas", "quad3CompareADCanvas", 10,10,800,600); 
  TCanvas* quad4CompareADCanvas = new TCanvas("quad4CompareADCanvas", "quad4CompareADCanvas", 10,10,800,600); 


  TCanvas phiRhoIntquad1Canvas("phiRhoIntquad1Canvas", "phiRhoIntquad1Canvas", 800,600); 
  TH2D phiRhoIntquad1Hist("phiRhoIntquad1Hist", "m(KK)^{2} m(#pi#pi)^{2} Quad 1; m(KK)^{2} GeV^{2}; m(#pi#pi)^{2} GeV^{2}", 100,0.9,2.5, 100,0.0,1.0); 
  TCanvas phiRhoIntquad2Canvas("phiRhoIntquad2Canvas", "phiRhoIntquad2Canvas", 800,600); 
  TH2D phiRhoIntquad2Hist("phiRhoIntquad2Hist", "m(KK)^{2} m(#pi#pi)^{2} Quad 2; m(KK)^{2} GeV^{2}; m(#pi#pi)^{2} GeV^{2}", 100,0.9,2.5, 100,0.0,1.0); 
  TCanvas phiRhoIntquad3Canvas("phiRhoIntquad3Canvas", "phiRhoIntquad3Canvas", 800,600); 
  TH2D phiRhoIntquad3Hist("phiRhoIntquad3Hist", "m(KK)^{2} m(#pi#pi)^{2} Quad 3; m(KK)^{2} GeV^{2}; m(#pi#pi)^{2} GeV^{2}", 100,0.9,2.5, 100,0.0,1.0); 
  TCanvas phiRhoIntquad4Canvas("phiRhoIntquad4Canvas", "phiRhoIntquad4Canvas", 800,600); 
  TH2D phiRhoIntquad4Hist("phiRhoIntquad4Hist", "m(KK)^{2} m(#pi#pi)^{2} Quad 4; m(KK)^{2} GeV^{2}; m(#pi#pi)^{2} GeV^{2}", 100,0.9,2.5, 100,0.0,1.0); 
 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Define reader
  TFile U1File("TestData.root");  
  TTreeReader U1reader("t1", &U1File);



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Read in variables 

  TTreeReaderValue<double> m12(U1reader, "m12"); 
  TTreeReaderValue<double> m34(U1reader, "m34");
  TTreeReaderValue<double> c12(U1reader, "c12");
  TTreeReaderValue<double> c34(U1reader, "c34");
  TTreeReaderValue<double> phi(U1reader, "phi");
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Begin per event code

  while(U1reader.Next())
  {
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//(Begin #phi(radians) code that is applicable to pAB and pCD pairs)

    if(m12> 1.02 && m12< 1.06 && m34> 0.1  && m34 < 0.7)
    //if(-1>0)
    {
      mAB2mCD2ResHist->Fill(m12, m34);
      phiHist2Hist.Fill(phi); 
      
      Double_t sin2ThetaA = sqrt(1- (2*c12*c12-1)*(2*c12*c12-1)); 
      Double_t sin2ThetaC = sqrt(1- (2*c34*c34-1)*(2*c34*c34 -1)); 
      int quadrant = 0;
      if (sin2ThetaA < 0 && sin2ThetaC < 0)
      {
        quadrant = 1;
      } 
      else if (sin2ThetaA < 0 && sin2ThetaC > 0)
      {
         quadrant = 2;
      } 
      else if (sin2ThetaA > 0 && sin2ThetaC < 0)
      {
         quadrant = 3;
      }
      else if (sin2ThetaA > 0 && sin2ThetaC > 0) 
      {
        quadrant = 4;
      }

      if(1 > 0)
      {
        if (quadrant == 1) 
        {
          quad1HistAB->Fill(phi);
        } 
          else if (quadrant == 2) 
          { 
            quad2HistAB->Fill(phi);
          } 
          else if (quadrant == 3)
          {
            quad3HistAB->Fill(phi);
          } 
          else if (quadrant == 4) 
          {
            quad4HistAB->Fill(phi);
          }
      }

     // else if(*D0_ID < 0)
     // {
     //   if (quadrant == 1)
     //   {
     //     quad1DSHistAB->Fill(phi);
     //   } 
     //   else if (quadrant == 2)
     //   {
     //     quad2DSHistAB->Fill(phi);
     //   } 
     //   else if (quadrant == 3)
     //   {
     //     quad3DSHistAB->Fill(phi);   
     //   }
     //   else if (quadrant == 4)
     //   {
     //     quad4DSHistAB->Fill(phi);
     //   }
     // }

       
      if(m34 > 0.48)
      {
      if(1 > 0)
      {
        if (quadrant == 1) 
        {
          quad1Phi1HistAB->Fill(phi);
        } 
          else if (quadrant == 2) 
          { 
            quad2Phi1HistAB->Fill(phi);
          } 
          else if (quadrant == 3)
          {
            quad3Phi1HistAB->Fill(phi);
          } 
          else if (quadrant == 4) 
          {
            quad4Phi1HistAB->Fill(phi);
          }
      }

    //  else if(*D0_ID < 0)
    //  {
    //    if (quadrant == 1)
    //    {
    //      quad1Phi1DSHistAB->Fill(phi);
    //    } 
    //    else if (quadrant == 2)
    //    {
    //      quad2Phi1DSHistAB->Fill(phi);
    //    } 
    //    else if (quadrant == 3)
    //    {
    //      quad3Phi1DSHistAB->Fill(phi);   
    //    }
    //    else if (quadrant == 4)
    //    {
    //     quad4Phi1DSHistAB->Fill(phi);
    //    }
    //  }
    }

     else
     {
      if(1 > 0)
      {
        if (quadrant == 1) 
        {
          quad1Phi2HistAB->Fill(phi);
        } 
          else if (quadrant == 2) 
          { 
            quad2Phi2HistAB->Fill(phi);
          } 
          else if (quadrant == 3)
          {
            quad3Phi2HistAB->Fill(phi);
          } 
          else if (quadrant == 4) 
          {
            quad4Phi2HistAB->Fill(phi);
          }
      }

      //else if(*D0_ID < 0)
      //{
      //  if (quadrant == 1)
      //  {
      //    quad1Phi2DSHistAB->Fill(phi);
      //  } 
      //  else if (quadrant == 2)
      //  {
      //   quad2Phi2DSHistAB->Fill(phi);
      //  } 
      //  else if (quadrant == 3)
      //  {
      //    quad3Phi2DSHistAB->Fill(phi);   
      //  }
      //  else if (quadrant == 4)
      //  {
      //    quad4Phi2DSHistAB->Fill(phi);
      //  }
      // }
     }
    
    }

   

  //(End #phi(radians) Code for pAB pCD)
}

  mAB2mCD2ResCanvas->cd(); 
  mAB2mCD2ResHist->GetXaxis()->CenterTitle();
  mAB2mCD2ResHist->GetYaxis()->CenterTitle();
  mAB2mCD2ResHist->Draw("colz"); 
  mAB2mCD2ResCanvas->SaveAs("mAB2mCD2ResHist.png");

  mAD2mBC2ResCanvas->cd(); 
  mAD2mBC2ResHist->GetXaxis()->CenterTitle();
  mAD2mBC2ResHist->GetYaxis()->CenterTitle();
  mAD2mBC2ResHist->Draw("colz"); 
  mAD2mBC2ResCanvas->SaveAs("mAD2mBC2ResHist.png");  

  phiCanvas.cd(); 
  phiHist.Draw(); 
  phiCanvas.SaveAs("phiHist.png"); 
 
  phiHist2Canvas.cd(); 
  phiHist2Hist.Draw(); 
  phiHist2Canvas.SaveAs("phi2Hist.png"); 


  phiRhoIntquad1Canvas.cd(); 
  phiRhoIntquad1Hist.GetXaxis()->CenterTitle();
  phiRhoIntquad1Hist.GetYaxis()->CenterTitle(); 
  phiRhoIntquad1Hist.Draw("colz"); 
  phiRhoIntquad1Canvas.SaveAs("phiRhoIntquad1.png"); 

  phiRhoIntquad2Canvas.cd(); 
  phiRhoIntquad2Hist.GetXaxis()->CenterTitle();
  phiRhoIntquad2Hist.GetYaxis()->CenterTitle();
  phiRhoIntquad2Hist.Draw("colz"); 
  phiRhoIntquad2Canvas.SaveAs("phiRhoIntquad2.png"); 

  phiRhoIntquad3Canvas.cd(); 
  phiRhoIntquad3Hist.GetXaxis()->CenterTitle();
  phiRhoIntquad3Hist.GetYaxis()->CenterTitle();
  phiRhoIntquad3Hist.Draw("colz"); 
  phiRhoIntquad3Canvas.SaveAs("phiRhoIntquad3.png"); 

  phiRhoIntquad4Canvas.cd(); 
  phiRhoIntquad4Hist.GetXaxis()->CenterTitle();
  phiRhoIntquad4Hist.GetYaxis()->CenterTitle();
  phiRhoIntquad4Hist.Draw("colz"); 
  phiRhoIntquad4Canvas.SaveAs("phiRhoIntquad4.png"); 


////////////////////////////////////////////////////////////////////////////////////////////////////////////

  mAB2mCD2Canvas->cd(); 
  mAB2mCD2Hist->Draw("colz"); 
  mAB2mCD2Hist->GetXaxis()->CenterTitle(); 
  mAB2mCD2Hist->GetYaxis()->CenterTitle(); 
  mAB2mCD2Hist->GetXaxis()->SetTitleOffset(1.25); 
  mAB2mCD2Hist->GetYaxis()->SetTitleOffset(1.35);
  mAB2mCD2Canvas->SaveAs("mAB2mCD2.png"); 

  mAB2mCD2LogCanvas->cd();
  mAB2mCD2LogCanvas->SetLogz(); 
  mAB2mCD2LogHist->Draw("colz"); 
  mAB2mCD2LogHist->GetXaxis()->CenterTitle(); 
  mAB2mCD2LogHist->GetYaxis()->CenterTitle(); 
  mAB2mCD2LogHist->GetXaxis()->SetTitleOffset(1.25); 
  mAB2mCD2LogHist->GetYaxis()->SetTitleOffset(1.35);
  mAB2mCD2LogCanvas->SaveAs("mAB2mCD2Log.png");

  mAB2mBC2Canvas->cd(); 
  mAB2mBC2Hist->Draw("colz"); 
  mAB2mBC2Hist->GetXaxis()->CenterTitle(); 
  mAB2mBC2Hist->GetYaxis()->CenterTitle(); 
  mAB2mBC2Hist->GetXaxis()->SetTitleOffset(1.25); 
  mAB2mBC2Hist->GetYaxis()->SetTitleOffset(1.35);
  mAB2mBC2Canvas->SaveAs("mAB2mBC2.png"); 

  mAB2mAD2Canvas->cd(); 
  mAB2mAD2Hist->Draw("colz"); 
  mAB2mAD2Hist->GetXaxis()->CenterTitle(); 
  mAB2mAD2Hist->GetYaxis()->CenterTitle(); 
  mAB2mAD2Hist->GetXaxis()->SetTitleOffset(1.25); 
  mAB2mAD2Hist->GetYaxis()->SetTitleOffset(1.35);
  mAB2mAD2Canvas->SaveAs("mAB2mAD2.png"); 

  mAB2mAD2LogCanvas->cd();
  mAB2mAD2LogCanvas->SetLogz(); 
  mAB2mAD2LogHist->Draw("colz"); 
  mAB2mAD2LogHist->GetXaxis()->CenterTitle(); 
  mAB2mAD2LogHist->GetYaxis()->CenterTitle(); 
  mAB2mAD2LogHist->GetXaxis()->SetTitleOffset(1.25); 
  mAB2mAD2LogHist->GetYaxis()->SetTitleOffset(1.35);
  mAB2mAD2LogCanvas->SaveAs("mAB2mAD2Log.png");

  mCD2mBC2Canvas->cd(); 
  mCD2mBC2Hist->Draw("colz"); 
  mCD2mBC2Hist->GetXaxis()->CenterTitle(); 
  mCD2mBC2Hist->GetYaxis()->CenterTitle(); 
  mCD2mBC2Hist->GetXaxis()->SetTitleOffset(1.25); 
  mCD2mBC2Hist->GetYaxis()->SetTitleOffset(1.35);
  mCD2mBC2Canvas->SaveAs("mCD2mBC2.png"); 

  mCD2mBC2LogCanvas->cd();
  mCD2mBC2LogCanvas->SetLogz(); 
  mCD2mBC2LogHist->Draw("colz"); 
  mCD2mBC2LogHist->GetXaxis()->CenterTitle(); 
  mCD2mBC2LogHist->GetYaxis()->CenterTitle(); 
  mCD2mBC2LogHist->GetXaxis()->SetTitleOffset(1.25); 
  mCD2mBC2LogHist->GetYaxis()->SetTitleOffset(1.35);
  mCD2mBC2LogCanvas->SaveAs("mCD2mBC2Log.png"); 
 

  mCD2mAD2Canvas->cd(); 
  mCD2mAD2Hist->Draw("colz"); 
  mCD2mAD2Hist->GetXaxis()->CenterTitle(); 
  mCD2mAD2Hist->GetYaxis()->CenterTitle(); 
  mCD2mAD2Hist->GetXaxis()->SetTitleOffset(1.25); 
  mCD2mAD2Hist->GetYaxis()->SetTitleOffset(1.35);
  mCD2mAD2Canvas->SaveAs("mCD2mAD2.png"); 

  mCD2mAD2LogCanvas->cd();
  mCD2mAD2LogCanvas->SetLogz(); 
  mCD2mAD2LogHist->Draw("colz"); 
  mCD2mAD2LogHist->GetXaxis()->CenterTitle(); 
  mCD2mAD2LogHist->GetYaxis()->CenterTitle(); 
  mCD2mAD2LogHist->GetXaxis()->SetTitleOffset(1.25); 
  mCD2mAD2LogHist->GetYaxis()->SetTitleOffset(1.35);
  mCD2mAD2LogCanvas->SaveAs("mCD2mAD2Log.png"); 

  mAB2mBC2LogCanvas->cd();
  mAB2mBC2LogCanvas->SetLogz(); 
  mAB2mBC2LogHist->Draw("colz"); 
  mAB2mBC2LogHist->GetXaxis()->CenterTitle(); 
  mAB2mBC2LogHist->GetYaxis()->CenterTitle(); 
  mAB2mBC2LogHist->GetXaxis()->SetTitleOffset(1.25); 
  mAB2mBC2LogHist->GetYaxis()->SetTitleOffset(1.35);
  mAB2mBC2LogCanvas->SaveAs("mAB2mBC2Log.png"); 
 
  mAD2mBC2Canvas->cd(); 
  mAD2mBC2Hist->Draw("colz"); 
  mAD2mBC2Hist->GetXaxis()->CenterTitle(); 
  mAD2mBC2Hist->GetYaxis()->CenterTitle(); 
  mAD2mBC2Hist->GetXaxis()->SetTitleOffset(1.25); 
  mAD2mBC2Hist->GetYaxis()->SetTitleOffset(1.35);
  mAD2mBC2Canvas->SaveAs("mAD2mBC2.png");

  mAD2mBC2LogCanvas->cd(); 
  mAD2mBC2LogCanvas->SetLogz(); 
  mAD2mBC2LogHist->Draw("colz"); 
  mAD2mBC2LogHist->GetXaxis()->CenterTitle(); 
  mAD2mBC2LogHist->GetYaxis()->CenterTitle(); 
  mAD2mBC2LogHist->GetXaxis()->SetTitleOffset(1.25); 
  mAD2mBC2LogHist->GetYaxis()->SetTitleOffset(1.35);
  mAD2mBC2LogCanvas->SaveAs("mAD2mBC2Log.png");

  ////////////////////////////////////////////////////////////////////////////////////////////////////

  mABCanvas->cd(); 
  mABHist->Draw(); 
  mABHist->GetXaxis()->CenterTitle(); 
  mABHist->GetYaxis()->CenterTitle(); 
  mABHist->GetXaxis()->SetTitleOffset(1.25); 
  mABHist->GetYaxis()->SetTitleOffset(1.35);
  mABCanvas->SaveAs("mAB.png"); 

  mADCanvas->cd(); 
  mADHist->Draw(); 
  mADHist->GetXaxis()->CenterTitle(); 
  mADHist->GetYaxis()->CenterTitle(); 
  mADHist->GetXaxis()->SetTitleOffset(1.25); 
  mADHist->GetYaxis()->SetTitleOffset(1.35);
  mADCanvas->SaveAs("mAD.png"); 

  mBCCanvas->cd(); 
  mBCHist->Draw(); 
  mBCHist->GetXaxis()->CenterTitle(); 
  mBCHist->GetYaxis()->CenterTitle(); 
  mBCHist->GetXaxis()->SetTitleOffset(1.25); 
  mBCHist->GetYaxis()->SetTitleOffset(1.35);
  mBCCanvas->SaveAs("mBC.png"); 

  mCDCanvas->cd(); 
  mCDHist->Draw(); 
  mCDHist->GetXaxis()->CenterTitle(); 
  mCDHist->GetYaxis()->CenterTitle(); 
  mCDHist->GetXaxis()->SetTitleOffset(1.25); 
  mCDHist->GetYaxis()->SetTitleOffset(1.35);
  mCDCanvas->SaveAs("mCD.png"); 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  double nD0; 
  double nD0bar; 

  TString quad1AsymmetryABStr; 
  TString quad2AsymmetryABStr; 
  TString quad3AsymmetryABStr; 
  TString quad4AsymmetryABStr; 

  TString quad1DSAsymmetryABStr; 
  TString quad2DSAsymmetryABStr; 
  TString quad3DSAsymmetryABStr; 
  TString quad4DSAsymmetryABStr; 

  TString quad1AsymmetryADStr;  
  TString quad2AsymmetryADStr;  
  TString quad3AsymmetryADStr;  
  TString quad4AsymmetryADStr;  

  TString quad1DSAsymmetryADStr;  
  TString quad2DSAsymmetryADStr;  
  TString quad3DSAsymmetryADStr;  
  TString quad4DSAsymmetryADStr;  

  TString quad1AsymmetryErrABStr; 
  TString quad2AsymmetryErrABStr; 
  TString quad3AsymmetryErrABStr; 
  TString quad4AsymmetryErrABStr; 

  TString quad1DSAsymmetryErrABStr; 
  TString quad2DSAsymmetryErrABStr; 
  TString quad3DSAsymmetryErrABStr; 
  TString quad4DSAsymmetryErrABStr; 

  TString quad1AsymmetryErrADStr;  
  TString quad2AsymmetryErrADStr;  
  TString quad3AsymmetryErrADStr;  
  TString quad4AsymmetryErrADStr;  

  TString quad1DSAsymmetryErrADStr;  
  TString quad2DSAsymmetryErrADStr;  
  TString quad3DSAsymmetryErrADStr;  
  TString quad4DSAsymmetryErrADStr;  

  TString quad1CPViolationABStr; 
  TString quad2CPViolationABStr; 
  TString quad3CPViolationABStr; 
  TString quad4CPViolationABStr; 

  TString quad1DSCPViolationABStr; 
  TString quad2DSCPViolationABStr; 
  TString quad3DSCPViolationABStr; 
  TString quad4DSCPViolationABStr; 

  TString quad1CPViolationErrABStr;
  TString quad2CPViolationErrABStr;
  TString quad3CPViolationErrABStr;
  TString quad4CPViolationErrABStr; 

  TString quad1DSCPViolationErrABStr;
  TString quad2DSCPViolationErrABStr;
  TString quad3DSCPViolationErrABStr;
  TString quad4DSCPViolationErrABStr; 

  TString quad1CPViolationADStr;  
  TString quad2CPViolationADStr;  
  TString quad3CPViolationADStr;  
  TString quad4CPViolationADStr;  

  TString quad1DSCPViolationADStr;  
  TString quad2DSCPViolationADStr;  
  TString quad3DSCPViolationADStr;  
  TString quad4DSCPViolationADStr;

  TString quad1CPViolationErrADStr; 
  TString quad2CPViolationErrADStr; 
  TString quad3CPViolationErrADStr; 
  TString quad4CPViolationErrADStr;  

  TString quad1DSCPViolationErrADStr; 
  TString quad2DSCPViolationErrADStr; 
  TString quad3DSCPViolationErrADStr; 
  TString quad4DSCPViolationErrADStr; 

  TString quad1D0ABAsymmetryStr; 
  TString quad2D0ABAsymmetryStr; 
  TString quad3D0ABAsymmetryStr; 
  TString quad4D0ABAsymmetryStr;

  TString quad1D0ABAsymmetryErrStr; 
  TString quad2D0ABAsymmetryErrStr; 
  TString quad3D0ABAsymmetryErrStr; 
  TString quad4D0ABAsymmetryErrStr; 

  TString quad1D0barABAsymmetryStr; 
  TString quad2D0barABAsymmetryStr; 
  TString quad3D0barABAsymmetryStr; 
  TString quad4D0barABAsymmetryStr;

  TString quad1D0barABAsymmetryErrStr; 
  TString quad2D0barABAsymmetryErrStr; 
  TString quad3D0barABAsymmetryErrStr; 
  TString quad4D0barABAsymmetryErrStr; 

  TString quad1D0ADAsymmetryStr; 
  TString quad2D0ADAsymmetryStr; 
  TString quad3D0ADAsymmetryStr; 
  TString quad4D0ADAsymmetryStr;

  TString quad1D0ADAsymmetryErrStr; 
  TString quad2D0ADAsymmetryErrStr; 
  TString quad3D0ADAsymmetryErrStr; 
  TString quad4D0ADAsymmetryErrStr; 

  TString quad1D0barADAsymmetryStr; 
  TString quad2D0barADAsymmetryStr; 
  TString quad3D0barADAsymmetryStr; 
  TString quad4D0barADAsymmetryStr;

  TString quad1D0barADAsymmetryErrStr; 
  TString quad2D0barADAsymmetryErrStr; 
  TString quad3D0barADAsymmetryErrStr; 
  TString quad4D0barADAsymmetryErrStr; 

  TString quad1AsymmetryAB2PhiStr;  
  TString quad2AsymmetryAB2PhiStr;  
  TString quad3AsymmetryAB2PhiStr;  
  TString quad4AsymmetryAB2PhiStr;  

  TString quad1DSAsymmetryAB2PhiStr;  
  TString quad2DSAsymmetryAB2PhiStr;  
  TString quad3DSAsymmetryAB2PhiStr;  
  TString quad4DSAsymmetryAB2PhiStr;  

  TString quad1AsymmetryAD2PhiStr;   
  TString quad2AsymmetryAD2PhiStr;   
  TString quad3AsymmetryAD2PhiStr;   
  TString quad4AsymmetryAD2PhiStr;   

  TString quad1DSAsymmetryAD2PhiStr;   
  TString quad2DSAsymmetryAD2PhiStr;   
  TString quad3DSAsymmetryAD2PhiStr;   
  TString quad4DSAsymmetryAD2PhiStr;   

  TString quad1AsymmetryErrAB2PhiStr;  
  TString quad2AsymmetryErrAB2PhiStr;  
  TString quad3AsymmetryErrAB2PhiStr;  
  TString quad4AsymmetryErrAB2PhiStr;  

  TString quad1DSAsymmetryErrAB2PhiStr;  
  TString quad2DSAsymmetryErrAB2PhiStr;  
  TString quad3DSAsymmetryErrAB2PhiStr;  
  TString quad4DSAsymmetryErrAB2PhiStr;  

  TString quad1AsymmetryErrAD2PhiStr;   
  TString quad2AsymmetryErrAD2PhiStr;   
  TString quad3AsymmetryErrAD2PhiStr;   
  TString quad4AsymmetryErrAD2PhiStr;   

  TString quad1DSAsymmetryErrAD2PhiStr;   
  TString quad2DSAsymmetryErrAD2PhiStr;   
  TString quad3DSAsymmetryErrAD2PhiStr;   
  TString quad4DSAsymmetryErrAD2PhiStr;   

  TString quad1CPViolationAB2PhiStr;  
  TString quad2CPViolationAB2PhiStr;  
  TString quad3CPViolationAB2PhiStr;  
  TString quad4CPViolationAB2PhiStr;  

  TString quad1DSCPViolationAB2PhiStr;  
  TString quad2DSCPViolationAB2PhiStr;  
  TString quad3DSCPViolationAB2PhiStr;  
  TString quad4DSCPViolationAB2PhiStr;  

  TString quad1CPViolationErrAB2PhiStr; 
  TString quad2CPViolationErrAB2PhiStr; 
  TString quad3CPViolationErrAB2PhiStr; 
  TString quad4CPViolationErrAB2PhiStr;  

  TString quad1DSCPViolationErrAB2PhiStr; 
  TString quad2DSCPViolationErrAB2PhiStr; 
  TString quad3DSCPViolationErrAB2PhiStr; 
  TString quad4DSCPViolationErrAB2PhiStr;  

  TString quad1CPViolationAD2PhiStr;   
  TString quad2CPViolationAD2PhiStr;   
  TString quad3CPViolationAD2PhiStr;   
  TString quad4CPViolationAD2PhiStr;   

  TString quad1DSCPViolationAD2PhiStr;   
  TString quad2DSCPViolationAD2PhiStr;   
  TString quad3DSCPViolationAD2PhiStr;   
  TString quad4DSCPViolationAD2PhiStr; 

  TString quad1CPViolationErrAD2PhiStr;  
  TString quad2CPViolationErrAD2PhiStr;  
  TString quad3CPViolationErrAD2PhiStr;  
  TString quad4CPViolationErrAD2PhiStr;   

  TString quad1DSCPViolationErrAD2PhiStr;  
  TString quad2DSCPViolationErrAD2PhiStr;  
  TString quad3DSCPViolationErrAD2PhiStr;  
  TString quad4DSCPViolationErrAD2PhiStr;  

  TString quad1D0ABAsymmetry2PhiStr;  
  TString quad2D0ABAsymmetry2PhiStr;  
  TString quad3D0ABAsymmetry2PhiStr;  
  TString quad4D0ABAsymmetry2PhiStr; 

  TString quad1D0ABAsymmetryErr2PhiStr;  
  TString quad2D0ABAsymmetryErr2PhiStr;  
  TString quad3D0ABAsymmetryErr2PhiStr;  
  TString quad4D0ABAsymmetryErr2PhiStr;  

  TString quad1D0barABAsymmetry2PhiStr;  
  TString quad2D0barABAsymmetry2PhiStr;  
  TString quad3D0barABAsymmetry2PhiStr;  
  TString quad4D0barABAsymmetry2PhiStr; 

  TString quad1D0barABAsymmetryErr2PhiStr;  
  TString quad2D0barABAsymmetryErr2PhiStr;  
  TString quad3D0barABAsymmetryErr2PhiStr;  
  TString quad4D0barABAsymmetryErr2PhiStr;  

  TString quad1D0ADAsymmetry2PhiStr;  
  TString quad2D0ADAsymmetry2PhiStr;  
  TString quad3D0ADAsymmetry2PhiStr;  
  TString quad4D0ADAsymmetry2PhiStr; 

  TString quad1D0ADAsymmetryErr2PhiStr;  
  TString quad2D0ADAsymmetryErr2PhiStr;  
  TString quad3D0ADAsymmetryErr2PhiStr;  
  TString quad4D0ADAsymmetryErr2PhiStr;  

  TString quad1D0barADAsymmetry2PhiStr;  
  TString quad2D0barADAsymmetry2PhiStr;  
  TString quad3D0barADAsymmetry2PhiStr;  
  TString quad4D0barADAsymmetry2PhiStr; 

  TString quad1D0barADAsymmetryErr2PhiStr;  
  TString quad2D0barADAsymmetryErr2PhiStr;  
  TString quad3D0barADAsymmetryErr2PhiStr;  
  TString quad4D0barADAsymmetryErr2PhiStr; 

  TString quad1Phi1AsymmetryABStr; 
  TString quad2Phi1AsymmetryABStr; 
  TString quad3Phi1AsymmetryABStr; 
  TString quad4Phi1AsymmetryABStr; 

  TString quad1Phi1DSAsymmetryABStr; 
  TString quad2Phi1DSAsymmetryABStr; 
  TString quad3Phi1DSAsymmetryABStr; 
  TString quad4Phi1DSAsymmetryABStr; 

  TString quad1Phi1AsymmetryADStr;  
  TString quad2Phi1AsymmetryADStr;  
  TString quad3Phi1AsymmetryADStr;  
  TString quad4Phi1AsymmetryADStr;  

  TString quad1Phi1DSAsymmetryADStr;  
  TString quad2Phi1DSAsymmetryADStr;  
  TString quad3Phi1DSAsymmetryADStr;  
  TString quad4Phi1DSAsymmetryADStr;  

  TString quad1Phi1AsymmetryErrABStr; 
  TString quad2Phi1AsymmetryErrABStr; 
  TString quad3Phi1AsymmetryErrABStr; 
  TString quad4Phi1AsymmetryErrABStr; 

  TString quad1Phi1DSAsymmetryErrABStr; 
  TString quad2Phi1DSAsymmetryErrABStr; 
  TString quad3Phi1DSAsymmetryErrABStr; 
  TString quad4Phi1DSAsymmetryErrABStr; 

  TString quad1Phi1AsymmetryErrADStr;  
  TString quad2Phi1AsymmetryErrADStr;  
  TString quad3Phi1AsymmetryErrADStr;  
  TString quad4Phi1AsymmetryErrADStr;  

  TString quad1Phi1DSAsymmetryErrADStr;  
  TString quad2Phi1DSAsymmetryErrADStr;  
  TString quad3Phi1DSAsymmetryErrADStr;  
  TString quad4Phi1DSAsymmetryErrADStr;  

  TString quad1Phi1CPViolationABStr; 
  TString quad2Phi1CPViolationABStr; 
  TString quad3Phi1CPViolationABStr; 
  TString quad4Phi1CPViolationABStr; 

  TString quad1Phi1DSCPViolationABStr; 
  TString quad2Phi1DSCPViolationABStr; 
  TString quad3Phi1DSCPViolationABStr; 
  TString quad4Phi1DSCPViolationABStr; 

  TString quad1Phi1CPViolationErrABStr;
  TString quad2Phi1CPViolationErrABStr;
  TString quad3Phi1CPViolationErrABStr;
  TString quad4Phi1CPViolationErrABStr; 

  TString quad1Phi1DSCPViolationErrABStr;
  TString quad2Phi1DSCPViolationErrABStr;
  TString quad3Phi1DSCPViolationErrABStr;
  TString quad4Phi1DSCPViolationErrABStr; 

  TString quad1Phi1CPViolationADStr;  
  TString quad2Phi1CPViolationADStr;  
  TString quad3Phi1CPViolationADStr;  
  TString quad4Phi1CPViolationADStr;  

  TString quad1Phi1DSCPViolationADStr;  
  TString quad2Phi1DSCPViolationADStr;  
  TString quad3Phi1DSCPViolationADStr;  
  TString quad4Phi1DSCPViolationADStr;

  TString quad1Phi1CPViolationErrADStr; 
  TString quad2Phi1CPViolationErrADStr; 
  TString quad3Phi1CPViolationErrADStr; 
  TString quad4Phi1CPViolationErrADStr;  

  TString quad1Phi1DSCPViolationErrADStr; 
  TString quad2Phi1DSCPViolationErrADStr; 
  TString quad3Phi1DSCPViolationErrADStr; 
  TString quad4Phi1DSCPViolationErrADStr; 

  TString quad1Phi1D0ABAsymmetryStr; 
  TString quad2Phi1D0ABAsymmetryStr; 
  TString quad3Phi1D0ABAsymmetryStr; 
  TString quad4Phi1D0ABAsymmetryStr;

  TString quad1Phi1D0ABAsymmetryErrStr; 
  TString quad2Phi1D0ABAsymmetryErrStr; 
  TString quad3Phi1D0ABAsymmetryErrStr; 
  TString quad4Phi1D0ABAsymmetryErrStr; 

  TString quad1Phi1D0barABAsymmetryStr; 
  TString quad2Phi1D0barABAsymmetryStr; 
  TString quad3Phi1D0barABAsymmetryStr; 
  TString quad4Phi1D0barABAsymmetryStr;

  TString quad1Phi1D0barABAsymmetryErrStr; 
  TString quad2Phi1D0barABAsymmetryErrStr; 
  TString quad3Phi1D0barABAsymmetryErrStr; 
  TString quad4Phi1D0barABAsymmetryErrStr; 

  TString quad1Phi1D0ADAsymmetryStr; 
  TString quad2Phi1D0ADAsymmetryStr; 
  TString quad3Phi1D0ADAsymmetryStr; 
  TString quad4Phi1D0ADAsymmetryStr;

  TString quad1Phi1D0ADAsymmetryErrStr; 
  TString quad2Phi1D0ADAsymmetryErrStr; 
  TString quad3Phi1D0ADAsymmetryErrStr; 
  TString quad4Phi1D0ADAsymmetryErrStr; 

  TString quad1Phi1D0barADAsymmetryStr; 
  TString quad2Phi1D0barADAsymmetryStr; 
  TString quad3Phi1D0barADAsymmetryStr; 
  TString quad4Phi1D0barADAsymmetryStr;

  TString quad1Phi1D0barADAsymmetryErrStr; 
  TString quad2Phi1D0barADAsymmetryErrStr; 
  TString quad3Phi1D0barADAsymmetryErrStr; 
  TString quad4Phi1D0barADAsymmetryErrStr; 

  TString quad1Phi1AsymmetryAB2PhiStr;  
  TString quad2Phi1AsymmetryAB2PhiStr;  
  TString quad3Phi1AsymmetryAB2PhiStr;  
  TString quad4Phi1AsymmetryAB2PhiStr;  

  TString quad1Phi1DSAsymmetryAB2PhiStr;  
  TString quad2Phi1DSAsymmetryAB2PhiStr;  
  TString quad3Phi1DSAsymmetryAB2PhiStr;  
  TString quad4Phi1DSAsymmetryAB2PhiStr;  

  TString quad1Phi1AsymmetryAD2PhiStr;   
  TString quad2Phi1AsymmetryAD2PhiStr;   
  TString quad3Phi1AsymmetryAD2PhiStr;   
  TString quad4Phi1AsymmetryAD2PhiStr;   

  TString quad1Phi1DSAsymmetryAD2PhiStr;   
  TString quad2Phi1DSAsymmetryAD2PhiStr;   
  TString quad3Phi1DSAsymmetryAD2PhiStr;   
  TString quad4Phi1DSAsymmetryAD2PhiStr;   

  TString quad1Phi1AsymmetryErrAB2PhiStr;  
  TString quad2Phi1AsymmetryErrAB2PhiStr;  
  TString quad3Phi1AsymmetryErrAB2PhiStr;  
  TString quad4Phi1AsymmetryErrAB2PhiStr;  

  TString quad1Phi1DSAsymmetryErrAB2PhiStr;  
  TString quad2Phi1DSAsymmetryErrAB2PhiStr;  
  TString quad3Phi1DSAsymmetryErrAB2PhiStr;  
  TString quad4Phi1DSAsymmetryErrAB2PhiStr;  

  TString quad1Phi1AsymmetryErrAD2PhiStr;   
  TString quad2Phi1AsymmetryErrAD2PhiStr;   
  TString quad3Phi1AsymmetryErrAD2PhiStr;   
  TString quad4Phi1AsymmetryErrAD2PhiStr;   

  TString quad1Phi1DSAsymmetryErrAD2PhiStr;   
  TString quad2Phi1DSAsymmetryErrAD2PhiStr;   
  TString quad3Phi1DSAsymmetryErrAD2PhiStr;   
  TString quad4Phi1DSAsymmetryErrAD2PhiStr;   

  TString quad1Phi1CPViolationAB2PhiStr;  
  TString quad2Phi1CPViolationAB2PhiStr;  
  TString quad3Phi1CPViolationAB2PhiStr;  
  TString quad4Phi1CPViolationAB2PhiStr;  

  TString quad1Phi1DSCPViolationAB2PhiStr;  
  TString quad2Phi1DSCPViolationAB2PhiStr;  
  TString quad3Phi1DSCPViolationAB2PhiStr;  
  TString quad4Phi1DSCPViolationAB2PhiStr;  

  TString quad1Phi1CPViolationErrAB2PhiStr; 
  TString quad2Phi1CPViolationErrAB2PhiStr; 
  TString quad3Phi1CPViolationErrAB2PhiStr; 
  TString quad4Phi1CPViolationErrAB2PhiStr;  

  TString quad1Phi1DSCPViolationErrAB2PhiStr; 
  TString quad2Phi1DSCPViolationErrAB2PhiStr; 
  TString quad3Phi1DSCPViolationErrAB2PhiStr; 
  TString quad4Phi1DSCPViolationErrAB2PhiStr;  

  TString quad1Phi1CPViolationAD2PhiStr;   
  TString quad2Phi1CPViolationAD2PhiStr;   
  TString quad3Phi1CPViolationAD2PhiStr;   
  TString quad4Phi1CPViolationAD2PhiStr;   

  TString quad1Phi1DSCPViolationAD2PhiStr;   
  TString quad2Phi1DSCPViolationAD2PhiStr;   
  TString quad3Phi1DSCPViolationAD2PhiStr;   
  TString quad4Phi1DSCPViolationAD2PhiStr; 

  TString quad1Phi1CPViolationErrAD2PhiStr;  
  TString quad2Phi1CPViolationErrAD2PhiStr;  
  TString quad3Phi1CPViolationErrAD2PhiStr;  
  TString quad4Phi1CPViolationErrAD2PhiStr;   

  TString quad1Phi1DSCPViolationErrAD2PhiStr;  
  TString quad2Phi1DSCPViolationErrAD2PhiStr;  
  TString quad3Phi1DSCPViolationErrAD2PhiStr;  
  TString quad4Phi1DSCPViolationErrAD2PhiStr;  

  TString quad1Phi1D0ABAsymmetry2PhiStr;  
  TString quad2Phi1D0ABAsymmetry2PhiStr;  
  TString quad3Phi1D0ABAsymmetry2PhiStr;  
  TString quad4Phi1D0ABAsymmetry2PhiStr; 

  TString quad1Phi1D0ABAsymmetryErr2PhiStr;  
  TString quad2Phi1D0ABAsymmetryErr2PhiStr;  
  TString quad3Phi1D0ABAsymmetryErr2PhiStr;  
  TString quad4Phi1D0ABAsymmetryErr2PhiStr;  

  TString quad1Phi1D0barABAsymmetry2PhiStr;  
  TString quad2Phi1D0barABAsymmetry2PhiStr;  
  TString quad3Phi1D0barABAsymmetry2PhiStr;  
  TString quad4Phi1D0barABAsymmetry2PhiStr; 

  TString quad1Phi1D0barABAsymmetryErr2PhiStr;  
  TString quad2Phi1D0barABAsymmetryErr2PhiStr;  
  TString quad3Phi1D0barABAsymmetryErr2PhiStr;  
  TString quad4Phi1D0barABAsymmetryErr2PhiStr;  

  TString quad1Phi1D0ADAsymmetry2PhiStr;  
  TString quad2Phi1D0ADAsymmetry2PhiStr;  
  TString quad3Phi1D0ADAsymmetry2PhiStr;  
  TString quad4Phi1D0ADAsymmetry2PhiStr; 

  TString quad1Phi1D0ADAsymmetryErr2PhiStr;  
  TString quad2Phi1D0ADAsymmetryErr2PhiStr;  
  TString quad3Phi1D0ADAsymmetryErr2PhiStr;  
  TString quad4Phi1D0ADAsymmetryErr2PhiStr;  

  TString quad1Phi1D0barADAsymmetry2PhiStr;  
  TString quad2Phi1D0barADAsymmetry2PhiStr;  
  TString quad3Phi1D0barADAsymmetry2PhiStr;  
  TString quad4Phi1D0barADAsymmetry2PhiStr; 

  TString quad1Phi1D0barADAsymmetryErr2PhiStr;  
  TString quad2Phi1D0barADAsymmetryErr2PhiStr;  
  TString quad3Phi1D0barADAsymmetryErr2PhiStr;  
  TString quad4Phi1D0barADAsymmetryErr2PhiStr;

  TString quad1Phi2AsymmetryABStr; 
  TString quad2Phi2AsymmetryABStr; 
  TString quad3Phi2AsymmetryABStr; 
  TString quad4Phi2AsymmetryABStr; 

  TString quad1Phi2DSAsymmetryABStr; 
  TString quad2Phi2DSAsymmetryABStr; 
  TString quad3Phi2DSAsymmetryABStr; 
  TString quad4Phi2DSAsymmetryABStr; 

  TString quad1Phi2AsymmetryADStr;  
  TString quad2Phi2AsymmetryADStr;  
  TString quad3Phi2AsymmetryADStr;  
  TString quad4Phi2AsymmetryADStr;  

  TString quad1Phi2DSAsymmetryADStr;  
  TString quad2Phi2DSAsymmetryADStr;  
  TString quad3Phi2DSAsymmetryADStr;  
  TString quad4Phi2DSAsymmetryADStr;  

  TString quad1Phi2AsymmetryErrABStr; 
  TString quad2Phi2AsymmetryErrABStr; 
  TString quad3Phi2AsymmetryErrABStr; 
  TString quad4Phi2AsymmetryErrABStr; 

  TString quad1Phi2DSAsymmetryErrABStr; 
  TString quad2Phi2DSAsymmetryErrABStr; 
  TString quad3Phi2DSAsymmetryErrABStr; 
  TString quad4Phi2DSAsymmetryErrABStr; 

  TString quad1Phi2AsymmetryErrADStr;  
  TString quad2Phi2AsymmetryErrADStr;  
  TString quad3Phi2AsymmetryErrADStr;  
  TString quad4Phi2AsymmetryErrADStr;  

  TString quad1Phi2DSAsymmetryErrADStr;  
  TString quad2Phi2DSAsymmetryErrADStr;  
  TString quad3Phi2DSAsymmetryErrADStr;  
  TString quad4Phi2DSAsymmetryErrADStr;  

  TString quad1Phi2CPViolationABStr; 
  TString quad2Phi2CPViolationABStr; 
  TString quad3Phi2CPViolationABStr; 
  TString quad4Phi2CPViolationABStr; 

  TString quad1Phi2DSCPViolationABStr; 
  TString quad2Phi2DSCPViolationABStr; 
  TString quad3Phi2DSCPViolationABStr; 
  TString quad4Phi2DSCPViolationABStr; 

  TString quad1Phi2CPViolationErrABStr;
  TString quad2Phi2CPViolationErrABStr;
  TString quad3Phi2CPViolationErrABStr;
  TString quad4Phi2CPViolationErrABStr; 

  TString quad1Phi2DSCPViolationErrABStr;
  TString quad2Phi2DSCPViolationErrABStr;
  TString quad3Phi2DSCPViolationErrABStr;
  TString quad4Phi2DSCPViolationErrABStr; 

  TString quad1Phi2CPViolationADStr;  
  TString quad2Phi2CPViolationADStr;  
  TString quad3Phi2CPViolationADStr;  
  TString quad4Phi2CPViolationADStr;  

  TString quad1Phi2DSCPViolationADStr;  
  TString quad2Phi2DSCPViolationADStr;  
  TString quad3Phi2DSCPViolationADStr;  
  TString quad4Phi2DSCPViolationADStr;

  TString quad1Phi2CPViolationErrADStr; 
  TString quad2Phi2CPViolationErrADStr; 
  TString quad3Phi2CPViolationErrADStr; 
  TString quad4Phi2CPViolationErrADStr;  

  TString quad1Phi2DSCPViolationErrADStr; 
  TString quad2Phi2DSCPViolationErrADStr; 
  TString quad3Phi2DSCPViolationErrADStr; 
  TString quad4Phi2DSCPViolationErrADStr; 

  TString quad1Phi2D0ABAsymmetryStr; 
  TString quad2Phi2D0ABAsymmetryStr; 
  TString quad3Phi2D0ABAsymmetryStr; 
  TString quad4Phi2D0ABAsymmetryStr;

  TString quad1Phi2D0ABAsymmetryErrStr; 
  TString quad2Phi2D0ABAsymmetryErrStr; 
  TString quad3Phi2D0ABAsymmetryErrStr; 
  TString quad4Phi2D0ABAsymmetryErrStr; 

  TString quad1Phi2D0barABAsymmetryStr; 
  TString quad2Phi2D0barABAsymmetryStr; 
  TString quad3Phi2D0barABAsymmetryStr; 
  TString quad4Phi2D0barABAsymmetryStr;

  TString quad1Phi2D0barABAsymmetryErrStr; 
  TString quad2Phi2D0barABAsymmetryErrStr; 
  TString quad3Phi2D0barABAsymmetryErrStr; 
  TString quad4Phi2D0barABAsymmetryErrStr; 

  TString quad1Phi2D0ADAsymmetryStr; 
  TString quad2Phi2D0ADAsymmetryStr; 
  TString quad3Phi2D0ADAsymmetryStr; 
  TString quad4Phi2D0ADAsymmetryStr;

  TString quad1Phi2D0ADAsymmetryErrStr; 
  TString quad2Phi2D0ADAsymmetryErrStr; 
  TString quad3Phi2D0ADAsymmetryErrStr; 
  TString quad4Phi2D0ADAsymmetryErrStr; 

  TString quad1Phi2D0barADAsymmetryStr; 
  TString quad2Phi2D0barADAsymmetryStr; 
  TString quad3Phi2D0barADAsymmetryStr; 
  TString quad4Phi2D0barADAsymmetryStr;

  TString quad1Phi2D0barADAsymmetryErrStr; 
  TString quad2Phi2D0barADAsymmetryErrStr; 
  TString quad3Phi2D0barADAsymmetryErrStr; 
  TString quad4Phi2D0barADAsymmetryErrStr; 

  TString quad1Phi2AsymmetryAB2PhiStr;  
  TString quad2Phi2AsymmetryAB2PhiStr;  
  TString quad3Phi2AsymmetryAB2PhiStr;  
  TString quad4Phi2AsymmetryAB2PhiStr;  

  TString quad1Phi2DSAsymmetryAB2PhiStr;  
  TString quad2Phi2DSAsymmetryAB2PhiStr;  
  TString quad3Phi2DSAsymmetryAB2PhiStr;  
  TString quad4Phi2DSAsymmetryAB2PhiStr;  

  TString quad1Phi2AsymmetryAD2PhiStr;   
  TString quad2Phi2AsymmetryAD2PhiStr;   
  TString quad3Phi2AsymmetryAD2PhiStr;   
  TString quad4Phi2AsymmetryAD2PhiStr;   

  TString quad1Phi2DSAsymmetryAD2PhiStr;   
  TString quad2Phi2DSAsymmetryAD2PhiStr;   
  TString quad3Phi2DSAsymmetryAD2PhiStr;   
  TString quad4Phi2DSAsymmetryAD2PhiStr;   

  TString quad1Phi2AsymmetryErrAB2PhiStr;  
  TString quad2Phi2AsymmetryErrAB2PhiStr;  
  TString quad3Phi2AsymmetryErrAB2PhiStr;  
  TString quad4Phi2AsymmetryErrAB2PhiStr;  

  TString quad1Phi2DSAsymmetryErrAB2PhiStr;  
  TString quad2Phi2DSAsymmetryErrAB2PhiStr;  
  TString quad3Phi2DSAsymmetryErrAB2PhiStr;  
  TString quad4Phi2DSAsymmetryErrAB2PhiStr;  

  TString quad1Phi2AsymmetryErrAD2PhiStr;   
  TString quad2Phi2AsymmetryErrAD2PhiStr;   
  TString quad3Phi2AsymmetryErrAD2PhiStr;   
  TString quad4Phi2AsymmetryErrAD2PhiStr;   

  TString quad1Phi2DSAsymmetryErrAD2PhiStr;   
  TString quad2Phi2DSAsymmetryErrAD2PhiStr;   
  TString quad3Phi2DSAsymmetryErrAD2PhiStr;   
  TString quad4Phi2DSAsymmetryErrAD2PhiStr;   

  TString quad1Phi2CPViolationAB2PhiStr;  
  TString quad2Phi2CPViolationAB2PhiStr;  
  TString quad3Phi2CPViolationAB2PhiStr;  
  TString quad4Phi2CPViolationAB2PhiStr;  

  TString quad1Phi2DSCPViolationAB2PhiStr;  
  TString quad2Phi2DSCPViolationAB2PhiStr;  
  TString quad3Phi2DSCPViolationAB2PhiStr;  
  TString quad4Phi2DSCPViolationAB2PhiStr;  

  TString quad1Phi2CPViolationErrAB2PhiStr; 
  TString quad2Phi2CPViolationErrAB2PhiStr; 
  TString quad3Phi2CPViolationErrAB2PhiStr; 
  TString quad4Phi2CPViolationErrAB2PhiStr;  

  TString quad1Phi2DSCPViolationErrAB2PhiStr; 
  TString quad2Phi2DSCPViolationErrAB2PhiStr; 
  TString quad3Phi2DSCPViolationErrAB2PhiStr; 
  TString quad4Phi2DSCPViolationErrAB2PhiStr;  

  TString quad1Phi2CPViolationAD2PhiStr;   
  TString quad2Phi2CPViolationAD2PhiStr;   
  TString quad3Phi2CPViolationAD2PhiStr;   
  TString quad4Phi2CPViolationAD2PhiStr;   

  TString quad1Phi2DSCPViolationAD2PhiStr;   
  TString quad2Phi2DSCPViolationAD2PhiStr;   
  TString quad3Phi2DSCPViolationAD2PhiStr;   
  TString quad4Phi2DSCPViolationAD2PhiStr; 

  TString quad1Phi2CPViolationErrAD2PhiStr;  
  TString quad2Phi2CPViolationErrAD2PhiStr;  
  TString quad3Phi2CPViolationErrAD2PhiStr;  
  TString quad4Phi2CPViolationErrAD2PhiStr;   

  TString quad1Phi2DSCPViolationErrAD2PhiStr;  
  TString quad2Phi2DSCPViolationErrAD2PhiStr;  
  TString quad3Phi2DSCPViolationErrAD2PhiStr;  
  TString quad4Phi2DSCPViolationErrAD2PhiStr;  

  TString quad1Phi2D0ABAsymmetry2PhiStr;  
  TString quad2Phi2D0ABAsymmetry2PhiStr;  
  TString quad3Phi2D0ABAsymmetry2PhiStr;  
  TString quad4Phi2D0ABAsymmetry2PhiStr; 

  TString quad1Phi2D0ABAsymmetryErr2PhiStr;  
  TString quad2Phi2D0ABAsymmetryErr2PhiStr;  
  TString quad3Phi2D0ABAsymmetryErr2PhiStr;  
  TString quad4Phi2D0ABAsymmetryErr2PhiStr;  

  TString quad1Phi2D0barABAsymmetry2PhiStr;  
  TString quad2Phi2D0barABAsymmetry2PhiStr;  
  TString quad3Phi2D0barABAsymmetry2PhiStr;  
  TString quad4Phi2D0barABAsymmetry2PhiStr; 

  TString quad1Phi2D0barABAsymmetryErr2PhiStr;  
  TString quad2Phi2D0barABAsymmetryErr2PhiStr;  
  TString quad3Phi2D0barABAsymmetryErr2PhiStr;  
  TString quad4Phi2D0barABAsymmetryErr2PhiStr;  

  TString quad1Phi2D0ADAsymmetry2PhiStr;  
  TString quad2Phi2D0ADAsymmetry2PhiStr;  
  TString quad3Phi2D0ADAsymmetry2PhiStr;  
  TString quad4Phi2D0ADAsymmetry2PhiStr; 

  TString quad1Phi2D0ADAsymmetryErr2PhiStr;  
  TString quad2Phi2D0ADAsymmetryErr2PhiStr;  
  TString quad3Phi2D0ADAsymmetryErr2PhiStr;  
  TString quad4Phi2D0ADAsymmetryErr2PhiStr;  

  TString quad1Phi2D0barADAsymmetry2PhiStr;  
  TString quad2Phi2D0barADAsymmetry2PhiStr;  
  TString quad3Phi2D0barADAsymmetry2PhiStr;  
  TString quad4Phi2D0barADAsymmetry2PhiStr; 

  TString quad1Phi2D0barADAsymmetryErr2PhiStr;  
  TString quad2Phi2D0barADAsymmetryErr2PhiStr;  
  TString quad3Phi2D0barADAsymmetryErr2PhiStr;  
  TString quad4Phi2D0barADAsymmetryErr2PhiStr;

  Double_t D0Sin2Phi_asymmetry; 
  Double_t D0Sin2Phi_asymmetryErr; 
  Double_t D0barSin2Phi_asymmetry; 
  Double_t D0barSin2Phi_asymmetryErr; 
  Double_t Sin2PhiCP_asymmetry; 
  Double_t Sin2PhiCP_asymmetryErr; 
  Double_t DiffSin2Phi; 


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 if(1>0)
{
  nD0 = quad1HistAB->GetEntries();
  nD0bar = quad1DSHistAB->GetEntries();
  Double_t scale = nD0/nD0bar;
  Double_t nD0Above = 0;
  Double_t nD0barAbove = 0;

  for (int bin=0; bin<51; bin++) 
  {
    nD0Above = nD0Above+quad1HistAB->GetBinContent(bin);
    nD0barAbove = nD0barAbove+quad1DSHistAB->GetBinContent(bin);
  }
   
  Double_t nD0Below = 0;
  Double_t nD0barBelow = 0;
  for (int bin=51; bin<101; bin++)
  {
    nD0Below = nD0Below+quad1HistAB->GetBinContent(bin);
    nD0barBelow = nD0barBelow+quad1DSHistAB->GetBinContent(bin);
  }

  Double_t D0AB_asymmetry = (nD0Above-nD0Below)/(nD0Above+nD0Below);
  Double_t D0AB_asymmetryErr = sqrt( (1.0 - D0AB_asymmetry*D0AB_asymmetry)/(nD0Above+nD0Below) );

  Double_t D0barAB_asymmetry = (nD0barAbove-nD0barBelow)/(nD0barAbove+nD0barBelow);
  Double_t D0barAB_asymmetryErr = sqrt( (1.0 - D0barAB_asymmetry*D0barAB_asymmetry)/(nD0barAbove+nD0barBelow) );

  quad1D0ABAsymmetryStr.Form("%6.2f\n", D0AB_asymmetry);
  quad1D0ABAsymmetryErrStr.Form("%6.2f\n", D0AB_asymmetryErr); 

  quad1D0barABAsymmetryStr.Form("%6.2f\n", D0barAB_asymmetry); 
  quad1D0barABAsymmetryErrStr.Form("%6.2f\n", D0barAB_asymmetryErr); 

  Double_t weight1 = 1./(D0AB_asymmetryErr*D0AB_asymmetryErr);
  Double_t weight2 = 1./(D0barAB_asymmetryErr*D0barAB_asymmetryErr);

  Double_t CPV = (D0AB_asymmetry*weight1+D0barAB_asymmetry*weight2)/(weight1+weight2);
  Double_t CPVErr = sqrt(1./(weight1+weight2));

  quad1AsymmetryABStr.Form("%6.4f\n",CPV); 
  quad1AsymmetryErrABStr.Form("%6.4f\n", CPVErr); 

  Double_t nD02PhiAbove = 0; 
  Double_t nD0bar2PhiAbove = 0;
  Double_t nD02PhiBelow = 0; 
  Double_t nD0bar2PhiBelow = 0;  

  for(int bin = 1; bin < 26; bin++)
  {
    nD02PhiAbove = nD02PhiAbove+quad1HistAB->GetBinContent(bin); 
    nD0bar2PhiAbove = nD0bar2PhiAbove+quad1DSHistAB->GetBinContent(bin); 
  }

  for(int bin = 26; bin < 51; bin++)
  {
    nD0bar2PhiBelow = nD0bar2PhiBelow + quad1DSHistAB->GetBinContent(bin); 
    nD02PhiBelow = nD02PhiBelow + quad1HistAB->GetBinContent(bin); 
  }

  for(int bin = 51; bin < 76; bin++)
  {
    nD02PhiAbove = nD02PhiAbove+quad1HistAB->GetBinContent(bin); 
    nD0bar2PhiAbove = nD0bar2PhiAbove+quad1DSHistAB->GetBinContent(bin); 
  }

  for(int bin = 76; bin <101; bin++)
  {
    nD0bar2PhiBelow = nD0bar2PhiBelow + quad1DSHistAB->GetBinContent(bin); 
    nD02PhiBelow = nD02PhiBelow + quad1HistAB->GetBinContent(bin); 
  }

    

  D0Sin2Phi_asymmetry = (nD02PhiAbove - nD02PhiBelow)/(nD02PhiAbove+nD02PhiBelow); 
  D0Sin2Phi_asymmetryErr = sqrt((1-D0Sin2Phi_asymmetry*D0Sin2Phi_asymmetry)/( nD02PhiAbove + nD02PhiBelow )); 

  D0barSin2Phi_asymmetry = (nD0bar2PhiAbove - nD0bar2PhiBelow)/(nD0bar2PhiAbove+nD0bar2PhiBelow); 
  D0barSin2Phi_asymmetryErr = sqrt((1-D0barSin2Phi_asymmetry*D0barSin2Phi_asymmetry)/( nD0bar2PhiAbove + nD0bar2PhiBelow ));

  weight1 = 1./(D0Sin2Phi_asymmetry*D0Sin2Phi_asymmetry); 
  weight2 = 1./(D0barSin2Phi_asymmetry*D0barSin2Phi_asymmetry); 

  Sin2PhiCP_asymmetry = (D0Sin2Phi_asymmetry*weight1+D0barSin2Phi_asymmetry*weight2)/(weight1+weight2);
  Sin2PhiCP_asymmetryErr = sqrt(1./(weight1+weight2));
  DiffSin2Phi = (D0Sin2Phi_asymmetry*weight1-D0barSin2Phi_asymmetry*weight2)/(weight1+weight2);

  TString D0Sin2Phi_asymmetryStr1; 
  TString D0Sin2Phi_asymmetryErrStr1; 
  TString D0barSin2Phi_asymmetryStr1; 
  TString D0barSin2Phi_asymmetryErrStr1; 
  TString Sin2PhiCP_asymmetryStr1; 
  TString Sin2PhiCP_asymmetryErrStr1; 
  TString DiffSin2PhiStr1; 

  D0Sin2Phi_asymmetryStr1.Form("%6.2f\n", D0Sin2Phi_asymmetry); 
  D0Sin2Phi_asymmetryErrStr1.Form("%6.2f\n", D0Sin2Phi_asymmetryErr); 
  D0barSin2Phi_asymmetryStr1.Form("%6.2f\n", D0barSin2Phi_asymmetry);
  D0barSin2Phi_asymmetryErrStr1.Form("%6.2f\n", D0barSin2Phi_asymmetryErr); 
  Sin2PhiCP_asymmetryStr1.Form("%6.4f\n", Sin2PhiCP_asymmetry); 
  Sin2PhiCP_asymmetryErrStr1.Form("%6.4f\n", Sin2PhiCP_asymmetryErr); 
  DiffSin2PhiStr1.Form("%6.2f\n", DiffSin2Phi); 

      

  quad1CanvasAB->cd(); 
  quad1HistAB->GetXaxis()->CenterTitle(); 
  quad1HistAB->GetYaxis()->CenterTitle(); 
  quad1HistAB->GetXaxis()->SetTitleOffset(1.25); 
  quad1HistAB->GetYaxis()->SetTitleOffset(1.35);
  quad1HistAB->Draw(); 
  quad1CanvasAB->SaveAs("quad1HistAB.png"); 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  

    nD0 = quad2HistAB->GetEntries();
    nD0bar = quad2DSHistAB->GetEntries();
      
    scale = nD0/nD0bar;
    nD0Above = 0;
    nD0barAbove = 0;

    for (int bin=0; bin<51; bin++) 
    {
        nD0Above = nD0Above+quad2HistAB->GetBinContent(bin);
        nD0barAbove = nD0barAbove+quad2DSHistAB->GetBinContent(bin);
    }
   
    nD0Below = 0;
    nD0barBelow = 0;
    for (int bin=51; bin<101; bin++)
    {
        nD0Below = nD0Below+quad2HistAB->GetBinContent(bin);
        nD0barBelow = nD0barBelow+quad2DSHistAB->GetBinContent(bin);
    }

    D0AB_asymmetry = (nD0Above-nD0Below)/(nD0Above+nD0Below);
    D0AB_asymmetryErr = sqrt( (1.0 - D0AB_asymmetry*D0AB_asymmetry)/(nD0Above+nD0Below) );

    D0barAB_asymmetry = (nD0barAbove-nD0barBelow)/(nD0barAbove+nD0barBelow);
    D0barAB_asymmetryErr = sqrt( (1.0 - D0barAB_asymmetry*D0barAB_asymmetry)/(nD0barAbove+nD0barBelow) );

    quad2D0ABAsymmetryStr.Form("%6.2f\n", D0AB_asymmetry);
    quad2D0ABAsymmetryErrStr.Form("%6.2f\n", D0AB_asymmetryErr); 

    quad2D0barABAsymmetryStr.Form("%6.2f\n", D0barAB_asymmetry); 
    quad2D0barABAsymmetryErrStr.Form("%6.2f\n", D0barAB_asymmetryErr); 

    //std::cout<< D0AB_asymmetry <<std::endl; 
    //std::cout<<D0barAB_asymmetry <<std::endl; 

    weight1 = 1./(D0AB_asymmetryErr*D0AB_asymmetryErr);
    weight2 = 1./(D0barAB_asymmetryErr*D0barAB_asymmetryErr);

  CPV = (D0AB_asymmetry*weight1+D0barAB_asymmetry*weight2)/(weight1+weight2);
  CPVErr = sqrt(1./(weight1+weight2));

  //std::cout<<"CPA " << CPV << std::endl; 

  quad2AsymmetryABStr.Form("%6.4f\n",CPV); 
  quad2AsymmetryErrABStr.Form("%6.4f\n", CPVErr); 

  for(int bin = 1; bin < 26; bin++)
    {
        nD02PhiAbove = nD02PhiAbove+quad2HistAB->GetBinContent(bin); 
        nD0bar2PhiAbove = nD0bar2PhiAbove+quad2DSHistAB->GetBinContent(bin); 
    }

    for(int bin = 26; bin < 51; bin++)
    {
        nD0bar2PhiBelow = nD0bar2PhiBelow + quad2DSHistAB->GetBinContent(bin); 
        nD02PhiBelow = nD02PhiBelow + quad2HistAB->GetBinContent(bin); 
    }

    for(int bin = 51; bin < 76; bin++)
    {
        nD02PhiAbove = nD02PhiAbove+quad2HistAB->GetBinContent(bin); 
        nD0bar2PhiAbove = nD0bar2PhiAbove+quad2DSHistAB->GetBinContent(bin); 
    }

    for(int bin = 76; bin <101; bin++)
    {
        nD0bar2PhiBelow = nD0bar2PhiBelow + quad2DSHistAB->GetBinContent(bin); 
        nD02PhiBelow = nD02PhiBelow + quad2HistAB->GetBinContent(bin); 
    }

    D0Sin2Phi_asymmetry = (nD02PhiAbove - nD02PhiBelow)/(nD02PhiAbove+nD02PhiBelow); 
    D0Sin2Phi_asymmetryErr = sqrt((1-D0Sin2Phi_asymmetry*D0Sin2Phi_asymmetry)/( nD02PhiAbove + nD02PhiBelow )); 

    D0barSin2Phi_asymmetry = (nD0bar2PhiAbove - nD0bar2PhiBelow)/(nD0bar2PhiAbove+nD0bar2PhiBelow); 
    D0barSin2Phi_asymmetryErr = sqrt((1-D0barSin2Phi_asymmetry*D0barSin2Phi_asymmetry)/( nD0bar2PhiAbove + nD0bar2PhiBelow ));

    weight1 = 1./(D0Sin2Phi_asymmetry*D0Sin2Phi_asymmetry); 
    weight2 = 1./(D0barSin2Phi_asymmetry*D0barSin2Phi_asymmetry); 

    Sin2PhiCP_asymmetry = (D0Sin2Phi_asymmetry*weight1+D0barSin2Phi_asymmetry*weight2)/(weight1+weight2);
    Sin2PhiCP_asymmetryErr = sqrt(1./(weight1+weight2));

    DiffSin2Phi = (D0Sin2Phi_asymmetry*weight1-D0barSin2Phi_asymmetry*weight2)/(weight1+weight2);

    TString D0Sin2Phi_asymmetryStr2; 
    TString D0Sin2Phi_asymmetryErrStr2; 
    TString D0barSin2Phi_asymmetryStr2; 
    TString D0barSin2Phi_asymmetryErrStr2; 
    TString Sin2PhiCP_asymmetryStr2; 
    TString Sin2PhiCP_asymmetryErrStr2; 
    TString DiffSin2PhiStr2;

    D0Sin2Phi_asymmetryStr2.Form("%6.2f\n", D0Sin2Phi_asymmetry); 
    D0Sin2Phi_asymmetryErrStr2.Form("%6.2f\n", D0Sin2Phi_asymmetryErr); 
    D0barSin2Phi_asymmetryStr2.Form("%6.2f\n", D0barSin2Phi_asymmetry);
    D0barSin2Phi_asymmetryErrStr2.Form("%6.2f\n", D0barSin2Phi_asymmetryErr); 
    Sin2PhiCP_asymmetryStr2.Form("%6.4f\n", Sin2PhiCP_asymmetry); 
    Sin2PhiCP_asymmetryErrStr2.Form("%6.4f\n", Sin2PhiCP_asymmetryErr); 
    DiffSin2PhiStr2.Form("%6.2f\n", DiffSin2Phi); 

  quad2CanvasAB->cd(); 
  quad2HistAB->GetXaxis()->CenterTitle(); 
  quad2HistAB->GetYaxis()->CenterTitle(); 
  quad2HistAB->GetXaxis()->SetTitleOffset(1.25); 
  quad2HistAB->GetYaxis()->SetTitleOffset(1.35);
  quad2HistAB->Draw(); 
  quad2CanvasAB->SaveAs("quad2HistAB.png"); 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    nD0 = quad3HistAB->GetEntries();
    nD0bar = quad3DSHistAB->GetEntries();
      
    scale = nD0/nD0bar;
    nD0Above = 0;
    nD0barAbove = 0;

    for (int bin=0; bin<51; bin++) 
    {
      nD0Above = nD0Above+quad3HistAB->GetBinContent(bin);
      nD0barAbove = nD0barAbove+quad3DSHistAB->GetBinContent(bin);
  }
   
  nD0Below = 0;
    nD0barBelow = 0;
    for (int bin=51; bin<101; bin++)
    {
        nD0Below = nD0Below+quad3HistAB->GetBinContent(bin);
      nD0barBelow = nD0barBelow+quad3DSHistAB->GetBinContent(bin);
    }

    D0AB_asymmetry = (nD0Above-nD0Below)/(nD0Above+nD0Below);
    D0AB_asymmetryErr = sqrt( (1.0 - D0AB_asymmetry*D0AB_asymmetry)/(nD0Above+nD0Below) );

    D0barAB_asymmetry = (nD0barAbove-nD0barBelow)/(nD0barAbove+nD0barBelow);
    D0barAB_asymmetryErr = sqrt( (1.0 - D0barAB_asymmetry*D0barAB_asymmetry)/(nD0barAbove+nD0barBelow) );

    quad3D0ABAsymmetryStr.Form("%6.2f\n", D0AB_asymmetry);
    quad3D0ABAsymmetryErrStr.Form("%6.2f\n", D0AB_asymmetryErr); 

    quad3D0barABAsymmetryStr.Form("%6.2f\n", D0barAB_asymmetry); 
    quad3D0barABAsymmetryErrStr.Form("%6.2f\n", D0barAB_asymmetryErr); 

    //std::cout<< D0AB_asymmetry <<std::endl; 
    //std::cout<<D0barAB_asymmetry <<std::endl; 

  weight1 = 1./(D0AB_asymmetryErr*D0AB_asymmetryErr);
    weight2 = 1./(D0barAB_asymmetryErr*D0barAB_asymmetryErr);

  CPV = (D0AB_asymmetry*weight1+D0barAB_asymmetry*weight2)/(weight1+weight2);
  CPVErr = sqrt(1./(weight1+weight2));

  quad3AsymmetryABStr.Form("%6.4f\n",CPV); 
  quad3AsymmetryErrABStr.Form("%6.4f\n", CPVErr); 

  //std::cout<<"CPA " << CPV << std::endl; 

  for(int bin = 1; bin < 26; bin++)
    {
        nD02PhiAbove = nD02PhiAbove+quad3HistAB->GetBinContent(bin); 
        nD0bar2PhiAbove = nD0bar2PhiAbove+quad3DSHistAB->GetBinContent(bin); 
    }

    for(int bin = 26; bin < 51; bin++)
    {
        nD0bar2PhiBelow = nD0bar2PhiBelow + quad3DSHistAB->GetBinContent(bin); 
        nD02PhiBelow = nD02PhiBelow + quad3HistAB->GetBinContent(bin); 
    }

    for(int bin = 51; bin < 76; bin++)
    {
        nD02PhiAbove = nD02PhiAbove+quad3HistAB->GetBinContent(bin); 
        nD0bar2PhiAbove = nD0bar2PhiAbove+quad3DSHistAB->GetBinContent(bin); 
    }

    for(int bin = 76; bin <101; bin++)
    {
        nD0bar2PhiBelow = nD0bar2PhiBelow + quad3DSHistAB->GetBinContent(bin); 
        nD02PhiBelow = nD02PhiBelow + quad3HistAB->GetBinContent(bin); 
    }

    D0Sin2Phi_asymmetry = (nD02PhiAbove - nD02PhiBelow)/(nD02PhiAbove+nD02PhiBelow); 
    D0Sin2Phi_asymmetryErr = sqrt((1-D0Sin2Phi_asymmetry*D0Sin2Phi_asymmetry)/( nD02PhiAbove + nD02PhiBelow )); 

    D0barSin2Phi_asymmetry = (nD0bar2PhiAbove - nD0bar2PhiBelow)/(nD0bar2PhiAbove+nD0bar2PhiBelow); 
    D0barSin2Phi_asymmetryErr = sqrt((1-D0barSin2Phi_asymmetry*D0barSin2Phi_asymmetry)/( nD0bar2PhiAbove + nD0bar2PhiBelow ));

    weight1 = 1./(D0Sin2Phi_asymmetry*D0Sin2Phi_asymmetry); 
    weight2 = 1./(D0barSin2Phi_asymmetry*D0barSin2Phi_asymmetry); 

    Sin2PhiCP_asymmetry = (D0Sin2Phi_asymmetry*weight1+D0barSin2Phi_asymmetry*weight2)/(weight1+weight2);
    Sin2PhiCP_asymmetryErr = sqrt(1./(weight1+weight2));

    DiffSin2Phi = (D0Sin2Phi_asymmetry*weight1-D0barSin2Phi_asymmetry*weight2)/(weight1+weight2);

    TString D0Sin2Phi_asymmetryStr3; 
    TString D0Sin2Phi_asymmetryErrStr3; 
    TString D0barSin2Phi_asymmetryStr3; 
    TString D0barSin2Phi_asymmetryErrStr3; 
    TString Sin2PhiCP_asymmetryStr3; 
    TString Sin2PhiCP_asymmetryErrStr3; 
    TString DiffSin2PhiStr3;

    D0Sin2Phi_asymmetryStr3.Form("%6.2f\n", D0Sin2Phi_asymmetry); 
    D0Sin2Phi_asymmetryErrStr3.Form("%6.2f\n", D0Sin2Phi_asymmetryErr); 
    D0barSin2Phi_asymmetryStr3.Form("%6.2f\n", D0barSin2Phi_asymmetry);
    D0barSin2Phi_asymmetryErrStr3.Form("%6.2f\n", D0barSin2Phi_asymmetryErr); 
    Sin2PhiCP_asymmetryStr3.Form("%6.4f\n", Sin2PhiCP_asymmetry); 
    Sin2PhiCP_asymmetryErrStr3.Form("%6.4f\n", Sin2PhiCP_asymmetryErr); 
    DiffSin2PhiStr3.Form("%6.2f\n", DiffSin2Phi); 

  quad3CanvasAB->cd(); 
  quad3HistAB->GetXaxis()->CenterTitle(); 
  quad3HistAB->GetYaxis()->CenterTitle(); 
  quad3HistAB->GetXaxis()->SetTitleOffset(1.25); 
  quad3HistAB->GetYaxis()->SetTitleOffset(1.35);
  quad3HistAB->Draw(); 
  quad3CanvasAB->SaveAs("quad3HistAB.png"); 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  nD0 = quad4HistAB->GetEntries();
    nD0bar = quad4DSHistAB->GetEntries();
  scale = nD0/nD0bar;
    nD0Above = 0;
    nD0barAbove = 0;

    for (int bin=0; bin<51; bin++) 
    {
      nD0Above = nD0Above+quad4HistAB->GetBinContent(bin);
      nD0barAbove = nD0barAbove+quad4DSHistAB->GetBinContent(bin);
  }
   
  nD0Below = 0;
    nD0barBelow = 0;
    for (int bin=51; bin<101; bin++)
    {
      nD0Below = nD0Below+quad4HistAB->GetBinContent(bin);
      nD0barBelow = nD0barBelow+quad4DSHistAB->GetBinContent(bin);
    }

    D0AB_asymmetry = (nD0Above-nD0Below)/(nD0Above+nD0Below);
    D0AB_asymmetryErr = sqrt( (1.0 - D0AB_asymmetry*D0AB_asymmetry)/(nD0Above+nD0Below) );

    D0barAB_asymmetry = (nD0barAbove-nD0barBelow)/(nD0barAbove+nD0barBelow);
    D0barAB_asymmetryErr = sqrt( (1.0 - D0barAB_asymmetry*D0barAB_asymmetry)/(nD0barAbove+nD0barBelow) );

    quad4D0ABAsymmetryStr.Form("%6.2f\n", D0AB_asymmetry);
    quad4D0ABAsymmetryErrStr.Form("%6.2f\n", D0AB_asymmetryErr); 

    quad4D0barABAsymmetryStr.Form("%6.2f\n", D0barAB_asymmetry); 
    quad4D0barABAsymmetryErrStr.Form("%6.2f\n", D0barAB_asymmetryErr); 

    //std::cout<< D0AB_asymmetry <<std::endl; 
    //std::cout<<D0barAB_asymmetry <<std::endl; 

  weight1 = 1./(D0AB_asymmetryErr*D0AB_asymmetryErr);
    weight2 = 1./(D0barAB_asymmetryErr*D0barAB_asymmetryErr);

  CPV = (D0AB_asymmetry*weight1+D0barAB_asymmetry*weight2)/(weight1+weight2);
  CPVErr = sqrt(1./(weight1+weight2));

  quad4AsymmetryABStr.Form("%6.4f\n",CPV); 
  quad4AsymmetryErrABStr.Form("%6.4f\n", CPVErr); 

  //std::cout<<"CPA " << CPV << std::endl; 

  for(int bin = 1; bin < 26; bin++)
    {
        nD02PhiAbove = nD02PhiAbove+quad4HistAB->GetBinContent(bin); 
        nD0bar2PhiAbove = nD0bar2PhiAbove+quad4DSHistAB->GetBinContent(bin); 
    }

    for(int bin = 26; bin < 51; bin++)
    {
        nD0bar2PhiBelow = nD0bar2PhiBelow + quad4DSHistAB->GetBinContent(bin); 
        nD02PhiBelow = nD02PhiBelow + quad4HistAB->GetBinContent(bin); 
    }

    for(int bin = 51; bin < 76; bin++)
    {
        nD02PhiAbove = nD02PhiAbove+quad4HistAB->GetBinContent(bin); 
        nD0bar2PhiAbove = nD0bar2PhiAbove+quad4DSHistAB->GetBinContent(bin); 
    }

    for(int bin = 76; bin <101; bin++)
    {
        nD0bar2PhiBelow = nD0bar2PhiBelow + quad4DSHistAB->GetBinContent(bin); 
        nD02PhiBelow = nD02PhiBelow + quad4HistAB->GetBinContent(bin); 
    }

    D0Sin2Phi_asymmetry = (nD02PhiAbove - nD02PhiBelow)/(nD02PhiAbove+nD02PhiBelow); 
    D0Sin2Phi_asymmetryErr = sqrt((1-D0Sin2Phi_asymmetry*D0Sin2Phi_asymmetry)/( nD02PhiAbove + nD02PhiBelow )); 

    D0barSin2Phi_asymmetry = (nD0bar2PhiAbove - nD0bar2PhiBelow)/(nD0bar2PhiAbove+nD0bar2PhiBelow); 
    D0barSin2Phi_asymmetryErr = sqrt((1-D0barSin2Phi_asymmetry*D0barSin2Phi_asymmetry)/( nD0bar2PhiAbove + nD0bar2PhiBelow ));

    weight1 = 1./(D0Sin2Phi_asymmetry*D0Sin2Phi_asymmetry); 
    weight2 = 1./(D0barSin2Phi_asymmetry*D0barSin2Phi_asymmetry); 

    Sin2PhiCP_asymmetry = (D0Sin2Phi_asymmetry*weight1+D0barSin2Phi_asymmetry*weight2)/(weight1+weight2);
    Sin2PhiCP_asymmetryErr = sqrt(1./(weight1+weight2));

    DiffSin2Phi = (D0Sin2Phi_asymmetry*weight1-D0barSin2Phi_asymmetry*weight2)/(weight1+weight2);

    TString D0Sin2Phi_asymmetryStr4; 
    TString D0Sin2Phi_asymmetryErrStr4; 
    TString D0barSin2Phi_asymmetryStr4; 
    TString D0barSin2Phi_asymmetryErrStr4; 
    TString Sin2PhiCP_asymmetryStr4; 
    TString Sin2PhiCP_asymmetryErrStr4; 
    TString DiffSin2PhiStr4;

    D0Sin2Phi_asymmetryStr4.Form("%6.2f\n", D0Sin2Phi_asymmetry); 
    D0Sin2Phi_asymmetryErrStr4.Form("%6.2f\n", D0Sin2Phi_asymmetryErr); 
    D0barSin2Phi_asymmetryStr4.Form("%6.2f\n", D0barSin2Phi_asymmetry);
    D0barSin2Phi_asymmetryErrStr4.Form("%6.2f\n", D0barSin2Phi_asymmetryErr); 
    Sin2PhiCP_asymmetryStr4.Form("%6.4f\n", Sin2PhiCP_asymmetry); 
    Sin2PhiCP_asymmetryErrStr4.Form("%6.4f\n", Sin2PhiCP_asymmetryErr); 
    DiffSin2PhiStr4.Form("%6.2f\n", DiffSin2Phi); 

  quad4CanvasAB->cd(); 
  quad4HistAB->GetXaxis()->CenterTitle(); 
  quad4HistAB->GetYaxis()->CenterTitle(); 
  quad4HistAB->GetXaxis()->SetTitleOffset(1.25); 
  quad4HistAB->GetYaxis()->SetTitleOffset(1.35);
  quad4HistAB->Draw(); 
  quad4CanvasAB->SaveAs("quad4HistAB.png"); 

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  quad1DSCanvasAB->cd(); 
  quad1DSHistAB->GetXaxis()->CenterTitle(); 
  quad1DSHistAB->GetYaxis()->CenterTitle(); 
  quad1DSHistAB->GetXaxis()->SetTitleOffset(1.25); 
  quad1DSHistAB->GetYaxis()->SetTitleOffset(1.35);
  quad1DSHistAB->Draw(); 
  quad1DSCanvasAB->SaveAs("quad1DSHistAB.png"); 

  

  quad2DSCanvasAB->cd(); 
  quad2DSHistAB->GetXaxis()->CenterTitle(); 
  quad2DSHistAB->GetYaxis()->CenterTitle(); 
  quad2DSHistAB->GetXaxis()->SetTitleOffset(1.25); 
  quad2DSHistAB->GetYaxis()->SetTitleOffset(1.35);
  quad2DSHistAB->Draw(); 
  quad2DSCanvasAB->SaveAs("quad2DSHistAB.png"); 

  

  quad3DSCanvasAB->cd();
  quad3DSHistAB->GetXaxis()->CenterTitle(); 
  quad3DSHistAB->GetYaxis()->CenterTitle(); 
  quad3DSHistAB->GetXaxis()->SetTitleOffset(1.25); 
  quad3DSHistAB->GetYaxis()->SetTitleOffset(1.35); 
  quad3DSHistAB->Draw(); 
  quad3DSCanvasAB->SaveAs("quad3DSHistAB.png"); 

  

  quad4DSCanvasAB->cd(); 
  quad4DSHistAB->GetXaxis()->CenterTitle(); 
  quad4DSHistAB->GetYaxis()->CenterTitle(); 
  quad4DSHistAB->GetXaxis()->SetTitleOffset(1.25); 
  quad4DSHistAB->GetYaxis()->SetTitleOffset(1.35);
  quad4DSHistAB->Draw(); 
  quad4DSCanvasAB->SaveAs("quad4DSHistAB.png"); 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double maxBin; 
  double binContent; 
  TLegend *leg = new TLegend(0.15,0.85,0.3,0.95);
  TPaveText *pt1 = new TPaveText(0.7,0.75,0.98,0.98,"NDC"); // NDC sets coords
  pt1->SetTextSize(0.025); 

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  quad1CompareABCanvas->cd(); 
  quad1HistAB->SetLineColor(kBlue); 
  quad1DSHistAB->SetLineColor(kRed); 
  quad1HistAB->GetYaxis()->SetLimits(0,1650);
  quad1HistAB->SetMinimum(0);
  maxBin = quad1HistAB->GetMaximumBin(); 
  binContent = 1.3*quad1HistAB->GetBinContent(maxBin); 
  quad1HistAB->SetMaximum(binContent); 
  TLine *line1A = new TLine(3.14,0,3.14, binContent); 
  leg->SetHeader("Legend"); 
  leg->AddEntry(quad1HistAB, "Dzero"); 
  leg->AddEntry(quad1DSHistAB, "Dzero-bar"); 
  pt1->AddText("PV D0:" + quad1D0ABAsymmetryStr + " #pm " + quad1D0ABAsymmetryErrStr);
  pt1->AddText("PV D0bar: " + quad1D0barABAsymmetryStr + " #pm " + quad1D0barABAsymmetryErrStr); 
  pt1->AddText("CPV: " + quad1AsymmetryABStr + " #pm" + quad1AsymmetryErrABStr);
  pt1->AddText("PV D0 Sin(2#phi): " + D0Sin2Phi_asymmetryStr1 + " #pm " + D0Sin2Phi_asymmetryErrStr1); 
  pt1->AddText("PV D0bar Sin(2#phi): " + D0barSin2Phi_asymmetryStr1 + " #pm " + D0barSin2Phi_asymmetryErrStr1); 
  pt1->AddText("CPV Sin(2#phi): " + Sin2PhiCP_asymmetryStr1 + " #pm " + Sin2PhiCP_asymmetryErrStr1); 
  pt1->AddText("D0-D0bar : " + DiffSin2PhiStr1); 
  quad1HistAB->SetTitle("Quad 1: " + DecayOne);
  quad1HistAB->Draw(); 
  quad1DSHistAB->Draw("same");
  line1A->Draw("same");
  leg->Draw("same"); 
  pt1->Draw("same");  
  quad1CompareABCanvas->SaveAs("Quad1CompAB.png"); 
  pt1->Clear(); 
  leg->Clear(); 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////   

  quad2CompareABCanvas->cd(); 
  quad2HistAB->SetLineColor(kBlue); 
  quad2DSHistAB->SetLineColor(kRed); 
  quad2HistAB->GetYaxis()->SetLimits(0,1650);
  quad2HistAB->SetMinimum(0);
  maxBin = quad2HistAB->GetMaximumBin(); 
  binContent = 1.3*quad2HistAB->GetBinContent(maxBin); 
  quad2HistAB->SetMaximum(binContent); 
  TLine *line2A = new TLine(3.14,0,3.14, binContent); 
  leg->SetHeader("Legend"); 
  leg->AddEntry(quad2HistAB, "Dzero"); 
  leg->AddEntry(quad2DSHistAB, "Dzero-bar"); 
  pt1->AddText("PV D0:" + quad2D0ABAsymmetryStr + " #pm " + quad2D0ABAsymmetryErrStr);
  pt1->AddText("PV D0bar: " + quad2D0barABAsymmetryStr + " #pm " + quad2D0barABAsymmetryErrStr); 
  pt1->AddText("CPV: " + quad2AsymmetryABStr + " #pm" + quad2AsymmetryErrABStr);
  pt1->AddText("PV D0 Sin(2#phi): " + D0Sin2Phi_asymmetryStr2 + " #pm " + D0Sin2Phi_asymmetryErrStr2); 
  pt1->AddText("PV D0bar Sin(2#phi): " + D0barSin2Phi_asymmetryStr2 + " #pm " + D0barSin2Phi_asymmetryErrStr2); 
  pt1->AddText("CPV Sin(2#phi): " + Sin2PhiCP_asymmetryStr2 + " #pm " + Sin2PhiCP_asymmetryErrStr2); 
  pt1->AddText("D0-D0bar : " + DiffSin2PhiStr2); 
  quad2HistAB->SetTitle("Quad 2: " + DecayOne);
  quad2HistAB->Draw(); 
  quad2DSHistAB->Draw("same");
  line2A->Draw("same");
  leg->Draw("same"); 
  pt1->Draw("same");  
  quad2CompareABCanvas->SaveAs("Quad2CompAB.png"); 
  pt1->Clear(); 
  leg->Clear(); 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  

  quad3CompareABCanvas->cd(); 
  quad3HistAB->SetLineColor(kBlue); 
  quad3DSHistAB->SetLineColor(kRed); 
  quad3HistAB->GetYaxis()->SetLimits(0,1650);
  quad3HistAB->SetMinimum(0);
  maxBin = quad3HistAB->GetMaximumBin(); 
  binContent = 1.3*quad3HistAB->GetBinContent(maxBin); 
  quad3HistAB->SetMaximum(binContent); 
  leg->SetHeader("Legend"); 
  leg->AddEntry(quad3HistAB, "Dzero"); 
  leg->AddEntry(quad3DSHistAB, "Dzero-bar"); 
  pt1->AddText("PV D0:" + quad3D0ABAsymmetryStr + " #pm " + quad3D0ABAsymmetryErrStr);
  pt1->AddText("PV D0bar: " + quad3D0barABAsymmetryStr + " #pm " + quad3D0barABAsymmetryErrStr); 
  pt1->AddText("CPV: " + quad3AsymmetryABStr + " #pm" + quad3AsymmetryErrABStr);
  pt1->AddText("PV D0 Sin(2#phi): " + D0Sin2Phi_asymmetryStr3 + " #pm " + D0Sin2Phi_asymmetryErrStr3); 
  pt1->AddText("PV D0bar Sin(2#phi): " + D0barSin2Phi_asymmetryStr3 + " #pm " + D0barSin2Phi_asymmetryErrStr3); 
  pt1->AddText("CPV Sin(2#phi): " + Sin2PhiCP_asymmetryStr3 + " #pm " + Sin2PhiCP_asymmetryErrStr3); 
  pt1->AddText("D0-D0bar : " + DiffSin2PhiStr3); 
  TLine *line3A = new TLine(3.14,0,3.14, binContent); 
  quad3HistAB->SetTitle("Quad 3: " + DecayOne);
  quad3HistAB->Draw(); 
  quad3DSHistAB->Draw("same");
  line3A->Draw("same");
  leg->Draw("same"); 
  pt1->Draw("same");  
  quad3CompareABCanvas->SaveAs("Quad3CompAB.png"); 
  leg->Clear("same"); 
  pt1->Clear("same");  
  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  quad4CompareABCanvas->cd(); 
  quad4HistAB->SetLineColor(kBlue); 
  quad4DSHistAB->SetLineColor(kRed); 
  quad4HistAB->GetYaxis()->SetLimits(0,1650);
  quad4HistAB->SetMinimum(0);
  maxBin = quad4HistAB->GetMaximumBin(); 
  binContent = 1.3*quad4HistAB->GetBinContent(maxBin); 
  quad4HistAB->SetMaximum(binContent); 
  TLine *line4A = new TLine(3.14,0,3.14, binContent); 
  leg->SetHeader("Legend"); 
  leg->AddEntry(quad4HistAB, "Dzero"); 
  leg->AddEntry(quad4DSHistAB, "Dzero-bar"); 
  pt1->AddText("PV D0:" + quad4D0ABAsymmetryStr + " #pm " + quad4D0ABAsymmetryErrStr);
  pt1->AddText("PV D0bar: " + quad4D0barABAsymmetryStr + " #pm " + quad4D0barABAsymmetryErrStr); 
  pt1->AddText("CPV: " + quad4AsymmetryABStr + " #pm" + quad4AsymmetryErrABStr);
  pt1->AddText("PV D0 Sin(2#phi): " + D0Sin2Phi_asymmetryStr4 + " #pm " + D0Sin2Phi_asymmetryErrStr4); 
  pt1->AddText("PV D0bar Sin(2#phi): " + D0barSin2Phi_asymmetryStr4 + " #pm " + D0barSin2Phi_asymmetryErrStr4); 
  pt1->AddText("CPV Sin(2#phi): " + Sin2PhiCP_asymmetryStr4 + " #pm " + Sin2PhiCP_asymmetryErrStr4); 
  pt1->AddText("D0-D0bar : " + DiffSin2PhiStr4); 
  quad4HistAB->SetTitle("Quad 4: " + DecayOne);
  quad4HistAB->Draw(); 
  quad4DSHistAB->Draw("same");
  line4A->Draw("same");
  leg->Draw("same"); 
  pt1->Draw("same");  
  quad4CompareABCanvas->SaveAs("Quad4CompAB.png"); 
  leg->Clear("same"); 
  pt1->Clear("same");  

}
  /////////////////////////////////PLZPLZPLZ///////////////////////////////////////////////////////////////////////////////////////
 
if(1>0)
{

  nD0 = quad1Phi1HistAB->GetEntries();
  nD0bar = quad1Phi1DSHistAB->GetEntries();
  Double_t scale = nD0/nD0bar;
  Double_t nD0Above = 0;
  Double_t nD0barAbove = 0;

  for (int bin=0; bin<51; bin++) 
  {
    nD0Above = nD0Above+quad1Phi1HistAB->GetBinContent(bin);
    nD0barAbove = nD0barAbove+quad1Phi1DSHistAB->GetBinContent(bin);
  }
   
  Double_t nD0Below = 0;
  Double_t nD0barBelow = 0;
  for (int bin=51; bin<101; bin++)
  {
    nD0Below = nD0Below+quad1Phi1HistAB->GetBinContent(bin);
    nD0barBelow = nD0barBelow+quad1Phi1DSHistAB->GetBinContent(bin);
  }

  Double_t D0AB_asymmetry = (nD0Above-nD0Below)/(nD0Above+nD0Below);
  Double_t D0AB_asymmetryErr = sqrt( (1.0 - D0AB_asymmetry*D0AB_asymmetry)/(nD0Above+nD0Below) );

  Double_t D0barAB_asymmetry = (nD0barAbove-nD0barBelow)/(nD0barAbove+nD0barBelow);
  Double_t D0barAB_asymmetryErr = sqrt( (1.0 - D0barAB_asymmetry*D0barAB_asymmetry)/(nD0barAbove+nD0barBelow) );

  quad1Phi1D0ABAsymmetryStr.Form("%6.2f\n", D0AB_asymmetry);
  quad1Phi1D0ABAsymmetryErrStr.Form("%6.2f\n", D0AB_asymmetryErr); 

  quad1Phi1D0barABAsymmetryStr.Form("%6.2f\n", D0barAB_asymmetry); 
  quad1Phi1D0barABAsymmetryErrStr.Form("%6.2f\n", D0barAB_asymmetryErr); 

  Double_t weight1 = 1./(D0AB_asymmetryErr*D0AB_asymmetryErr);
  Double_t weight2 = 1./(D0barAB_asymmetryErr*D0barAB_asymmetryErr);

  Double_t CPV = (D0AB_asymmetry*weight1+D0barAB_asymmetry*weight2)/(weight1+weight2);
  Double_t CPVErr = sqrt(1./(weight1+weight2));

  quad1Phi1AsymmetryABStr.Form("%6.4f\n",CPV); 
  quad1Phi1AsymmetryErrABStr.Form("%6.4f\n", CPVErr); 

  Double_t nD02PhiAbove = 0; 
  Double_t nD0bar2PhiAbove = 0;
  Double_t nD02PhiBelow = 0; 
  Double_t nD0bar2PhiBelow = 0;  

  for(int bin = 1; bin < 26; bin++)
  {
    nD02PhiAbove = nD02PhiAbove+quad1Phi1HistAB->GetBinContent(bin); 
    nD0bar2PhiAbove = nD0bar2PhiAbove+quad1Phi1DSHistAB->GetBinContent(bin); 
  }

  for(int bin = 26; bin < 51; bin++)
  {
    nD0bar2PhiBelow = nD0bar2PhiBelow + quad1Phi1DSHistAB->GetBinContent(bin); 
    nD02PhiBelow = nD02PhiBelow + quad1Phi1HistAB->GetBinContent(bin); 
  }

  for(int bin = 51; bin < 76; bin++)
  {
    nD02PhiAbove = nD02PhiAbove+quad1Phi1HistAB->GetBinContent(bin); 
    nD0bar2PhiAbove = nD0bar2PhiAbove+quad1Phi1DSHistAB->GetBinContent(bin); 
  }

  for(int bin = 76; bin <101; bin++)
  {
    nD0bar2PhiBelow = nD0bar2PhiBelow + quad1Phi1DSHistAB->GetBinContent(bin); 
    nD02PhiBelow = nD02PhiBelow + quad1Phi1HistAB->GetBinContent(bin); 
  }

    

  D0Sin2Phi_asymmetry = (nD02PhiAbove - nD02PhiBelow)/(nD02PhiAbove+nD02PhiBelow); 
  D0Sin2Phi_asymmetryErr = sqrt((1-D0Sin2Phi_asymmetry*D0Sin2Phi_asymmetry)/( nD02PhiAbove + nD02PhiBelow )); 

  D0barSin2Phi_asymmetry = (nD0bar2PhiAbove - nD0bar2PhiBelow)/(nD0bar2PhiAbove+nD0bar2PhiBelow); 
  D0barSin2Phi_asymmetryErr = sqrt((1-D0barSin2Phi_asymmetry*D0barSin2Phi_asymmetry)/( nD0bar2PhiAbove + nD0bar2PhiBelow ));

  weight1 = 1./(D0Sin2Phi_asymmetry*D0Sin2Phi_asymmetry); 
  weight2 = 1./(D0barSin2Phi_asymmetry*D0barSin2Phi_asymmetry); 

  Sin2PhiCP_asymmetry = (D0Sin2Phi_asymmetry*weight1+D0barSin2Phi_asymmetry*weight2)/(weight1+weight2);
  Sin2PhiCP_asymmetryErr = sqrt(1./(weight1+weight2));
  DiffSin2Phi = (D0Sin2Phi_asymmetry*weight1-D0barSin2Phi_asymmetry*weight2)/(weight1+weight2);

  TString D0Sin2Phi_asymmetryStr1; 
  TString D0Sin2Phi_asymmetryErrStr1; 
  TString D0barSin2Phi_asymmetryStr1; 
  TString D0barSin2Phi_asymmetryErrStr1; 
  TString Sin2PhiCP_asymmetryStr1; 
  TString Sin2PhiCP_asymmetryErrStr1; 
  TString DiffSin2PhiStr1; 

  D0Sin2Phi_asymmetryStr1.Form("%6.2f\n", D0Sin2Phi_asymmetry); 
  D0Sin2Phi_asymmetryErrStr1.Form("%6.2f\n", D0Sin2Phi_asymmetryErr); 
  D0barSin2Phi_asymmetryStr1.Form("%6.2f\n", D0barSin2Phi_asymmetry);
  D0barSin2Phi_asymmetryErrStr1.Form("%6.2f\n", D0barSin2Phi_asymmetryErr); 
  Sin2PhiCP_asymmetryStr1.Form("%6.4f\n", Sin2PhiCP_asymmetry); 
  Sin2PhiCP_asymmetryErrStr1.Form("%6.4f\n", Sin2PhiCP_asymmetryErr); 
  DiffSin2PhiStr1.Form("%6.2f\n", DiffSin2Phi); 

      

  quad1Phi1CanvasAB->cd(); 
  quad1Phi1HistAB->GetXaxis()->CenterTitle(); 
  quad1Phi1HistAB->GetYaxis()->CenterTitle(); 
  quad1Phi1HistAB->GetXaxis()->SetTitleOffset(1.25); 
  quad1Phi1HistAB->GetYaxis()->SetTitleOffset(1.35);
  quad1Phi1HistAB->Draw(); 
  quad1Phi1CanvasAB->SaveAs("quad1Phi1HistAB.png"); 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  

    nD0 = quad2Phi1HistAB->GetEntries();
    nD0bar = quad2Phi1DSHistAB->GetEntries();
      
    scale = nD0/nD0bar;
    nD0Above = 0;
    nD0barAbove = 0;

    for (int bin=0; bin<51; bin++) 
    {
        nD0Above = nD0Above+quad2Phi1HistAB->GetBinContent(bin);
        nD0barAbove = nD0barAbove+quad2Phi1DSHistAB->GetBinContent(bin);
    }
   
    nD0Below = 0;
    nD0barBelow = 0;
    for (int bin=51; bin<101; bin++)
    {
        nD0Below = nD0Below+quad2Phi1HistAB->GetBinContent(bin);
        nD0barBelow = nD0barBelow+quad2Phi1DSHistAB->GetBinContent(bin);
    }

    D0AB_asymmetry = (nD0Above-nD0Below)/(nD0Above+nD0Below);
    D0AB_asymmetryErr = sqrt( (1.0 - D0AB_asymmetry*D0AB_asymmetry)/(nD0Above+nD0Below) );

    D0barAB_asymmetry = (nD0barAbove-nD0barBelow)/(nD0barAbove+nD0barBelow);
    D0barAB_asymmetryErr = sqrt( (1.0 - D0barAB_asymmetry*D0barAB_asymmetry)/(nD0barAbove+nD0barBelow) );

    quad2Phi1D0ABAsymmetryStr.Form("%6.2f\n", D0AB_asymmetry);
    quad2Phi1D0ABAsymmetryErrStr.Form("%6.2f\n", D0AB_asymmetryErr); 

    quad2Phi1D0barABAsymmetryStr.Form("%6.2f\n", D0barAB_asymmetry); 
    quad2Phi1D0barABAsymmetryErrStr.Form("%6.2f\n", D0barAB_asymmetryErr); 

    //std::cout<< D0AB_asymmetry <<std::endl; 
    //std::cout<<D0barAB_asymmetry <<std::endl; 

    weight1 = 1./(D0AB_asymmetryErr*D0AB_asymmetryErr);
    weight2 = 1./(D0barAB_asymmetryErr*D0barAB_asymmetryErr);

  CPV = (D0AB_asymmetry*weight1+D0barAB_asymmetry*weight2)/(weight1+weight2);
  CPVErr = sqrt(1./(weight1+weight2));

  //std::cout<<"CPA " << CPV << std::endl; 

  quad2Phi1AsymmetryABStr.Form("%6.4f\n",CPV); 
  quad2Phi1AsymmetryErrABStr.Form("%6.4f\n", CPVErr); 

  for(int bin = 1; bin < 26; bin++)
    {
        nD02PhiAbove = nD02PhiAbove+quad2Phi1HistAB->GetBinContent(bin); 
        nD0bar2PhiAbove = nD0bar2PhiAbove+quad2Phi1DSHistAB->GetBinContent(bin); 
    }

    for(int bin = 26; bin < 51; bin++)
    {
        nD0bar2PhiBelow = nD0bar2PhiBelow + quad2Phi1DSHistAB->GetBinContent(bin); 
        nD02PhiBelow = nD02PhiBelow + quad2Phi1HistAB->GetBinContent(bin); 
    }

    for(int bin = 51; bin < 76; bin++)
    {
        nD02PhiAbove = nD02PhiAbove+quad2Phi1HistAB->GetBinContent(bin); 
        nD0bar2PhiAbove = nD0bar2PhiAbove+quad2Phi1DSHistAB->GetBinContent(bin); 
    }

    for(int bin = 76; bin <101; bin++)
    {
        nD0bar2PhiBelow = nD0bar2PhiBelow + quad2Phi1DSHistAB->GetBinContent(bin); 
        nD02PhiBelow = nD02PhiBelow + quad2Phi1HistAB->GetBinContent(bin); 
    }

    D0Sin2Phi_asymmetry = (nD02PhiAbove - nD02PhiBelow)/(nD02PhiAbove+nD02PhiBelow); 
    D0Sin2Phi_asymmetryErr = sqrt((1-D0Sin2Phi_asymmetry*D0Sin2Phi_asymmetry)/( nD02PhiAbove + nD02PhiBelow )); 

    D0barSin2Phi_asymmetry = (nD0bar2PhiAbove - nD0bar2PhiBelow)/(nD0bar2PhiAbove+nD0bar2PhiBelow); 
    D0barSin2Phi_asymmetryErr = sqrt((1-D0barSin2Phi_asymmetry*D0barSin2Phi_asymmetry)/( nD0bar2PhiAbove + nD0bar2PhiBelow ));

    weight1 = 1./(D0Sin2Phi_asymmetry*D0Sin2Phi_asymmetry); 
    weight2 = 1./(D0barSin2Phi_asymmetry*D0barSin2Phi_asymmetry); 

    Sin2PhiCP_asymmetry = (D0Sin2Phi_asymmetry*weight1+D0barSin2Phi_asymmetry*weight2)/(weight1+weight2);
    Sin2PhiCP_asymmetryErr = sqrt(1./(weight1+weight2));

    DiffSin2Phi = (D0Sin2Phi_asymmetry*weight1-D0barSin2Phi_asymmetry*weight2)/(weight1+weight2);

    TString D0Sin2Phi_asymmetryStr2; 
    TString D0Sin2Phi_asymmetryErrStr2; 
    TString D0barSin2Phi_asymmetryStr2; 
    TString D0barSin2Phi_asymmetryErrStr2; 
    TString Sin2PhiCP_asymmetryStr2; 
    TString Sin2PhiCP_asymmetryErrStr2; 
    TString DiffSin2PhiStr2;

    D0Sin2Phi_asymmetryStr2.Form("%6.2f\n", D0Sin2Phi_asymmetry); 
    D0Sin2Phi_asymmetryErrStr2.Form("%6.2f\n", D0Sin2Phi_asymmetryErr); 
    D0barSin2Phi_asymmetryStr2.Form("%6.2f\n", D0barSin2Phi_asymmetry);
    D0barSin2Phi_asymmetryErrStr2.Form("%6.2f\n", D0barSin2Phi_asymmetryErr); 
    Sin2PhiCP_asymmetryStr2.Form("%6.4f\n", Sin2PhiCP_asymmetry); 
    Sin2PhiCP_asymmetryErrStr2.Form("%6.4f\n", Sin2PhiCP_asymmetryErr); 
    DiffSin2PhiStr2.Form("%6.2f\n", DiffSin2Phi); 

  quad2Phi1CanvasAB->cd(); 
  quad2Phi1HistAB->GetXaxis()->CenterTitle(); 
  quad2Phi1HistAB->GetYaxis()->CenterTitle(); 
  quad2Phi1HistAB->GetXaxis()->SetTitleOffset(1.25); 
  quad2Phi1HistAB->GetYaxis()->SetTitleOffset(1.35);
  quad2Phi1HistAB->Draw(); 
  quad2Phi1CanvasAB->SaveAs("quad2Phi1HistAB.png"); 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    nD0 = quad3Phi1HistAB->GetEntries();
    nD0bar = quad3Phi1DSHistAB->GetEntries();
      
    scale = nD0/nD0bar;
    nD0Above = 0;
    nD0barAbove = 0;

    for (int bin=0; bin<51; bin++) 
    {
      nD0Above = nD0Above+quad3Phi1HistAB->GetBinContent(bin);
      nD0barAbove = nD0barAbove+quad3Phi1DSHistAB->GetBinContent(bin);
  }
   
  nD0Below = 0;
    nD0barBelow = 0;
    for (int bin=51; bin<101; bin++)
    {
        nD0Below = nD0Below+quad3Phi1HistAB->GetBinContent(bin);
      nD0barBelow = nD0barBelow+quad3Phi1DSHistAB->GetBinContent(bin);
    }

    D0AB_asymmetry = (nD0Above-nD0Below)/(nD0Above+nD0Below);
    D0AB_asymmetryErr = sqrt( (1.0 - D0AB_asymmetry*D0AB_asymmetry)/(nD0Above+nD0Below) );

    D0barAB_asymmetry = (nD0barAbove-nD0barBelow)/(nD0barAbove+nD0barBelow);
    D0barAB_asymmetryErr = sqrt( (1.0 - D0barAB_asymmetry*D0barAB_asymmetry)/(nD0barAbove+nD0barBelow) );

    quad3Phi1D0ABAsymmetryStr.Form("%6.2f\n", D0AB_asymmetry);
    quad3Phi1D0ABAsymmetryErrStr.Form("%6.2f\n", D0AB_asymmetryErr); 

    quad3Phi1D0barABAsymmetryStr.Form("%6.2f\n", D0barAB_asymmetry); 
    quad3Phi1D0barABAsymmetryErrStr.Form("%6.2f\n", D0barAB_asymmetryErr); 

    //std::cout<< D0AB_asymmetry <<std::endl; 
    //std::cout<<D0barAB_asymmetry <<std::endl; 

  weight1 = 1./(D0AB_asymmetryErr*D0AB_asymmetryErr);
    weight2 = 1./(D0barAB_asymmetryErr*D0barAB_asymmetryErr);

  CPV = (D0AB_asymmetry*weight1+D0barAB_asymmetry*weight2)/(weight1+weight2);
  CPVErr = sqrt(1./(weight1+weight2));

  quad3Phi1AsymmetryABStr.Form("%6.4f\n",CPV); 
  quad3Phi1AsymmetryErrABStr.Form("%6.4f\n", CPVErr); 

  //std::cout<<"CPA " << CPV << std::endl; 

  for(int bin = 1; bin < 26; bin++)
    {
        nD02PhiAbove = nD02PhiAbove+quad3Phi1HistAB->GetBinContent(bin); 
        nD0bar2PhiAbove = nD0bar2PhiAbove+quad3Phi1DSHistAB->GetBinContent(bin); 
    }

    for(int bin = 26; bin < 51; bin++)
    {
        nD0bar2PhiBelow = nD0bar2PhiBelow + quad3Phi1DSHistAB->GetBinContent(bin); 
        nD02PhiBelow = nD02PhiBelow + quad3Phi1HistAB->GetBinContent(bin); 
    }

    for(int bin = 51; bin < 76; bin++)
    {
        nD02PhiAbove = nD02PhiAbove+quad3Phi1HistAB->GetBinContent(bin); 
        nD0bar2PhiAbove = nD0bar2PhiAbove+quad3Phi1DSHistAB->GetBinContent(bin); 
    }

    for(int bin = 76; bin <101; bin++)
    {
        nD0bar2PhiBelow = nD0bar2PhiBelow + quad3Phi1DSHistAB->GetBinContent(bin); 
        nD02PhiBelow = nD02PhiBelow + quad3Phi1HistAB->GetBinContent(bin); 
    }

    D0Sin2Phi_asymmetry = (nD02PhiAbove - nD02PhiBelow)/(nD02PhiAbove+nD02PhiBelow); 
    D0Sin2Phi_asymmetryErr = sqrt((1-D0Sin2Phi_asymmetry*D0Sin2Phi_asymmetry)/( nD02PhiAbove + nD02PhiBelow )); 

    D0barSin2Phi_asymmetry = (nD0bar2PhiAbove - nD0bar2PhiBelow)/(nD0bar2PhiAbove+nD0bar2PhiBelow); 
    D0barSin2Phi_asymmetryErr = sqrt((1-D0barSin2Phi_asymmetry*D0barSin2Phi_asymmetry)/( nD0bar2PhiAbove + nD0bar2PhiBelow ));

    weight1 = 1./(D0Sin2Phi_asymmetry*D0Sin2Phi_asymmetry); 
    weight2 = 1./(D0barSin2Phi_asymmetry*D0barSin2Phi_asymmetry); 

    Sin2PhiCP_asymmetry = (D0Sin2Phi_asymmetry*weight1+D0barSin2Phi_asymmetry*weight2)/(weight1+weight2);
    Sin2PhiCP_asymmetryErr = sqrt(1./(weight1+weight2));

    DiffSin2Phi = (D0Sin2Phi_asymmetry*weight1-D0barSin2Phi_asymmetry*weight2)/(weight1+weight2);

    TString D0Sin2Phi_asymmetryStr3; 
    TString D0Sin2Phi_asymmetryErrStr3; 
    TString D0barSin2Phi_asymmetryStr3; 
    TString D0barSin2Phi_asymmetryErrStr3; 
    TString Sin2PhiCP_asymmetryStr3; 
    TString Sin2PhiCP_asymmetryErrStr3; 
    TString DiffSin2PhiStr3;

    D0Sin2Phi_asymmetryStr3.Form("%6.2f\n", D0Sin2Phi_asymmetry); 
    D0Sin2Phi_asymmetryErrStr3.Form("%6.2f\n", D0Sin2Phi_asymmetryErr); 
    D0barSin2Phi_asymmetryStr3.Form("%6.2f\n", D0barSin2Phi_asymmetry);
    D0barSin2Phi_asymmetryErrStr3.Form("%6.2f\n", D0barSin2Phi_asymmetryErr); 
    Sin2PhiCP_asymmetryStr3.Form("%6.4f\n", Sin2PhiCP_asymmetry); 
    Sin2PhiCP_asymmetryErrStr3.Form("%6.4f\n", Sin2PhiCP_asymmetryErr); 
    DiffSin2PhiStr3.Form("%6.2f\n", DiffSin2Phi); 

  quad3Phi1CanvasAB->cd(); 
  quad3Phi1HistAB->GetXaxis()->CenterTitle(); 
  quad3Phi1HistAB->GetYaxis()->CenterTitle(); 
  quad3Phi1HistAB->GetXaxis()->SetTitleOffset(1.25); 
  quad3Phi1HistAB->GetYaxis()->SetTitleOffset(1.35);
  quad3Phi1HistAB->Draw(); 
  quad3Phi1CanvasAB->SaveAs("quad3Phi1HistAB.png"); 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  nD0 = quad4Phi1HistAB->GetEntries();
    nD0bar = quad4Phi1DSHistAB->GetEntries();
  scale = nD0/nD0bar;
    nD0Above = 0;
    nD0barAbove = 0;

    for (int bin=0; bin<51; bin++) 
    {
      nD0Above = nD0Above+quad4Phi1HistAB->GetBinContent(bin);
      nD0barAbove = nD0barAbove+quad4Phi1DSHistAB->GetBinContent(bin);
  }
   
  nD0Below = 0;
    nD0barBelow = 0;
    for (int bin=51; bin<101; bin++)
    {
      nD0Below = nD0Below+quad4Phi1HistAB->GetBinContent(bin);
      nD0barBelow = nD0barBelow+quad4Phi1DSHistAB->GetBinContent(bin);
    }

    D0AB_asymmetry = (nD0Above-nD0Below)/(nD0Above+nD0Below);
    D0AB_asymmetryErr = sqrt( (1.0 - D0AB_asymmetry*D0AB_asymmetry)/(nD0Above+nD0Below) );

    D0barAB_asymmetry = (nD0barAbove-nD0barBelow)/(nD0barAbove+nD0barBelow);
    D0barAB_asymmetryErr = sqrt( (1.0 - D0barAB_asymmetry*D0barAB_asymmetry)/(nD0barAbove+nD0barBelow) );

    quad4Phi1D0ABAsymmetryStr.Form("%6.2f\n", D0AB_asymmetry);
    quad4Phi1D0ABAsymmetryErrStr.Form("%6.2f\n", D0AB_asymmetryErr); 

    quad4Phi1D0barABAsymmetryStr.Form("%6.2f\n", D0barAB_asymmetry); 
    quad4Phi1D0barABAsymmetryErrStr.Form("%6.2f\n", D0barAB_asymmetryErr); 

    //std::cout<< D0AB_asymmetry <<std::endl; 
    //std::cout<<D0barAB_asymmetry <<std::endl; 

  weight1 = 1./(D0AB_asymmetryErr*D0AB_asymmetryErr);
    weight2 = 1./(D0barAB_asymmetryErr*D0barAB_asymmetryErr);

  CPV = (D0AB_asymmetry*weight1+D0barAB_asymmetry*weight2)/(weight1+weight2);
  CPVErr = sqrt(1./(weight1+weight2));

  quad4Phi1AsymmetryABStr.Form("%6.4f\n",CPV); 
  quad4Phi1AsymmetryErrABStr.Form("%6.4f\n", CPVErr); 

  //std::cout<<"CPA " << CPV << std::endl; 

  for(int bin = 1; bin < 26; bin++)
    {
        nD02PhiAbove = nD02PhiAbove+quad4Phi1HistAB->GetBinContent(bin); 
        nD0bar2PhiAbove = nD0bar2PhiAbove+quad4Phi1DSHistAB->GetBinContent(bin); 
    }

    for(int bin = 26; bin < 51; bin++)
    {
        nD0bar2PhiBelow = nD0bar2PhiBelow + quad4Phi1DSHistAB->GetBinContent(bin); 
        nD02PhiBelow = nD02PhiBelow + quad4Phi1HistAB->GetBinContent(bin); 
    }

    for(int bin = 51; bin < 76; bin++)
    {
        nD02PhiAbove = nD02PhiAbove+quad4Phi1HistAB->GetBinContent(bin); 
        nD0bar2PhiAbove = nD0bar2PhiAbove+quad4Phi1DSHistAB->GetBinContent(bin); 
    }

    for(int bin = 76; bin <101; bin++)
    {
        nD0bar2PhiBelow = nD0bar2PhiBelow + quad4Phi1DSHistAB->GetBinContent(bin); 
        nD02PhiBelow = nD02PhiBelow + quad4Phi1HistAB->GetBinContent(bin); 
    }

    D0Sin2Phi_asymmetry = (nD02PhiAbove - nD02PhiBelow)/(nD02PhiAbove+nD02PhiBelow); 
    D0Sin2Phi_asymmetryErr = sqrt((1-D0Sin2Phi_asymmetry*D0Sin2Phi_asymmetry)/( nD02PhiAbove + nD02PhiBelow )); 

    D0barSin2Phi_asymmetry = (nD0bar2PhiAbove - nD0bar2PhiBelow)/(nD0bar2PhiAbove+nD0bar2PhiBelow); 
    D0barSin2Phi_asymmetryErr = sqrt((1-D0barSin2Phi_asymmetry*D0barSin2Phi_asymmetry)/( nD0bar2PhiAbove + nD0bar2PhiBelow ));

    weight1 = 1./(D0Sin2Phi_asymmetry*D0Sin2Phi_asymmetry); 
    weight2 = 1./(D0barSin2Phi_asymmetry*D0barSin2Phi_asymmetry); 

    Sin2PhiCP_asymmetry = (D0Sin2Phi_asymmetry*weight1+D0barSin2Phi_asymmetry*weight2)/(weight1+weight2);
    Sin2PhiCP_asymmetryErr = sqrt(1./(weight1+weight2));

    DiffSin2Phi = (D0Sin2Phi_asymmetry*weight1-D0barSin2Phi_asymmetry*weight2)/(weight1+weight2);

    TString D0Sin2Phi_asymmetryStr4; 
    TString D0Sin2Phi_asymmetryErrStr4; 
    TString D0barSin2Phi_asymmetryStr4; 
    TString D0barSin2Phi_asymmetryErrStr4; 
    TString Sin2PhiCP_asymmetryStr4; 
    TString Sin2PhiCP_asymmetryErrStr4; 
    TString DiffSin2PhiStr4;

    D0Sin2Phi_asymmetryStr4.Form("%6.2f\n", D0Sin2Phi_asymmetry); 
    D0Sin2Phi_asymmetryErrStr4.Form("%6.2f\n", D0Sin2Phi_asymmetryErr); 
    D0barSin2Phi_asymmetryStr4.Form("%6.2f\n", D0barSin2Phi_asymmetry);
    D0barSin2Phi_asymmetryErrStr4.Form("%6.2f\n", D0barSin2Phi_asymmetryErr); 
    Sin2PhiCP_asymmetryStr4.Form("%6.4f\n", Sin2PhiCP_asymmetry); 
    Sin2PhiCP_asymmetryErrStr4.Form("%6.4f\n", Sin2PhiCP_asymmetryErr); 
    DiffSin2PhiStr4.Form("%6.2f\n", DiffSin2Phi); 

  quad4Phi1CanvasAB->cd(); 
  quad4Phi1HistAB->GetXaxis()->CenterTitle(); 
  quad4Phi1HistAB->GetYaxis()->CenterTitle(); 
  quad4Phi1HistAB->GetXaxis()->SetTitleOffset(1.25); 
  quad4Phi1HistAB->GetYaxis()->SetTitleOffset(1.35);
  quad4Phi1HistAB->Draw(); 
  quad4Phi1CanvasAB->SaveAs("quad4Phi1HistAB.png"); 

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  quad1Phi1DSCanvasAB->cd(); 
  quad1Phi1DSHistAB->GetXaxis()->CenterTitle(); 
  quad1Phi1DSHistAB->GetYaxis()->CenterTitle(); 
  quad1Phi1DSHistAB->GetXaxis()->SetTitleOffset(1.25); 
  quad1Phi1DSHistAB->GetYaxis()->SetTitleOffset(1.35);
  quad1Phi1DSHistAB->Draw(); 
  quad1Phi1DSCanvasAB->SaveAs("quad1Phi1DSHistAB.png"); 

  

  quad2Phi1DSCanvasAB->cd(); 
  quad2Phi1DSHistAB->GetXaxis()->CenterTitle(); 
  quad2Phi1DSHistAB->GetYaxis()->CenterTitle(); 
  quad2Phi1DSHistAB->GetXaxis()->SetTitleOffset(1.25); 
  quad2Phi1DSHistAB->GetYaxis()->SetTitleOffset(1.35);
  quad2Phi1DSHistAB->Draw(); 
  quad2Phi1DSCanvasAB->SaveAs("quad2Phi1DSHistAB.png"); 

  

  quad3Phi1DSCanvasAB->cd();
  quad3Phi1DSHistAB->GetXaxis()->CenterTitle(); 
  quad3Phi1DSHistAB->GetYaxis()->CenterTitle(); 
  quad3Phi1DSHistAB->GetXaxis()->SetTitleOffset(1.25); 
  quad3Phi1DSHistAB->GetYaxis()->SetTitleOffset(1.35); 
  quad3Phi1DSHistAB->Draw(); 
  quad3Phi1DSCanvasAB->SaveAs("quad3Phi1DSHistAB.png"); 

  

  quad4Phi1DSCanvasAB->cd(); 
  quad4Phi1DSHistAB->GetXaxis()->CenterTitle(); 
  quad4Phi1DSHistAB->GetYaxis()->CenterTitle(); 
  quad4Phi1DSHistAB->GetXaxis()->SetTitleOffset(1.25); 
  quad4Phi1DSHistAB->GetYaxis()->SetTitleOffset(1.35);
  quad4Phi1DSHistAB->Draw(); 
  quad4Phi1DSCanvasAB->SaveAs("quad4Phi1DSHistAB.png"); 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double maxBin; 
  double binContent; 
  TLegend *leg = new TLegend(0.15,0.85,0.3,0.95);
  TPaveText *pt1 = new TPaveText(0.7,0.75,0.98,0.98,"NDC"); // NDC sets coords
  pt1->SetTextSize(0.025); 

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  quad1Phi1CompareABCanvas->cd(); 
  quad1Phi1HistAB->SetLineColor(kBlue); 
  quad1Phi1DSHistAB->SetLineColor(kRed); 
  quad1Phi1HistAB->GetYaxis()->SetLimits(0,1650);
  quad1Phi1HistAB->SetMinimum(0);
  maxBin = quad1Phi1HistAB->GetMaximumBin(); 
  binContent = 1.3*quad1Phi1HistAB->GetBinContent(maxBin); 
  quad1Phi1HistAB->SetMaximum(binContent); 
  TLine *line1A = new TLine(3.14,0,3.14, binContent); 
  leg->SetHeader("Legend"); 
  leg->AddEntry(quad1Phi1HistAB, "Dzero"); 
  leg->AddEntry(quad1Phi1DSHistAB, "Dzero-bar"); 
  pt1->AddText("PV D0:" + quad1Phi1D0ABAsymmetryStr + " #pm " + quad1Phi1D0ABAsymmetryErrStr);
  pt1->AddText("PV D0bar: " + quad1Phi1D0barABAsymmetryStr + " #pm " + quad1Phi1D0barABAsymmetryErrStr); 
  pt1->AddText("CPV: " + quad1Phi1AsymmetryABStr + " #pm" + quad1Phi1AsymmetryErrABStr);
  pt1->AddText("PV D0 Sin(2#phi): " + D0Sin2Phi_asymmetryStr1 + " #pm " + D0Sin2Phi_asymmetryErrStr1); 
  pt1->AddText("PV D0bar Sin(2#phi): " + D0barSin2Phi_asymmetryStr1 + " #pm " + D0barSin2Phi_asymmetryErrStr1); 
  pt1->AddText("CPV Sin(2#phi): " + Sin2PhiCP_asymmetryStr1 + " #pm " + Sin2PhiCP_asymmetryErrStr1); 
  pt1->AddText("D0-D0bar : " + DiffSin2PhiStr1); 
  quad1Phi1HistAB->SetTitle("Quad 1: " + DecayOne);
  quad1Phi1HistAB->Draw(); 
  quad1Phi1DSHistAB->Draw("same");
  line1A->Draw("same");
  leg->Draw("same"); 
  pt1->Draw("same");  
  quad1Phi1CompareABCanvas->SaveAs("quad1Phi1CompAB.png"); 
  pt1->Clear(); 
  leg->Clear(); 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////   

  quad2Phi1CompareABCanvas->cd(); 
  quad2Phi1HistAB->SetLineColor(kBlue); 
  quad2Phi1DSHistAB->SetLineColor(kRed); 
  quad2Phi1HistAB->GetYaxis()->SetLimits(0,1650);
  quad2Phi1HistAB->SetMinimum(0);
  maxBin = quad2Phi1HistAB->GetMaximumBin(); 
  binContent = 1.3*quad2Phi1HistAB->GetBinContent(maxBin); 
  quad2Phi1HistAB->SetMaximum(binContent); 
  TLine *line2A = new TLine(3.14,0,3.14, binContent); 
  leg->SetHeader("Legend"); 
  leg->AddEntry(quad2Phi1HistAB, "Dzero"); 
  leg->AddEntry(quad2Phi1DSHistAB, "Dzero-bar"); 
  pt1->AddText("PV D0:" + quad2Phi1D0ABAsymmetryStr + " #pm " + quad2Phi1D0ABAsymmetryErrStr);
  pt1->AddText("PV D0bar: " + quad2Phi1D0barABAsymmetryStr + " #pm " + quad2Phi1D0barABAsymmetryErrStr); 
  pt1->AddText("CPV: " + quad2Phi1AsymmetryABStr + " #pm" + quad2Phi1AsymmetryErrABStr);
  pt1->AddText("PV D0 Sin(2#phi): " + D0Sin2Phi_asymmetryStr2 + " #pm " + D0Sin2Phi_asymmetryErrStr2); 
  pt1->AddText("PV D0bar Sin(2#phi): " + D0barSin2Phi_asymmetryStr2 + " #pm " + D0barSin2Phi_asymmetryErrStr2); 
  pt1->AddText("CPV Sin(2#phi): " + Sin2PhiCP_asymmetryStr2 + " #pm " + Sin2PhiCP_asymmetryErrStr2); 
  pt1->AddText("D0-D0bar : " + DiffSin2PhiStr2); 
  quad2Phi1HistAB->SetTitle("Quad 2: " + DecayOne);
  quad2Phi1HistAB->Draw(); 
  quad2Phi1DSHistAB->Draw("same");
  line2A->Draw("same");
  leg->Draw("same"); 
  pt1->Draw("same");  
  quad2Phi1CompareABCanvas->SaveAs("quad2Phi1CompAB.png"); 
  pt1->Clear(); 
  leg->Clear(); 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  

  quad3Phi1CompareABCanvas->cd(); 
  quad3Phi1HistAB->SetLineColor(kBlue); 
  quad3Phi1DSHistAB->SetLineColor(kRed); 
  quad3Phi1HistAB->GetYaxis()->SetLimits(0,1650);
  quad3Phi1HistAB->SetMinimum(0);
  maxBin = quad3Phi1HistAB->GetMaximumBin(); 
  binContent = 1.3*quad3Phi1HistAB->GetBinContent(maxBin); 
  quad3Phi1HistAB->SetMaximum(binContent); 
  leg->SetHeader("Legend"); 
  leg->AddEntry(quad3Phi1HistAB, "Dzero"); 
  leg->AddEntry(quad3Phi1DSHistAB, "Dzero-bar"); 
  pt1->AddText("PV D0:" + quad3Phi1D0ABAsymmetryStr + " #pm " + quad3Phi1D0ABAsymmetryErrStr);
  pt1->AddText("PV D0bar: " + quad3Phi1D0barABAsymmetryStr + " #pm " + quad3Phi1D0barABAsymmetryErrStr); 
  pt1->AddText("CPV: " + quad3Phi1AsymmetryABStr + " #pm" + quad3Phi1AsymmetryErrABStr);
  pt1->AddText("PV D0 Sin(2#phi): " + D0Sin2Phi_asymmetryStr3 + " #pm " + D0Sin2Phi_asymmetryErrStr3); 
  pt1->AddText("PV D0bar Sin(2#phi): " + D0barSin2Phi_asymmetryStr3 + " #pm " + D0barSin2Phi_asymmetryErrStr3); 
  pt1->AddText("CPV Sin(2#phi): " + Sin2PhiCP_asymmetryStr3 + " #pm " + Sin2PhiCP_asymmetryErrStr3); 
  pt1->AddText("D0-D0bar : " + DiffSin2PhiStr3); 
  TLine *line3A = new TLine(3.14,0,3.14, binContent); 
  quad3Phi1HistAB->SetTitle("Quad 3: " + DecayOne);
  quad3Phi1HistAB->Draw(); 
  quad3Phi1DSHistAB->Draw("same");
  line3A->Draw("same");
  leg->Draw("same"); 
  pt1->Draw("same");  
  quad3Phi1CompareABCanvas->SaveAs("quad3Phi1CompAB.png"); 
  leg->Clear("same"); 
  pt1->Clear("same");  
  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  quad4Phi1CompareABCanvas->cd(); 
  quad4Phi1HistAB->SetLineColor(kBlue); 
  quad4Phi1DSHistAB->SetLineColor(kRed); 
  quad4Phi1HistAB->GetYaxis()->SetLimits(0,1650);
  quad4Phi1HistAB->SetMinimum(0);
  maxBin = quad4Phi1HistAB->GetMaximumBin(); 
  binContent = 1.3*quad4Phi1HistAB->GetBinContent(maxBin); 
  quad4Phi1HistAB->SetMaximum(binContent); 
  TLine *line4A = new TLine(3.14,0,3.14, binContent); 
  leg->SetHeader("Legend"); 
  leg->AddEntry(quad4Phi1HistAB, "Dzero"); 
  leg->AddEntry(quad4Phi1DSHistAB, "Dzero-bar"); 
  pt1->AddText("PV D0:" + quad4Phi1D0ABAsymmetryStr + " #pm " + quad4Phi1D0ABAsymmetryErrStr);
  pt1->AddText("PV D0bar: " + quad4Phi1D0barABAsymmetryStr + " #pm " + quad4Phi1D0barABAsymmetryErrStr); 
  pt1->AddText("CPV: " + quad4Phi1AsymmetryABStr + " #pm" + quad4Phi1AsymmetryErrABStr);
  pt1->AddText("PV D0 Sin(2#phi): " + D0Sin2Phi_asymmetryStr4 + " #pm " + D0Sin2Phi_asymmetryErrStr4); 
  pt1->AddText("PV D0bar Sin(2#phi): " + D0barSin2Phi_asymmetryStr4 + " #pm " + D0barSin2Phi_asymmetryErrStr4); 
  pt1->AddText("CPV Sin(2#phi): " + Sin2PhiCP_asymmetryStr4 + " #pm " + Sin2PhiCP_asymmetryErrStr4); 
  pt1->AddText("D0-D0bar : " + DiffSin2PhiStr4); 
  quad4Phi1HistAB->SetTitle("Quad 4: " + DecayOne);
  quad4Phi1HistAB->Draw(); 
  quad4Phi1DSHistAB->Draw("same");
  line4A->Draw("same");
  leg->Draw("same"); 
  pt1->Draw("same");  
  quad4Phi1CompareABCanvas->SaveAs("quad4Phi1CompAB.png"); 
  leg->Clear("same"); 
  pt1->Clear("same");  
}


  nD0 = quad1Phi2HistAB->GetEntries();
  nD0bar = quad1Phi2DSHistAB->GetEntries();
  Double_t scale = nD0/nD0bar;
  Double_t nD0Above = 0;
  Double_t nD0barAbove = 0;

  for (int bin=0; bin<51; bin++) 
  {
    nD0Above = nD0Above+quad1Phi2HistAB->GetBinContent(bin);
    nD0barAbove = nD0barAbove+quad1Phi2DSHistAB->GetBinContent(bin);
  }
   
  Double_t nD0Below = 0;
  Double_t nD0barBelow = 0;
  for (int bin=51; bin<101; bin++)
  {
    nD0Below = nD0Below+quad1Phi2HistAB->GetBinContent(bin);
    nD0barBelow = nD0barBelow+quad1Phi2DSHistAB->GetBinContent(bin);
  }

  Double_t D0AB_asymmetry = (nD0Above-nD0Below)/(nD0Above+nD0Below);
  Double_t D0AB_asymmetryErr = sqrt( (1.0 - D0AB_asymmetry*D0AB_asymmetry)/(nD0Above+nD0Below) );

  Double_t D0barAB_asymmetry = (nD0barAbove-nD0barBelow)/(nD0barAbove+nD0barBelow);
  Double_t D0barAB_asymmetryErr = sqrt( (1.0 - D0barAB_asymmetry*D0barAB_asymmetry)/(nD0barAbove+nD0barBelow) );

  quad1Phi2D0ABAsymmetryStr.Form("%6.2f\n", D0AB_asymmetry);
  quad1Phi2D0ABAsymmetryErrStr.Form("%6.2f\n", D0AB_asymmetryErr); 

  quad1Phi2D0barABAsymmetryStr.Form("%6.2f\n", D0barAB_asymmetry); 
  quad1Phi2D0barABAsymmetryErrStr.Form("%6.2f\n", D0barAB_asymmetryErr); 

  Double_t weight1 = 1./(D0AB_asymmetryErr*D0AB_asymmetryErr);
  Double_t weight2 = 1./(D0barAB_asymmetryErr*D0barAB_asymmetryErr);

  Double_t CPV = (D0AB_asymmetry*weight1+D0barAB_asymmetry*weight2)/(weight1+weight2);
  Double_t CPVErr = sqrt(1./(weight1+weight2));

  quad1Phi2AsymmetryABStr.Form("%6.4f\n",CPV); 
  quad1Phi2AsymmetryErrABStr.Form("%6.4f\n", CPVErr); 

  Double_t nD02PhiAbove = 0; 
  Double_t nD0bar2PhiAbove = 0;
  Double_t nD02PhiBelow = 0; 
  Double_t nD0bar2PhiBelow = 0;  

  for(int bin = 1; bin < 26; bin++)
  {
    nD02PhiAbove = nD02PhiAbove+quad1Phi2HistAB->GetBinContent(bin); 
    nD0bar2PhiAbove = nD0bar2PhiAbove+quad1Phi2DSHistAB->GetBinContent(bin); 
  }

  for(int bin = 26; bin < 51; bin++)
  {
    nD0bar2PhiBelow = nD0bar2PhiBelow + quad1Phi2DSHistAB->GetBinContent(bin); 
    nD02PhiBelow = nD02PhiBelow + quad1Phi2HistAB->GetBinContent(bin); 
  }

  for(int bin = 51; bin < 76; bin++)
  {
    nD02PhiAbove = nD02PhiAbove+quad1Phi2HistAB->GetBinContent(bin); 
    nD0bar2PhiAbove = nD0bar2PhiAbove+quad1Phi2DSHistAB->GetBinContent(bin); 
  }

  for(int bin = 76; bin <101; bin++)
  {
    nD0bar2PhiBelow = nD0bar2PhiBelow + quad1Phi2DSHistAB->GetBinContent(bin); 
    nD02PhiBelow = nD02PhiBelow + quad1Phi2HistAB->GetBinContent(bin); 
  }

    

  D0Sin2Phi_asymmetry = (nD02PhiAbove - nD02PhiBelow)/(nD02PhiAbove+nD02PhiBelow); 
  D0Sin2Phi_asymmetryErr = sqrt((1-D0Sin2Phi_asymmetry*D0Sin2Phi_asymmetry)/( nD02PhiAbove + nD02PhiBelow )); 

  D0barSin2Phi_asymmetry = (nD0bar2PhiAbove - nD0bar2PhiBelow)/(nD0bar2PhiAbove+nD0bar2PhiBelow); 
  D0barSin2Phi_asymmetryErr = sqrt((1-D0barSin2Phi_asymmetry*D0barSin2Phi_asymmetry)/( nD0bar2PhiAbove + nD0bar2PhiBelow ));

  weight1 = 1./(D0Sin2Phi_asymmetry*D0Sin2Phi_asymmetry); 
  weight2 = 1./(D0barSin2Phi_asymmetry*D0barSin2Phi_asymmetry); 

  Sin2PhiCP_asymmetry = (D0Sin2Phi_asymmetry*weight1+D0barSin2Phi_asymmetry*weight2)/(weight1+weight2);
  Sin2PhiCP_asymmetryErr = sqrt(1./(weight1+weight2));
  DiffSin2Phi = (D0Sin2Phi_asymmetry*weight1-D0barSin2Phi_asymmetry*weight2)/(weight1+weight2);

  TString D0Sin2Phi_asymmetryStr1; 
  TString D0Sin2Phi_asymmetryErrStr1; 
  TString D0barSin2Phi_asymmetryStr1; 
  TString D0barSin2Phi_asymmetryErrStr1; 
  TString Sin2PhiCP_asymmetryStr1; 
  TString Sin2PhiCP_asymmetryErrStr1; 
  TString DiffSin2PhiStr1; 

  D0Sin2Phi_asymmetryStr1.Form("%6.2f\n", D0Sin2Phi_asymmetry); 
  D0Sin2Phi_asymmetryErrStr1.Form("%6.2f\n", D0Sin2Phi_asymmetryErr); 
  D0barSin2Phi_asymmetryStr1.Form("%6.2f\n", D0barSin2Phi_asymmetry);
  D0barSin2Phi_asymmetryErrStr1.Form("%6.2f\n", D0barSin2Phi_asymmetryErr); 
  Sin2PhiCP_asymmetryStr1.Form("%6.4f\n", Sin2PhiCP_asymmetry); 
  Sin2PhiCP_asymmetryErrStr1.Form("%6.4f\n", Sin2PhiCP_asymmetryErr); 
  DiffSin2PhiStr1.Form("%6.2f\n", DiffSin2Phi); 

      

  quad1Phi2CanvasAB->cd(); 
  quad1Phi2HistAB->GetXaxis()->CenterTitle(); 
  quad1Phi2HistAB->GetYaxis()->CenterTitle(); 
  quad1Phi2HistAB->GetXaxis()->SetTitleOffset(1.25); 
  quad1Phi2HistAB->GetYaxis()->SetTitleOffset(1.35);
  quad1Phi2HistAB->Draw(); 
  quad1Phi2CanvasAB->SaveAs("quad1Phi2HistAB.png"); 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  

    nD0 = quad2Phi2HistAB->GetEntries();
    nD0bar = quad2Phi2DSHistAB->GetEntries();
      
    scale = nD0/nD0bar;
    nD0Above = 0;
    nD0barAbove = 0;

    for (int bin=0; bin<51; bin++) 
    {
        nD0Above = nD0Above+quad2Phi2HistAB->GetBinContent(bin);
        nD0barAbove = nD0barAbove+quad2Phi2DSHistAB->GetBinContent(bin);
    }
   
    nD0Below = 0;
    nD0barBelow = 0;
    for (int bin=51; bin<101; bin++)
    {
        nD0Below = nD0Below+quad2Phi2HistAB->GetBinContent(bin);
        nD0barBelow = nD0barBelow+quad2Phi2DSHistAB->GetBinContent(bin);
    }

    D0AB_asymmetry = (nD0Above-nD0Below)/(nD0Above+nD0Below);
    D0AB_asymmetryErr = sqrt( (1.0 - D0AB_asymmetry*D0AB_asymmetry)/(nD0Above+nD0Below) );

    D0barAB_asymmetry = (nD0barAbove-nD0barBelow)/(nD0barAbove+nD0barBelow);
    D0barAB_asymmetryErr = sqrt( (1.0 - D0barAB_asymmetry*D0barAB_asymmetry)/(nD0barAbove+nD0barBelow) );

    quad2Phi2D0ABAsymmetryStr.Form("%6.2f\n", D0AB_asymmetry);
    quad2Phi2D0ABAsymmetryErrStr.Form("%6.2f\n", D0AB_asymmetryErr); 

    quad2Phi2D0barABAsymmetryStr.Form("%6.2f\n", D0barAB_asymmetry); 
    quad2Phi2D0barABAsymmetryErrStr.Form("%6.2f\n", D0barAB_asymmetryErr); 

    //std::cout<< D0AB_asymmetry <<std::endl; 
    //std::cout<<D0barAB_asymmetry <<std::endl; 

    weight1 = 1./(D0AB_asymmetryErr*D0AB_asymmetryErr);
    weight2 = 1./(D0barAB_asymmetryErr*D0barAB_asymmetryErr);

  CPV = (D0AB_asymmetry*weight1+D0barAB_asymmetry*weight2)/(weight1+weight2);
  CPVErr = sqrt(1./(weight1+weight2));

  //std::cout<<"CPA " << CPV << std::endl; 

  quad2Phi2AsymmetryABStr.Form("%6.4f\n",CPV); 
  quad2Phi2AsymmetryErrABStr.Form("%6.4f\n", CPVErr); 

  for(int bin = 1; bin < 26; bin++)
    {
        nD02PhiAbove = nD02PhiAbove+quad2Phi2HistAB->GetBinContent(bin); 
        nD0bar2PhiAbove = nD0bar2PhiAbove+quad2Phi2DSHistAB->GetBinContent(bin); 
    }

    for(int bin = 26; bin < 51; bin++)
    {
        nD0bar2PhiBelow = nD0bar2PhiBelow + quad2Phi2DSHistAB->GetBinContent(bin); 
        nD02PhiBelow = nD02PhiBelow + quad2Phi2HistAB->GetBinContent(bin); 
    }

    for(int bin = 51; bin < 76; bin++)
    {
        nD02PhiAbove = nD02PhiAbove+quad2Phi2HistAB->GetBinContent(bin); 
        nD0bar2PhiAbove = nD0bar2PhiAbove+quad2Phi2DSHistAB->GetBinContent(bin); 
    }

    for(int bin = 76; bin <101; bin++)
    {
        nD0bar2PhiBelow = nD0bar2PhiBelow + quad2Phi2DSHistAB->GetBinContent(bin); 
        nD02PhiBelow = nD02PhiBelow + quad2Phi2HistAB->GetBinContent(bin); 
    }

    D0Sin2Phi_asymmetry = (nD02PhiAbove - nD02PhiBelow)/(nD02PhiAbove+nD02PhiBelow); 
    D0Sin2Phi_asymmetryErr = sqrt((1-D0Sin2Phi_asymmetry*D0Sin2Phi_asymmetry)/( nD02PhiAbove + nD02PhiBelow )); 

    D0barSin2Phi_asymmetry = (nD0bar2PhiAbove - nD0bar2PhiBelow)/(nD0bar2PhiAbove+nD0bar2PhiBelow); 
    D0barSin2Phi_asymmetryErr = sqrt((1-D0barSin2Phi_asymmetry*D0barSin2Phi_asymmetry)/( nD0bar2PhiAbove + nD0bar2PhiBelow ));

    weight1 = 1./(D0Sin2Phi_asymmetry*D0Sin2Phi_asymmetry); 
    weight2 = 1./(D0barSin2Phi_asymmetry*D0barSin2Phi_asymmetry); 

    Sin2PhiCP_asymmetry = (D0Sin2Phi_asymmetry*weight1+D0barSin2Phi_asymmetry*weight2)/(weight1+weight2);
    Sin2PhiCP_asymmetryErr = sqrt(1./(weight1+weight2));

    DiffSin2Phi = (D0Sin2Phi_asymmetry*weight1-D0barSin2Phi_asymmetry*weight2)/(weight1+weight2);

    TString D0Sin2Phi_asymmetryStr2; 
    TString D0Sin2Phi_asymmetryErrStr2; 
    TString D0barSin2Phi_asymmetryStr2; 
    TString D0barSin2Phi_asymmetryErrStr2; 
    TString Sin2PhiCP_asymmetryStr2; 
    TString Sin2PhiCP_asymmetryErrStr2; 
    TString DiffSin2PhiStr2;

    D0Sin2Phi_asymmetryStr2.Form("%6.2f\n", D0Sin2Phi_asymmetry); 
    D0Sin2Phi_asymmetryErrStr2.Form("%6.2f\n", D0Sin2Phi_asymmetryErr); 
    D0barSin2Phi_asymmetryStr2.Form("%6.2f\n", D0barSin2Phi_asymmetry);
    D0barSin2Phi_asymmetryErrStr2.Form("%6.2f\n", D0barSin2Phi_asymmetryErr); 
    Sin2PhiCP_asymmetryStr2.Form("%6.4f\n", Sin2PhiCP_asymmetry); 
    Sin2PhiCP_asymmetryErrStr2.Form("%6.4f\n", Sin2PhiCP_asymmetryErr); 
    DiffSin2PhiStr2.Form("%6.2f\n", DiffSin2Phi); 

  quad2Phi2CanvasAB->cd(); 
  quad2Phi2HistAB->GetXaxis()->CenterTitle(); 
  quad2Phi2HistAB->GetYaxis()->CenterTitle(); 
  quad2Phi2HistAB->GetXaxis()->SetTitleOffset(1.25); 
  quad2Phi2HistAB->GetYaxis()->SetTitleOffset(1.35);
  quad2Phi2HistAB->Draw(); 
  quad2Phi2CanvasAB->SaveAs("quad2Phi2HistAB.png"); 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    nD0 = quad3Phi2HistAB->GetEntries();
    nD0bar = quad3Phi2DSHistAB->GetEntries();
      
    scale = nD0/nD0bar;
    nD0Above = 0;
    nD0barAbove = 0;

    for (int bin=0; bin<51; bin++) 
    {
      nD0Above = nD0Above+quad3Phi2HistAB->GetBinContent(bin);
      nD0barAbove = nD0barAbove+quad3Phi2DSHistAB->GetBinContent(bin);
  }
   
  nD0Below = 0;
    nD0barBelow = 0;
    for (int bin=51; bin<101; bin++)
    {
        nD0Below = nD0Below+quad3Phi2HistAB->GetBinContent(bin);
      nD0barBelow = nD0barBelow+quad3Phi2DSHistAB->GetBinContent(bin);
    }

    D0AB_asymmetry = (nD0Above-nD0Below)/(nD0Above+nD0Below);
    D0AB_asymmetryErr = sqrt( (1.0 - D0AB_asymmetry*D0AB_asymmetry)/(nD0Above+nD0Below) );

    D0barAB_asymmetry = (nD0barAbove-nD0barBelow)/(nD0barAbove+nD0barBelow);
    D0barAB_asymmetryErr = sqrt( (1.0 - D0barAB_asymmetry*D0barAB_asymmetry)/(nD0barAbove+nD0barBelow) );

    quad3Phi2D0ABAsymmetryStr.Form("%6.2f\n", D0AB_asymmetry);
    quad3Phi2D0ABAsymmetryErrStr.Form("%6.2f\n", D0AB_asymmetryErr); 

    quad3Phi2D0barABAsymmetryStr.Form("%6.2f\n", D0barAB_asymmetry); 
    quad3Phi2D0barABAsymmetryErrStr.Form("%6.2f\n", D0barAB_asymmetryErr); 

    //std::cout<< D0AB_asymmetry <<std::endl; 
    //std::cout<<D0barAB_asymmetry <<std::endl; 

  weight1 = 1./(D0AB_asymmetryErr*D0AB_asymmetryErr);
    weight2 = 1./(D0barAB_asymmetryErr*D0barAB_asymmetryErr);

  CPV = (D0AB_asymmetry*weight1+D0barAB_asymmetry*weight2)/(weight1+weight2);
  CPVErr = sqrt(1./(weight1+weight2));

  quad3Phi2AsymmetryABStr.Form("%6.4f\n",CPV); 
  quad3Phi2AsymmetryErrABStr.Form("%6.4f\n", CPVErr); 

  //std::cout<<"CPA " << CPV << std::endl; 

  for(int bin = 1; bin < 26; bin++)
    {
        nD02PhiAbove = nD02PhiAbove+quad3Phi2HistAB->GetBinContent(bin); 
        nD0bar2PhiAbove = nD0bar2PhiAbove+quad3Phi2DSHistAB->GetBinContent(bin); 
    }

    for(int bin = 26; bin < 51; bin++)
    {
        nD0bar2PhiBelow = nD0bar2PhiBelow + quad3Phi2DSHistAB->GetBinContent(bin); 
        nD02PhiBelow = nD02PhiBelow + quad3Phi2HistAB->GetBinContent(bin); 
    }

    for(int bin = 51; bin < 76; bin++)
    {
        nD02PhiAbove = nD02PhiAbove+quad3Phi2HistAB->GetBinContent(bin); 
        nD0bar2PhiAbove = nD0bar2PhiAbove+quad3Phi2DSHistAB->GetBinContent(bin); 
    }

    for(int bin = 76; bin <101; bin++)
    {
        nD0bar2PhiBelow = nD0bar2PhiBelow + quad3Phi2DSHistAB->GetBinContent(bin); 
        nD02PhiBelow = nD02PhiBelow + quad3Phi2HistAB->GetBinContent(bin); 
    }

    D0Sin2Phi_asymmetry = (nD02PhiAbove - nD02PhiBelow)/(nD02PhiAbove+nD02PhiBelow); 
    D0Sin2Phi_asymmetryErr = sqrt((1-D0Sin2Phi_asymmetry*D0Sin2Phi_asymmetry)/( nD02PhiAbove + nD02PhiBelow )); 

    D0barSin2Phi_asymmetry = (nD0bar2PhiAbove - nD0bar2PhiBelow)/(nD0bar2PhiAbove+nD0bar2PhiBelow); 
    D0barSin2Phi_asymmetryErr = sqrt((1-D0barSin2Phi_asymmetry*D0barSin2Phi_asymmetry)/( nD0bar2PhiAbove + nD0bar2PhiBelow ));

    weight1 = 1./(D0Sin2Phi_asymmetry*D0Sin2Phi_asymmetry); 
    weight2 = 1./(D0barSin2Phi_asymmetry*D0barSin2Phi_asymmetry); 

    Sin2PhiCP_asymmetry = (D0Sin2Phi_asymmetry*weight1+D0barSin2Phi_asymmetry*weight2)/(weight1+weight2);
    Sin2PhiCP_asymmetryErr = sqrt(1./(weight1+weight2));

    DiffSin2Phi = (D0Sin2Phi_asymmetry*weight1-D0barSin2Phi_asymmetry*weight2)/(weight1+weight2);

    TString D0Sin2Phi_asymmetryStr3; 
    TString D0Sin2Phi_asymmetryErrStr3; 
    TString D0barSin2Phi_asymmetryStr3; 
    TString D0barSin2Phi_asymmetryErrStr3; 
    TString Sin2PhiCP_asymmetryStr3; 
    TString Sin2PhiCP_asymmetryErrStr3; 
    TString DiffSin2PhiStr3;

    D0Sin2Phi_asymmetryStr3.Form("%6.2f\n", D0Sin2Phi_asymmetry); 
    D0Sin2Phi_asymmetryErrStr3.Form("%6.2f\n", D0Sin2Phi_asymmetryErr); 
    D0barSin2Phi_asymmetryStr3.Form("%6.2f\n", D0barSin2Phi_asymmetry);
    D0barSin2Phi_asymmetryErrStr3.Form("%6.2f\n", D0barSin2Phi_asymmetryErr); 
    Sin2PhiCP_asymmetryStr3.Form("%6.4f\n", Sin2PhiCP_asymmetry); 
    Sin2PhiCP_asymmetryErrStr3.Form("%6.4f\n", Sin2PhiCP_asymmetryErr); 
    DiffSin2PhiStr3.Form("%6.2f\n", DiffSin2Phi); 

  quad3Phi2CanvasAB->cd(); 
  quad3Phi2HistAB->GetXaxis()->CenterTitle(); 
  quad3Phi2HistAB->GetYaxis()->CenterTitle(); 
  quad3Phi2HistAB->GetXaxis()->SetTitleOffset(1.25); 
  quad3Phi2HistAB->GetYaxis()->SetTitleOffset(1.35);
  quad3Phi2HistAB->Draw(); 
  quad3Phi2CanvasAB->SaveAs("quad3Phi2HistAB.png"); 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  nD0 = quad4Phi2HistAB->GetEntries();
    nD0bar = quad4Phi2DSHistAB->GetEntries();
  scale = nD0/nD0bar;
    nD0Above = 0;
    nD0barAbove = 0;

    for (int bin=0; bin<51; bin++) 
    {
      nD0Above = nD0Above+quad4Phi2HistAB->GetBinContent(bin);
      nD0barAbove = nD0barAbove+quad4Phi2DSHistAB->GetBinContent(bin);
  }
   
  nD0Below = 0;
    nD0barBelow = 0;
    for (int bin=51; bin<101; bin++)
    {
      nD0Below = nD0Below+quad4Phi2HistAB->GetBinContent(bin);
      nD0barBelow = nD0barBelow+quad4Phi2DSHistAB->GetBinContent(bin);
    }

    D0AB_asymmetry = (nD0Above-nD0Below)/(nD0Above+nD0Below);
    D0AB_asymmetryErr = sqrt( (1.0 - D0AB_asymmetry*D0AB_asymmetry)/(nD0Above+nD0Below) );

    D0barAB_asymmetry = (nD0barAbove-nD0barBelow)/(nD0barAbove+nD0barBelow);
    D0barAB_asymmetryErr = sqrt( (1.0 - D0barAB_asymmetry*D0barAB_asymmetry)/(nD0barAbove+nD0barBelow) );

    quad4Phi2D0ABAsymmetryStr.Form("%6.2f\n", D0AB_asymmetry);
    quad4Phi2D0ABAsymmetryErrStr.Form("%6.2f\n", D0AB_asymmetryErr); 

    quad4Phi2D0barABAsymmetryStr.Form("%6.2f\n", D0barAB_asymmetry); 
    quad4Phi2D0barABAsymmetryErrStr.Form("%6.2f\n", D0barAB_asymmetryErr); 

    //std::cout<< D0AB_asymmetry <<std::endl; 
    //std::cout<<D0barAB_asymmetry <<std::endl; 

  weight1 = 1./(D0AB_asymmetryErr*D0AB_asymmetryErr);
    weight2 = 1./(D0barAB_asymmetryErr*D0barAB_asymmetryErr);

  CPV = (D0AB_asymmetry*weight1+D0barAB_asymmetry*weight2)/(weight1+weight2);
  CPVErr = sqrt(1./(weight1+weight2));

  quad4Phi2AsymmetryABStr.Form("%6.4f\n",CPV); 
  quad4Phi2AsymmetryErrABStr.Form("%6.4f\n", CPVErr); 

  //std::cout<<"CPA " << CPV << std::endl; 

  for(int bin = 1; bin < 26; bin++)
    {
        nD02PhiAbove = nD02PhiAbove+quad4Phi2HistAB->GetBinContent(bin); 
        nD0bar2PhiAbove = nD0bar2PhiAbove+quad4Phi2DSHistAB->GetBinContent(bin); 
    }

    for(int bin = 26; bin < 51; bin++)
    {
        nD0bar2PhiBelow = nD0bar2PhiBelow + quad4Phi2DSHistAB->GetBinContent(bin); 
        nD02PhiBelow = nD02PhiBelow + quad4Phi2HistAB->GetBinContent(bin); 
    }

    for(int bin = 51; bin < 76; bin++)
    {
        nD02PhiAbove = nD02PhiAbove+quad4Phi2HistAB->GetBinContent(bin); 
        nD0bar2PhiAbove = nD0bar2PhiAbove+quad4Phi2DSHistAB->GetBinContent(bin); 
    }

    for(int bin = 76; bin <101; bin++)
    {
        nD0bar2PhiBelow = nD0bar2PhiBelow + quad4Phi2DSHistAB->GetBinContent(bin); 
        nD02PhiBelow = nD02PhiBelow + quad4Phi2HistAB->GetBinContent(bin); 
    }

    D0Sin2Phi_asymmetry = (nD02PhiAbove - nD02PhiBelow)/(nD02PhiAbove+nD02PhiBelow); 
    D0Sin2Phi_asymmetryErr = sqrt((1-D0Sin2Phi_asymmetry*D0Sin2Phi_asymmetry)/( nD02PhiAbove + nD02PhiBelow )); 

    D0barSin2Phi_asymmetry = (nD0bar2PhiAbove - nD0bar2PhiBelow)/(nD0bar2PhiAbove+nD0bar2PhiBelow); 
    D0barSin2Phi_asymmetryErr = sqrt((1-D0barSin2Phi_asymmetry*D0barSin2Phi_asymmetry)/( nD0bar2PhiAbove + nD0bar2PhiBelow ));

    weight1 = 1./(D0Sin2Phi_asymmetry*D0Sin2Phi_asymmetry); 
    weight2 = 1./(D0barSin2Phi_asymmetry*D0barSin2Phi_asymmetry); 

    Sin2PhiCP_asymmetry = (D0Sin2Phi_asymmetry*weight1+D0barSin2Phi_asymmetry*weight2)/(weight1+weight2);
    Sin2PhiCP_asymmetryErr = sqrt(1./(weight1+weight2));

    DiffSin2Phi = (D0Sin2Phi_asymmetry*weight1-D0barSin2Phi_asymmetry*weight2)/(weight1+weight2);

    TString D0Sin2Phi_asymmetryStr4; 
    TString D0Sin2Phi_asymmetryErrStr4; 
    TString D0barSin2Phi_asymmetryStr4; 
    TString D0barSin2Phi_asymmetryErrStr4; 
    TString Sin2PhiCP_asymmetryStr4; 
    TString Sin2PhiCP_asymmetryErrStr4; 
    TString DiffSin2PhiStr4;

    D0Sin2Phi_asymmetryStr4.Form("%6.2f\n", D0Sin2Phi_asymmetry); 
    D0Sin2Phi_asymmetryErrStr4.Form("%6.2f\n", D0Sin2Phi_asymmetryErr); 
    D0barSin2Phi_asymmetryStr4.Form("%6.2f\n", D0barSin2Phi_asymmetry);
    D0barSin2Phi_asymmetryErrStr4.Form("%6.2f\n", D0barSin2Phi_asymmetryErr); 
    Sin2PhiCP_asymmetryStr4.Form("%6.4f\n", Sin2PhiCP_asymmetry); 
    Sin2PhiCP_asymmetryErrStr4.Form("%6.4f\n", Sin2PhiCP_asymmetryErr); 
    DiffSin2PhiStr4.Form("%6.2f\n", DiffSin2Phi); 

  quad4Phi2CanvasAB->cd(); 
  quad4Phi2HistAB->GetXaxis()->CenterTitle(); 
  quad4Phi2HistAB->GetYaxis()->CenterTitle(); 
  quad4Phi2HistAB->GetXaxis()->SetTitleOffset(1.25); 
  quad4Phi2HistAB->GetYaxis()->SetTitleOffset(1.35);
  quad4Phi2HistAB->Draw(); 
  quad4Phi2CanvasAB->SaveAs("quad4Phi2HistAB.png"); 

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  quad1Phi2DSCanvasAB->cd(); 
  quad1Phi2DSHistAB->GetXaxis()->CenterTitle(); 
  quad1Phi2DSHistAB->GetYaxis()->CenterTitle(); 
  quad1Phi2DSHistAB->GetXaxis()->SetTitleOffset(1.25); 
  quad1Phi2DSHistAB->GetYaxis()->SetTitleOffset(1.35);
  quad1Phi2DSHistAB->Draw(); 
  quad1Phi2DSCanvasAB->SaveAs("quad1Phi2DSHistAB.png"); 

  

  quad2Phi2DSCanvasAB->cd(); 
  quad2Phi2DSHistAB->GetXaxis()->CenterTitle(); 
  quad2Phi2DSHistAB->GetYaxis()->CenterTitle(); 
  quad2Phi2DSHistAB->GetXaxis()->SetTitleOffset(1.25); 
  quad2Phi2DSHistAB->GetYaxis()->SetTitleOffset(1.35);
  quad2Phi2DSHistAB->Draw(); 
  quad2Phi2DSCanvasAB->SaveAs("quad2Phi2DSHistAB.png"); 

  

  quad3Phi2DSCanvasAB->cd();
  quad3Phi2DSHistAB->GetXaxis()->CenterTitle(); 
  quad3Phi2DSHistAB->GetYaxis()->CenterTitle(); 
  quad3Phi2DSHistAB->GetXaxis()->SetTitleOffset(1.25); 
  quad3Phi2DSHistAB->GetYaxis()->SetTitleOffset(1.35); 
  quad3Phi2DSHistAB->Draw(); 
  quad3Phi2DSCanvasAB->SaveAs("quad3Phi2DSHistAB.png"); 

  

  quad4Phi2DSCanvasAB->cd(); 
  quad4Phi2DSHistAB->GetXaxis()->CenterTitle(); 
  quad4Phi2DSHistAB->GetYaxis()->CenterTitle(); 
  quad4Phi2DSHistAB->GetXaxis()->SetTitleOffset(1.25); 
  quad4Phi2DSHistAB->GetYaxis()->SetTitleOffset(1.35);
  quad4Phi2DSHistAB->Draw(); 
  quad4Phi2DSCanvasAB->SaveAs("quad4Phi2DSHistAB.png"); 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double maxBin; 
  double binContent; 
  TLegend *leg = new TLegend(0.15,0.85,0.3,0.95);
  TPaveText *pt1 = new TPaveText(0.7,0.75,0.98,0.98,"NDC"); // NDC sets coords
  pt1->SetTextSize(0.025); 

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  quad1Phi2CompareABCanvas->cd(); 
  quad1Phi2HistAB->SetLineColor(kBlue); 
  quad1Phi2DSHistAB->SetLineColor(kRed); 
  quad1Phi2HistAB->GetYaxis()->SetLimits(0,1650);
  quad1Phi2HistAB->SetMinimum(0);
  maxBin = quad1Phi2HistAB->GetMaximumBin(); 
  binContent = 1.3*quad1Phi2HistAB->GetBinContent(maxBin); 
  quad1Phi2HistAB->SetMaximum(binContent); 
  TLine *line1A = new TLine(3.14,0,3.14, binContent); 
  leg->SetHeader("Legend"); 
  leg->AddEntry(quad1Phi2HistAB, "Dzero"); 
  leg->AddEntry(quad1Phi2DSHistAB, "Dzero-bar"); 
  pt1->AddText("PV D0:" + quad1Phi2D0ABAsymmetryStr + " #pm " + quad1Phi2D0ABAsymmetryErrStr);
  pt1->AddText("PV D0bar: " + quad1Phi2D0barABAsymmetryStr + " #pm " + quad1Phi2D0barABAsymmetryErrStr); 
  pt1->AddText("CPV: " + quad1Phi2AsymmetryABStr + " #pm" + quad1Phi2AsymmetryErrABStr);
  pt1->AddText("PV D0 Sin(2#phi): " + D0Sin2Phi_asymmetryStr1 + " #pm " + D0Sin2Phi_asymmetryErrStr1); 
  pt1->AddText("PV D0bar Sin(2#phi): " + D0barSin2Phi_asymmetryStr1 + " #pm " + D0barSin2Phi_asymmetryErrStr1); 
  pt1->AddText("CPV Sin(2#phi): " + Sin2PhiCP_asymmetryStr1 + " #pm " + Sin2PhiCP_asymmetryErrStr1); 
  pt1->AddText("D0-D0bar : " + DiffSin2PhiStr1); 
  quad1Phi2HistAB->SetTitle("Quad 1: " + DecayOne);
  quad1Phi2HistAB->Draw(); 
  quad1Phi2DSHistAB->Draw("same");
  line1A->Draw("same");
  leg->Draw("same"); 
  pt1->Draw("same");  
  quad1Phi2CompareABCanvas->SaveAs("quad1Phi2CompAB.png"); 
  pt1->Clear(); 
  leg->Clear(); 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////   

  quad2Phi2CompareABCanvas->cd(); 
  quad2Phi2HistAB->SetLineColor(kBlue); 
  quad2Phi2DSHistAB->SetLineColor(kRed); 
  quad2Phi2HistAB->GetYaxis()->SetLimits(0,1650);
  quad2Phi2HistAB->SetMinimum(0);
  maxBin = quad2Phi2HistAB->GetMaximumBin(); 
  binContent = 1.3*quad2Phi2HistAB->GetBinContent(maxBin); 
  quad2Phi2HistAB->SetMaximum(binContent); 
  TLine *line2A = new TLine(3.14,0,3.14, binContent); 
  leg->SetHeader("Legend"); 
  leg->AddEntry(quad2Phi2HistAB, "Dzero"); 
  leg->AddEntry(quad2Phi2DSHistAB, "Dzero-bar"); 
  pt1->AddText("PV D0:" + quad2Phi2D0ABAsymmetryStr + " #pm " + quad2Phi2D0ABAsymmetryErrStr);
  pt1->AddText("PV D0bar: " + quad2Phi2D0barABAsymmetryStr + " #pm " + quad2Phi2D0barABAsymmetryErrStr); 
  pt1->AddText("CPV: " + quad2Phi2AsymmetryABStr + " #pm" + quad2Phi2AsymmetryErrABStr);
  pt1->AddText("PV D0 Sin(2#phi): " + D0Sin2Phi_asymmetryStr2 + " #pm " + D0Sin2Phi_asymmetryErrStr2); 
  pt1->AddText("PV D0bar Sin(2#phi): " + D0barSin2Phi_asymmetryStr2 + " #pm " + D0barSin2Phi_asymmetryErrStr2); 
  pt1->AddText("CPV Sin(2#phi): " + Sin2PhiCP_asymmetryStr2 + " #pm " + Sin2PhiCP_asymmetryErrStr2); 
  pt1->AddText("D0-D0bar : " + DiffSin2PhiStr2); 
  quad2Phi2HistAB->SetTitle("Quad 2: " + DecayOne);
  quad2Phi2HistAB->Draw(); 
  quad2Phi2DSHistAB->Draw("same");
  line2A->Draw("same");
  leg->Draw("same"); 
  pt1->Draw("same");  
  quad2Phi2CompareABCanvas->SaveAs("quad2Phi2CompAB.png"); 
  pt1->Clear(); 
  leg->Clear(); 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  

  quad3Phi2CompareABCanvas->cd(); 
  quad3Phi2HistAB->SetLineColor(kBlue); 
  quad3Phi2DSHistAB->SetLineColor(kRed); 
  quad3Phi2HistAB->GetYaxis()->SetLimits(0,1650);
  quad3Phi2HistAB->SetMinimum(0);
  maxBin = quad3Phi2HistAB->GetMaximumBin(); 
  binContent = 1.3*quad3Phi2HistAB->GetBinContent(maxBin); 
  quad3Phi2HistAB->SetMaximum(binContent); 
  leg->SetHeader("Legend"); 
  leg->AddEntry(quad3Phi2HistAB, "Dzero"); 
  leg->AddEntry(quad3Phi2DSHistAB, "Dzero-bar"); 
  pt1->AddText("PV D0:" + quad3Phi2D0ABAsymmetryStr + " #pm " + quad3Phi2D0ABAsymmetryErrStr);
  pt1->AddText("PV D0bar: " + quad3Phi2D0barABAsymmetryStr + " #pm " + quad3Phi2D0barABAsymmetryErrStr); 
  pt1->AddText("CPV: " + quad3Phi2AsymmetryABStr + " #pm" + quad3Phi2AsymmetryErrABStr);
  pt1->AddText("PV D0 Sin(2#phi): " + D0Sin2Phi_asymmetryStr3 + " #pm " + D0Sin2Phi_asymmetryErrStr3); 
  pt1->AddText("PV D0bar Sin(2#phi): " + D0barSin2Phi_asymmetryStr3 + " #pm " + D0barSin2Phi_asymmetryErrStr3); 
  pt1->AddText("CPV Sin(2#phi): " + Sin2PhiCP_asymmetryStr3 + " #pm " + Sin2PhiCP_asymmetryErrStr3); 
  pt1->AddText("D0-D0bar : " + DiffSin2PhiStr3); 
  TLine *line3A = new TLine(3.14,0,3.14, binContent); 
  quad3Phi2HistAB->SetTitle("Quad 3: " + DecayOne);
  quad3Phi2HistAB->Draw(); 
  quad3Phi2DSHistAB->Draw("same");
  line3A->Draw("same");
  leg->Draw("same"); 
  pt1->Draw("same");  
  quad3Phi2CompareABCanvas->SaveAs("quad3Phi2CompAB.png"); 
  leg->Clear("same"); 
  pt1->Clear("same");  
  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  quad4Phi2CompareABCanvas->cd(); 
  quad4Phi2HistAB->SetLineColor(kBlue); 
  quad4Phi2DSHistAB->SetLineColor(kRed); 
  quad4Phi2HistAB->GetYaxis()->SetLimits(0,1650);
  quad4Phi2HistAB->SetMinimum(0);
  maxBin = quad4Phi2HistAB->GetMaximumBin(); 
  binContent = 1.3*quad4Phi2HistAB->GetBinContent(maxBin); 
  quad4Phi2HistAB->SetMaximum(binContent); 
  TLine *line4A = new TLine(3.14,0,3.14, binContent); 
  leg->SetHeader("Legend"); 
  leg->AddEntry(quad4Phi2HistAB, "Dzero"); 
  leg->AddEntry(quad4Phi2DSHistAB, "Dzero-bar"); 
  pt1->AddText("PV D0:" + quad4Phi2D0ABAsymmetryStr + " #pm " + quad4Phi2D0ABAsymmetryErrStr);
  pt1->AddText("PV D0bar: " + quad4Phi2D0barABAsymmetryStr + " #pm " + quad4Phi2D0barABAsymmetryErrStr); 
  pt1->AddText("CPV: " + quad4Phi2AsymmetryABStr + " #pm" + quad4Phi2AsymmetryErrABStr);
  pt1->AddText("PV D0 Sin(2#phi): " + D0Sin2Phi_asymmetryStr4 + " #pm " + D0Sin2Phi_asymmetryErrStr4); 
  pt1->AddText("PV D0bar Sin(2#phi): " + D0barSin2Phi_asymmetryStr4 + " #pm " + D0barSin2Phi_asymmetryErrStr4); 
  pt1->AddText("CPV Sin(2#phi): " + Sin2PhiCP_asymmetryStr4 + " #pm " + Sin2PhiCP_asymmetryErrStr4); 
  pt1->AddText("D0-D0bar : " + DiffSin2PhiStr4); 
  quad4Phi2HistAB->SetTitle("Quad 4: " + DecayOne);
  quad4Phi2HistAB->Draw(); 
  quad4Phi2DSHistAB->Draw("same");
  line4A->Draw("same");
  leg->Draw("same"); 
  pt1->Draw("same");  
  quad4Phi2CompareABCanvas->SaveAs("quad4Phi2CompAB.png"); 
  leg->Clear("same"); 
  pt1->Clear("same");  

  
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    nD0 = quad1HistAD->GetEntries();
    nD0bar = quad1DSHistAD->GetEntries();
    scale = nD0/nD0bar;
    nD0Above = 0;
    nD0barAbove = 0;

    for (int bin=0; bin<51; bin++) 
    {
      nD0Above = nD0Above+quad1HistAD->GetBinContent(bin);
      nD0barAbove = nD0barAbove+quad1DSHistAD->GetBinContent(bin);
    }
   
    nD0Below = 0;
    nD0barBelow = 0;
    for (int bin=51; bin<101; bin++)
    {
      nD0Below = nD0Below+quad1HistAD->GetBinContent(bin);
      nD0barBelow = nD0barBelow+quad1DSHistAD->GetBinContent(bin);
    }

    Double_t D0AD_asymmetry = (nD0Above-nD0Below)/(nD0Above+nD0Below);
    Double_t D0AD_asymmetryErr = sqrt( (1.0 - D0AD_asymmetry*D0AD_asymmetry)/(nD0Above+nD0Below) );

    Double_t D0barAD_asymmetry = (nD0barAbove-nD0barBelow)/(nD0barAbove+nD0barBelow);
    Double_t D0barAD_asymmetryErr = sqrt( (1.0 - D0barAD_asymmetry*D0barAD_asymmetry)/(nD0barAbove+nD0barBelow) );

    quad1D0ADAsymmetryStr.Form("%6.2f\n", D0AD_asymmetry);
    quad1D0ADAsymmetryErrStr.Form("%6.2f\n", D0AD_asymmetryErr); 

    quad1D0barADAsymmetryStr.Form("%6.2f\n", D0barAD_asymmetry); 
    quad1D0barADAsymmetryErrStr.Form("%6.2f\n", D0barAD_asymmetryErr); 

    //std::cout<< D0AD_asymmetry <<std::endl; 
    //std::cout<<D0barAD_asymmetry <<std::endl; 

    weight1 = 1./(D0AD_asymmetryErr*D0AD_asymmetryErr);
    weight2 = 1./(D0barAD_asymmetryErr*D0barAD_asymmetryErr);

  CPV = (D0AD_asymmetry*weight1+D0barAD_asymmetry*weight2)/(weight1+weight2);
  CPVErr = sqrt(1./(weight1+weight2));

  quad1AsymmetryADStr.Form("%6.4f\n",CPV); 
  quad1AsymmetryErrADStr.Form("%6.4f\n", CPVErr); 

  //std::cout<<"CPA " << CPV << std::endl; 

  for(int bin = 1; bin < 26; bin++)
    {
        nD02PhiAbove = nD02PhiAbove+quad1HistAD->GetBinContent(bin); 
        nD0bar2PhiAbove = nD0bar2PhiAbove+quad1DSHistAD->GetBinContent(bin); 
    }

    for(int bin = 26; bin < 51; bin++)
    {
        nD0bar2PhiBelow = nD0bar2PhiBelow + quad1DSHistAD->GetBinContent(bin); 
        nD02PhiBelow = nD02PhiBelow + quad1HistAD->GetBinContent(bin); 
    }

    for(int bin = 51; bin < 76; bin++)
    {
        nD02PhiAbove = nD02PhiAbove+quad1HistAD->GetBinContent(bin); 
        nD0bar2PhiAbove = nD0bar2PhiAbove+quad1DSHistAD->GetBinContent(bin); 
    }

    for(int bin = 76; bin <101; bin++)
    {
        nD0bar2PhiBelow = nD0bar2PhiBelow + quad1DSHistAD->GetBinContent(bin); 
        nD02PhiBelow = nD02PhiBelow + quad1HistAD->GetBinContent(bin); 
    }

    

    D0Sin2Phi_asymmetry = (nD02PhiAbove - nD02PhiBelow)/(nD02PhiAbove+nD02PhiBelow); 
    D0Sin2Phi_asymmetryErr = sqrt((1-D0Sin2Phi_asymmetry*D0Sin2Phi_asymmetry)/( nD02PhiAbove + nD02PhiBelow )); 

    D0barSin2Phi_asymmetry = (nD0bar2PhiAbove - nD0bar2PhiBelow)/(nD0bar2PhiAbove+nD0bar2PhiBelow); 
    D0barSin2Phi_asymmetryErr = sqrt((1-D0barSin2Phi_asymmetry*D0barSin2Phi_asymmetry)/( nD0bar2PhiAbove + nD0bar2PhiBelow ));

    weight1 = 1./(D0Sin2Phi_asymmetry*D0Sin2Phi_asymmetry); 
    weight2 = 1./(D0barSin2Phi_asymmetry*D0barSin2Phi_asymmetry); 

    Sin2PhiCP_asymmetry = (D0Sin2Phi_asymmetry*weight1+D0barSin2Phi_asymmetry*weight2)/(weight1+weight2);
    Sin2PhiCP_asymmetryErr = sqrt(1./(weight1+weight2));

    DiffSin2Phi = (D0Sin2Phi_asymmetry*weight1-D0barSin2Phi_asymmetry*weight2)/(weight1+weight2);

    D0Sin2Phi_asymmetryStr1; 
    D0Sin2Phi_asymmetryErrStr1; 
    D0barSin2Phi_asymmetryStr1; 
    D0barSin2Phi_asymmetryErrStr1; 
    Sin2PhiCP_asymmetryStr1; 
    Sin2PhiCP_asymmetryErrStr1; 
    DiffSin2PhiStr1; 

    D0Sin2Phi_asymmetryStr1.Form("%6.2f\n", D0Sin2Phi_asymmetry); 
    D0Sin2Phi_asymmetryErrStr1.Form("%6.2f\n", D0Sin2Phi_asymmetryErr); 
    D0barSin2Phi_asymmetryStr1.Form("%6.2f\n", D0barSin2Phi_asymmetry);
    D0barSin2Phi_asymmetryErrStr1.Form("%6.2f\n", D0barSin2Phi_asymmetryErr); 
    Sin2PhiCP_asymmetryStr1.Form("%6.4f\n", Sin2PhiCP_asymmetry); 
    Sin2PhiCP_asymmetryErrStr1.Form("%6.4f\n", Sin2PhiCP_asymmetryErr); 
    DiffSin2PhiStr1.Form("%6.2f\n", DiffSin2Phi); 
      

  quad1CanvasAD->cd(); 
  quad1HistAD->GetXaxis()->CenterTitle(); 
  quad1HistAD->GetYaxis()->CenterTitle(); 
  quad1HistAD->GetXaxis()->SetTitleOffset(1.25); 
  quad1HistAD->GetYaxis()->SetTitleOffset(1.35);
  quad1HistAD->Draw(); 
  quad1CanvasAD->SaveAs("quad1HistAD.png"); 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  

    nD0 = quad2HistAD->GetEntries();
    nD0bar = quad2DSHistAD->GetEntries();
      
    scale = nD0/nD0bar;
    nD0Above = 0;
    nD0barAbove = 0;

    for (int bin=0; bin<51; bin++) 
    {
        nD0Above = nD0Above+quad2HistAD->GetBinContent(bin);
        nD0barAbove = nD0barAbove+quad2DSHistAD->GetBinContent(bin);
    }
   
    nD0Below = 0;
    nD0barBelow = 0;
    for (int bin=51; bin<101; bin++)
    {
        nD0Below = nD0Below+quad2HistAD->GetBinContent(bin);
        nD0barBelow = nD0barBelow+quad2DSHistAD->GetBinContent(bin);
    }

    D0AD_asymmetry = (nD0Above-nD0Below)/(nD0Above+nD0Below);
    D0AD_asymmetryErr = sqrt( (1.0 - D0AD_asymmetry*D0AD_asymmetry)/(nD0Above+nD0Below) );

    D0barAD_asymmetry = (nD0barAbove-nD0barBelow)/(nD0barAbove+nD0barBelow);
    D0barAD_asymmetryErr = sqrt( (1.0 - D0barAD_asymmetry*D0barAD_asymmetry)/(nD0barAbove+nD0barBelow) );

    quad2D0ADAsymmetryStr.Form("%6.2f\n", D0AD_asymmetry);
    quad2D0ADAsymmetryErrStr.Form("%6.2f\n", D0AD_asymmetryErr); 

    quad2D0barADAsymmetryStr.Form("%6.2f\n", D0barAD_asymmetry); 
    quad2D0barADAsymmetryErrStr.Form("%6.2f\n", D0barAD_asymmetryErr); 

    //std::cout<< D0AD_asymmetry <<std::endl; 
    //std::cout<<D0barAD_asymmetry <<std::endl; 

    weight1 = 1./(D0AD_asymmetryErr*D0AD_asymmetryErr);
    weight2 = 1./(D0barAD_asymmetryErr*D0barAD_asymmetryErr);

  CPV = (D0AD_asymmetry*weight1+D0barAD_asymmetry*weight2)/(weight1+weight2);
  CPVErr = sqrt(1./(weight1+weight2));

  //std::cout<<"CPA " << CPV << std::endl; 

  quad2AsymmetryADStr.Form("%6.4f\n",CPV); 
  quad2AsymmetryErrADStr.Form("%6.4f\n", CPVErr); 

  for(int bin = 1; bin < 26; bin++)
    {
        nD02PhiAbove = nD02PhiAbove+quad2HistAD->GetBinContent(bin); 
        nD0bar2PhiAbove = nD0bar2PhiAbove+quad2DSHistAD->GetBinContent(bin); 
    }

    for(int bin = 26; bin < 51; bin++)
    {
        nD0bar2PhiBelow = nD0bar2PhiBelow + quad2DSHistAD->GetBinContent(bin); 
        nD02PhiBelow = nD02PhiBelow + quad2HistAD->GetBinContent(bin); 
    }

    for(int bin = 51; bin < 76; bin++)
    {
        nD02PhiAbove = nD02PhiAbove+quad2HistAD->GetBinContent(bin); 
        nD0bar2PhiAbove = nD0bar2PhiAbove+quad2DSHistAD->GetBinContent(bin); 
    }

    for(int bin = 76; bin <101; bin++)
    {
        nD0bar2PhiBelow = nD0bar2PhiBelow + quad2DSHistAD->GetBinContent(bin); 
        nD02PhiBelow = nD02PhiBelow + quad2HistAD->GetBinContent(bin); 
    }

    

    D0Sin2Phi_asymmetry = (nD02PhiAbove - nD02PhiBelow)/(nD02PhiAbove+nD02PhiBelow); 
    D0Sin2Phi_asymmetryErr = sqrt((1-D0Sin2Phi_asymmetry*D0Sin2Phi_asymmetry)/( nD02PhiAbove + nD02PhiBelow )); 

    D0barSin2Phi_asymmetry = (nD0bar2PhiAbove - nD0bar2PhiBelow)/(nD0bar2PhiAbove+nD0bar2PhiBelow); 
    D0barSin2Phi_asymmetryErr = sqrt((1-D0barSin2Phi_asymmetry*D0barSin2Phi_asymmetry)/( nD0bar2PhiAbove + nD0bar2PhiBelow ));

    weight1 = 1./(D0Sin2Phi_asymmetry*D0Sin2Phi_asymmetry); 
    weight2 = 1./(D0barSin2Phi_asymmetry*D0barSin2Phi_asymmetry); 

    Sin2PhiCP_asymmetry = (D0Sin2Phi_asymmetry*weight1+D0barSin2Phi_asymmetry*weight2)/(weight1+weight2);
    Sin2PhiCP_asymmetryErr = sqrt(1./(weight1+weight2));

    DiffSin2Phi = (D0Sin2Phi_asymmetry*weight1-D0barSin2Phi_asymmetry*weight2)/(weight1+weight2);

    D0Sin2Phi_asymmetryStr2; 
    D0Sin2Phi_asymmetryErrStr2; 
    D0barSin2Phi_asymmetryStr2; 
    D0barSin2Phi_asymmetryErrStr2; 
    Sin2PhiCP_asymmetryStr2; 
    Sin2PhiCP_asymmetryErrStr2; 
    DiffSin2PhiStr2; 

    D0Sin2Phi_asymmetryStr2.Form("%6.2f\n", D0Sin2Phi_asymmetry); 
    D0Sin2Phi_asymmetryErrStr2.Form("%6.2f\n", D0Sin2Phi_asymmetryErr); 
    D0barSin2Phi_asymmetryStr2.Form("%6.2f\n", D0barSin2Phi_asymmetry);
    D0barSin2Phi_asymmetryErrStr2.Form("%6.2f\n", D0barSin2Phi_asymmetryErr); 
    Sin2PhiCP_asymmetryStr2.Form("%6.4f\n", Sin2PhiCP_asymmetry); 
    Sin2PhiCP_asymmetryErrStr2.Form("%6.4f\n", Sin2PhiCP_asymmetryErr); 
    DiffSin2PhiStr2.Form("%6.2f\n", DiffSin2Phi); 

  quad2CanvasAD->cd(); 
  quad2HistAD->GetXaxis()->CenterTitle(); 
  quad2HistAD->GetYaxis()->CenterTitle(); 
  quad2HistAD->GetXaxis()->SetTitleOffset(1.25); 
  quad2HistAD->GetYaxis()->SetTitleOffset(1.35);
  quad2HistAD->Draw(); 
  quad2CanvasAD->SaveAs("quad2HistAD.png"); 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    nD0 = quad3HistAD->GetEntries();
    nD0bar = quad3DSHistAD->GetEntries();
      
    scale = nD0/nD0bar;
    nD0Above = 0;
    nD0barAbove = 0;

    for (int bin=0; bin<51; bin++) 
    {
      nD0Above = nD0Above+quad3HistAD->GetBinContent(bin);
      nD0barAbove = nD0barAbove+quad3DSHistAD->GetBinContent(bin);
  }
   
  nD0Below = 0;
    nD0barBelow = 0;
    for (int bin=51; bin<101; bin++)
    {
        nD0Below = nD0Below+quad3HistAD->GetBinContent(bin);
      nD0barBelow = nD0barBelow+quad3DSHistAD->GetBinContent(bin);
    }

    D0AD_asymmetry = (nD0Above-nD0Below)/(nD0Above+nD0Below);
    D0AD_asymmetryErr = sqrt( (1.0 - D0AD_asymmetry*D0AD_asymmetry)/(nD0Above+nD0Below) );

    D0barAD_asymmetry = (nD0barAbove-nD0barBelow)/(nD0barAbove+nD0barBelow);
    D0barAD_asymmetryErr = sqrt( (1.0 - D0barAD_asymmetry*D0barAD_asymmetry)/(nD0barAbove+nD0barBelow) );

    quad3D0ADAsymmetryStr.Form("%6.2f\n", D0AD_asymmetry);
    quad3D0ADAsymmetryErrStr.Form("%6.2f\n", D0AD_asymmetryErr); 

    quad3D0barADAsymmetryStr.Form("%6.2f\n", D0barAD_asymmetry); 
    quad3D0barADAsymmetryErrStr.Form("%6.2f\n", D0barAD_asymmetryErr); 

    //std::cout<< D0AD_asymmetry <<std::endl; 
    //std::cout<<D0barAD_asymmetry <<std::endl; 

  weight1 = 1./(D0AD_asymmetryErr*D0AD_asymmetryErr);
    weight2 = 1./(D0barAD_asymmetryErr*D0barAD_asymmetryErr);

  CPV = (D0AD_asymmetry*weight1+D0barAD_asymmetry*weight2)/(weight1+weight2);
  CPVErr = sqrt(1./(weight1+weight2));

  quad3AsymmetryADStr.Form("%6.4f\n",CPV); 
  quad3AsymmetryErrADStr.Form("%6.4f\n", CPVErr); 

  //std::cout<<"CPA " << CPV << std::endl; 

  for(int bin = 1; bin < 26; bin++)
    {
        nD02PhiAbove = nD02PhiAbove+quad3HistAD->GetBinContent(bin); 
        nD0bar2PhiAbove = nD0bar2PhiAbove+quad3DSHistAD->GetBinContent(bin); 
    }

    for(int bin = 26; bin < 51; bin++)
    {
        nD0bar2PhiBelow = nD0bar2PhiBelow + quad3DSHistAD->GetBinContent(bin); 
        nD02PhiBelow = nD02PhiBelow + quad3HistAD->GetBinContent(bin); 
    }

    for(int bin = 51; bin < 76; bin++)
    {
        nD02PhiAbove = nD02PhiAbove+quad3HistAD->GetBinContent(bin); 
        nD0bar2PhiAbove = nD0bar2PhiAbove+quad3DSHistAD->GetBinContent(bin); 
    }

    for(int bin = 76; bin <101; bin++)
    {
        nD0bar2PhiBelow = nD0bar2PhiBelow + quad3DSHistAD->GetBinContent(bin); 
        nD02PhiBelow = nD02PhiBelow + quad3HistAD->GetBinContent(bin); 
    }

    

    D0Sin2Phi_asymmetry = (nD02PhiAbove - nD02PhiBelow)/(nD02PhiAbove+nD02PhiBelow); 
    D0Sin2Phi_asymmetryErr = sqrt((1-D0Sin2Phi_asymmetry*D0Sin2Phi_asymmetry)/( nD02PhiAbove + nD02PhiBelow )); 

    D0barSin2Phi_asymmetry = (nD0bar2PhiAbove - nD0bar2PhiBelow)/(nD0bar2PhiAbove+nD0bar2PhiBelow); 
    D0barSin2Phi_asymmetryErr = sqrt((1-D0barSin2Phi_asymmetry*D0barSin2Phi_asymmetry)/( nD0bar2PhiAbove + nD0bar2PhiBelow ));

    weight1 = 1./(D0Sin2Phi_asymmetry*D0Sin2Phi_asymmetry); 
    weight2 = 1./(D0barSin2Phi_asymmetry*D0barSin2Phi_asymmetry); 

    Sin2PhiCP_asymmetry = (D0Sin2Phi_asymmetry*weight1+D0barSin2Phi_asymmetry*weight2)/(weight1+weight2);
    Sin2PhiCP_asymmetryErr = sqrt(1./(weight1+weight2));

    DiffSin2Phi = (D0Sin2Phi_asymmetry*weight1-D0barSin2Phi_asymmetry*weight2)/(weight1+weight2);

    D0Sin2Phi_asymmetryStr3; 
    D0Sin2Phi_asymmetryErrStr3; 
    D0barSin2Phi_asymmetryStr3; 
    D0barSin2Phi_asymmetryErrStr3; 
    Sin2PhiCP_asymmetryStr3; 
    Sin2PhiCP_asymmetryErrStr3; 
    DiffSin2PhiStr3; 

    D0Sin2Phi_asymmetryStr3.Form("%6.2f\n", D0Sin2Phi_asymmetry); 
    D0Sin2Phi_asymmetryErrStr3.Form("%6.2f\n", D0Sin2Phi_asymmetryErr); 
    D0barSin2Phi_asymmetryStr3.Form("%6.2f\n", D0barSin2Phi_asymmetry);
    D0barSin2Phi_asymmetryErrStr3.Form("%6.2f\n", D0barSin2Phi_asymmetryErr); 
    Sin2PhiCP_asymmetryStr3.Form("%6.4f\n", Sin2PhiCP_asymmetry); 
    Sin2PhiCP_asymmetryErrStr3.Form("%6.4f\n", Sin2PhiCP_asymmetryErr); 
    DiffSin2PhiStr3.Form("%6.2f\n", DiffSin2Phi); 

  quad3CanvasAD->cd(); 
  quad3HistAD->GetXaxis()->CenterTitle(); 
  quad3HistAD->GetYaxis()->CenterTitle(); 
  quad3HistAD->GetXaxis()->SetTitleOffset(1.25); 
  quad3HistAD->GetYaxis()->SetTitleOffset(1.35);
  quad3HistAD->Draw(); 
  quad3CanvasAD->SaveAs("quad3HistAD.png"); 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  nD0 = quad4HistAD->GetEntries();
    nD0bar = quad4DSHistAD->GetEntries();
  scale = nD0/nD0bar;
    nD0Above = 0;
    nD0barAbove = 0;

    for (int bin=0; bin<51; bin++) 
    {
      nD0Above = nD0Above+quad4HistAD->GetBinContent(bin);
      nD0barAbove = nD0barAbove+quad4DSHistAD->GetBinContent(bin);
  }
   
  nD0Below = 0;
    nD0barBelow = 0;
    for (int bin=51; bin<101; bin++)
    {
      nD0Below = nD0Below+quad4HistAD->GetBinContent(bin);
      nD0barBelow = nD0barBelow+quad4DSHistAD->GetBinContent(bin);
    }

    D0AD_asymmetry = (nD0Above-nD0Below)/(nD0Above+nD0Below);
    D0AD_asymmetryErr = sqrt( (1.0 - D0AD_asymmetry*D0AD_asymmetry)/(nD0Above+nD0Below) );

    D0barAD_asymmetry = (nD0barAbove-nD0barBelow)/(nD0barAbove+nD0barBelow);
    D0barAD_asymmetryErr = sqrt( (1.0 - D0barAD_asymmetry*D0barAD_asymmetry)/(nD0barAbove+nD0barBelow) );

    quad4D0ADAsymmetryStr.Form("%6.2f\n", D0AD_asymmetry);
    quad4D0ADAsymmetryErrStr.Form("%6.2f\n", D0AD_asymmetryErr); 

    quad4D0barADAsymmetryStr.Form("%6.2f\n", D0barAD_asymmetry); 
    quad4D0barADAsymmetryErrStr.Form("%6.2f\n", D0barAD_asymmetryErr); 

    //std::cout<< D0AD_asymmetry <<std::endl; 
    //std::cout<<D0barAD_asymmetry <<std::endl; 

  weight1 = 1./(D0AD_asymmetryErr*D0AD_asymmetryErr);
    weight2 = 1./(D0barAD_asymmetryErr*D0barAD_asymmetryErr);

  CPV = (D0AD_asymmetry*weight1+D0barAD_asymmetry*weight2)/(weight1+weight2);
  CPVErr = sqrt(1./(weight1+weight2));

  quad4AsymmetryADStr.Form("%6.4f\n",CPV); 
  quad4AsymmetryErrADStr.Form("%6.4f\n", CPVErr); 

  //std::cout<<"CPA " << CPV << std::endl; 

  for(int bin = 1; bin < 26; bin++)
    {
        nD02PhiAbove = nD02PhiAbove+quad4HistAD->GetBinContent(bin); 
        nD0bar2PhiAbove = nD0bar2PhiAbove+quad4DSHistAD->GetBinContent(bin); 
    }

    for(int bin = 26; bin < 51; bin++)
    {
        nD0bar2PhiBelow = nD0bar2PhiBelow + quad4DSHistAD->GetBinContent(bin); 
        nD02PhiBelow = nD02PhiBelow + quad4HistAD->GetBinContent(bin); 
    }

    for(int bin = 51; bin < 76; bin++)
    {
        nD02PhiAbove = nD02PhiAbove+quad4HistAD->GetBinContent(bin); 
        nD0bar2PhiAbove = nD0bar2PhiAbove+quad4DSHistAD->GetBinContent(bin); 
    }

    for(int bin = 76; bin <101; bin++)
    {
        nD0bar2PhiBelow = nD0bar2PhiBelow + quad4DSHistAD->GetBinContent(bin); 
        nD02PhiBelow = nD02PhiBelow + quad4HistAD->GetBinContent(bin); 
    }

    

    D0Sin2Phi_asymmetry = (nD02PhiAbove - nD02PhiBelow)/(nD02PhiAbove+nD02PhiBelow); 
    D0Sin2Phi_asymmetryErr = sqrt((1-D0Sin2Phi_asymmetry*D0Sin2Phi_asymmetry)/( nD02PhiAbove + nD02PhiBelow )); 

    D0barSin2Phi_asymmetry = (nD0bar2PhiAbove - nD0bar2PhiBelow)/(nD0bar2PhiAbove+nD0bar2PhiBelow); 
    D0barSin2Phi_asymmetryErr = sqrt((1-D0barSin2Phi_asymmetry*D0barSin2Phi_asymmetry)/( nD0bar2PhiAbove + nD0bar2PhiBelow ));

    weight1 = 1./(D0Sin2Phi_asymmetry*D0Sin2Phi_asymmetry); 
    weight2 = 1./(D0barSin2Phi_asymmetry*D0barSin2Phi_asymmetry); 

    Sin2PhiCP_asymmetry = (D0Sin2Phi_asymmetry*weight1+D0barSin2Phi_asymmetry*weight2)/(weight1+weight2);
    Sin2PhiCP_asymmetryErr = sqrt(1./(weight1+weight2));

    DiffSin2Phi = (D0Sin2Phi_asymmetry*weight1-D0barSin2Phi_asymmetry*weight2)/(weight1+weight2);

    D0Sin2Phi_asymmetryStr4; 
    D0Sin2Phi_asymmetryErrStr4; 
    D0barSin2Phi_asymmetryStr4; 
    D0barSin2Phi_asymmetryErrStr4; 
    Sin2PhiCP_asymmetryStr4; 
    Sin2PhiCP_asymmetryErrStr4; 
    DiffSin2PhiStr4; 

    D0Sin2Phi_asymmetryStr4.Form("%6.2f\n", D0Sin2Phi_asymmetry); 
    D0Sin2Phi_asymmetryErrStr4.Form("%6.2f\n", D0Sin2Phi_asymmetryErr); 
    D0barSin2Phi_asymmetryStr4.Form("%6.2f\n", D0barSin2Phi_asymmetry);
    D0barSin2Phi_asymmetryErrStr4.Form("%6.2f\n", D0barSin2Phi_asymmetryErr); 
    Sin2PhiCP_asymmetryStr4.Form("%6.4f\n", Sin2PhiCP_asymmetry); 
    Sin2PhiCP_asymmetryErrStr4.Form("%6.4f\n", Sin2PhiCP_asymmetryErr); 
    DiffSin2PhiStr4.Form("%6.2f\n", DiffSin2Phi); 

  quad4CanvasAD->cd(); 
  quad4HistAD->GetXaxis()->CenterTitle(); 
  quad4HistAD->GetYaxis()->CenterTitle(); 
  quad4HistAD->GetXaxis()->SetTitleOffset(1.25); 
  quad4HistAD->GetYaxis()->SetTitleOffset(1.35);
  quad4HistAD->Draw(); 
  quad4CanvasAD->SaveAs("quad4HistAD.png"); 

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  quad1DSCanvasAD->cd(); 
  quad1DSHistAD->GetXaxis()->CenterTitle(); 
  quad1DSHistAD->GetYaxis()->CenterTitle(); 
  quad1DSHistAD->GetXaxis()->SetTitleOffset(1.25); 
  quad1DSHistAD->GetYaxis()->SetTitleOffset(1.35);
  quad1DSHistAD->Draw(); 
  quad1DSCanvasAD->SaveAs("quad1DSHistAD.png"); 

  

  quad2DSCanvasAD->cd(); 
  quad2DSHistAD->GetXaxis()->CenterTitle(); 
  quad2DSHistAD->GetYaxis()->CenterTitle(); 
  quad2DSHistAD->GetXaxis()->SetTitleOffset(1.25); 
  quad2DSHistAD->GetYaxis()->SetTitleOffset(1.35);
  quad2DSHistAD->Draw(); 
  quad2DSCanvasAD->SaveAs("quad2DSHistAD.png"); 

  

  quad3DSCanvasAD->cd();
  quad3DSHistAD->GetXaxis()->CenterTitle(); 
  quad3DSHistAD->GetYaxis()->CenterTitle(); 
  quad3DSHistAD->GetXaxis()->SetTitleOffset(1.25); 
  quad3DSHistAD->GetYaxis()->SetTitleOffset(1.35); 
  quad3DSHistAD->Draw(); 
  quad3DSCanvasAD->SaveAs("quad3DSHistAD.png"); 

  

  quad4DSCanvasAD->cd(); 
  quad4DSHistAD->GetXaxis()->CenterTitle(); 
  quad4DSHistAD->GetYaxis()->CenterTitle(); 
  quad4DSHistAD->GetXaxis()->SetTitleOffset(1.25); 
  quad4DSHistAD->GetYaxis()->SetTitleOffset(1.35);
  quad4DSHistAD->Draw(); 
  quad4DSCanvasAD->SaveAs("quad4DSHistAD.png"); 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  

  quad1CompareADCanvas->cd(); 
  quad1HistAD->SetLineColor(kBlue); 
  quad1DSHistAD->SetLineColor(kRed); 
  quad1HistAD->GetYaxis()->SetLimits(0,1650);
  quad1HistAD->SetMinimum(0);
  maxBin = quad1HistAD->GetMaximumBin(); 
  binContent = 1.3*quad1HistAD->GetBinContent(maxBin); 
  quad1HistAD->SetMaximum(binContent); 
  TLine *line1 = new TLine(3.14,0,3.14, binContent); 
  leg->SetHeader("Legend"); 
  leg->AddEntry(quad1HistAD, "Dzero"); 
  leg->AddEntry(quad1DSHistAD, "Dzero-bar"); 
  pt1->AddText("PV D0:" + quad1D0ADAsymmetryStr + " #pm " + quad1D0ADAsymmetryErrStr);
  pt1->AddText("PV D0bar: " + quad1D0barADAsymmetryStr + " #pm " + quad1D0barADAsymmetryErrStr); 
  pt1->AddText("CPV: " + quad1AsymmetryADStr + " #pm" + quad1AsymmetryErrADStr);
  pt1->AddText("PV D0 Sin(2#phi): " + D0Sin2Phi_asymmetryStr1 + " #pm " + D0Sin2Phi_asymmetryErrStr1); 
  pt1->AddText("PV D0bar Sin(2#phi): " + D0barSin2Phi_asymmetryStr1 + " #pm " + D0barSin2Phi_asymmetryErrStr1); 
  pt1->AddText("CPV Sin(2#phi): " + Sin2PhiCP_asymmetryStr1 + " #pm " + Sin2PhiCP_asymmetryErrStr1); 
  pt1->AddText("D0-D0bar : " + DiffSin2PhiStr1); 
  quad1HistAD->SetTitle("Quad 1: " + DecayTwo);
  quad1HistAD->Draw(); 
  quad1DSHistAD->Draw("same");
  line1->Draw("same");
  leg->Draw("same"); 
  pt1->Draw("same");  
  quad1CompareADCanvas->SaveAs("Quad1CompAD.png"); 
  pt1->Clear(); 
  leg->Clear(); 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////   

  quad2CompareADCanvas->cd(); 
  quad2HistAD->SetLineColor(kBlue); 
  quad2DSHistAD->SetLineColor(kRed); 
  quad2HistAD->GetYaxis()->SetLimits(0,1650);
  quad2HistAD->SetMinimum(0);
  maxBin = quad2HistAD->GetMaximumBin(); 
  binContent = 1.3*quad2HistAD->GetBinContent(maxBin); 
  quad2HistAD->SetMaximum(binContent); 
  TLine *line2 = new TLine(3.14,0,3.14, binContent); 
  leg->SetHeader("Legend"); 
  leg->AddEntry(quad2HistAD, "Dzero"); 
  leg->AddEntry(quad2DSHistAD, "Dzero-bar"); 
  pt1->AddText("PV D0:" + quad2D0ADAsymmetryStr + " #pm " + quad2D0ADAsymmetryErrStr);
  pt1->AddText("PV D0bar: " + quad2D0barADAsymmetryStr + " #pm " + quad2D0barADAsymmetryErrStr); 
  pt1->AddText("CPV: " + quad2AsymmetryADStr + " #pm" + quad2AsymmetryErrADStr);
  pt1->AddText("PV D0 Sin(2#phi): " + D0Sin2Phi_asymmetryStr2 + " #pm " + D0Sin2Phi_asymmetryErrStr2); 
    pt1->AddText("PV D0bar Sin(2#phi): " + D0barSin2Phi_asymmetryStr2 + " #pm " + D0barSin2Phi_asymmetryErrStr2); 
    pt1->AddText("CPV Sin(2#phi): " + Sin2PhiCP_asymmetryStr2 + " #pm " + Sin2PhiCP_asymmetryErrStr2); 
    pt1->AddText("D0-D0bar : " + DiffSin2PhiStr2); 
    quad2HistAD->SetTitle("Quad 2: " + DecayTwo);
  quad2HistAD->Draw(); 
  quad2DSHistAD->Draw("same");
  line2->Draw("same");
  leg->Draw("same"); 
  pt1->Draw("same");  
  quad2CompareADCanvas->SaveAs("Quad2CompAD.png"); 
  pt1->Clear(); 
  leg->Clear(); 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  

  quad3CompareADCanvas->cd(); 
  quad3HistAD->SetLineColor(kBlue); 
  quad3DSHistAD->SetLineColor(kRed); 
  quad3HistAD->GetYaxis()->SetLimits(0,1650);
  quad3HistAD->SetMinimum(0);
  maxBin = quad3HistAD->GetMaximumBin(); 
  binContent = 1.3*quad3HistAD->GetBinContent(maxBin); 
  quad3HistAD->SetMaximum(binContent); 
  leg->SetHeader("Legend"); 
  leg->AddEntry(quad3HistAD, "Dzero"); 
  leg->AddEntry(quad3DSHistAD, "Dzero-bar"); 
  pt1->AddText("PV D0:" + quad3D0ADAsymmetryStr + " #pm " + quad3D0ADAsymmetryErrStr);
  pt1->AddText("PV D0bar: " + quad3D0barADAsymmetryStr + " #pm " + quad3D0barADAsymmetryErrStr); 
  pt1->AddText("CPV: " + quad3AsymmetryADStr + " #pm" + quad3AsymmetryErrADStr);
  pt1->AddText("PV D0 Sin(2#phi): " + D0Sin2Phi_asymmetryStr3 + " #pm " + D0Sin2Phi_asymmetryErrStr3); 
    pt1->AddText("PV D0bar Sin(2#phi): " + D0barSin2Phi_asymmetryStr3 + " #pm " + D0barSin2Phi_asymmetryErrStr3); 
    pt1->AddText("CPV Sin(2#phi): " + Sin2PhiCP_asymmetryStr3 + " #pm " + Sin2PhiCP_asymmetryErrStr3); 
    pt1->AddText("D0-D0bar : " + DiffSin2PhiStr3); 
  TLine *line3 = new TLine(3.14,0,3.14, binContent); 
  quad3HistAD->SetTitle("Quad 3: " + DecayTwo);
  quad3HistAD->Draw(); 
  quad3DSHistAD->Draw("same");
  line3->Draw("same");
  leg->Draw("same"); 
  pt1->Draw("same");  
  quad3CompareADCanvas->SaveAs("Quad3CompAD.png"); 
  leg->Clear("same"); 
  pt1->Clear("same");  
  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  quad4CompareADCanvas->cd(); 
  quad4HistAD->SetLineColor(kBlue); 
  quad4DSHistAD->SetLineColor(kRed); 
  quad4HistAD->GetYaxis()->SetLimits(0,1650);
  quad4HistAD->SetMinimum(0);
  maxBin = quad4HistAD->GetMaximumBin(); 
  binContent = 1.3*quad4HistAD->GetBinContent(maxBin); 
  quad4HistAD->SetMaximum(binContent); 
  TLine *line4 = new TLine(3.14,0,3.14, binContent); 
  leg->SetHeader("Legend"); 
  leg->AddEntry(quad4HistAD, "Dzero"); 
  leg->AddEntry(quad4DSHistAD, "Dzero-bar"); 
  pt1->AddText("PV D0:" + quad4D0ADAsymmetryStr + " #pm " + quad4D0ADAsymmetryErrStr);
  pt1->AddText("PV D0bar: " + quad4D0barADAsymmetryStr + " #pm " + quad4D0barADAsymmetryErrStr); 
  pt1->AddText("CPV: " + quad4AsymmetryADStr + " #pm" + quad4AsymmetryErrADStr);
  pt1->AddText("PV D0 Sin(2#phi): " + D0Sin2Phi_asymmetryStr4 + " #pm " + D0Sin2Phi_asymmetryErrStr4); 
    pt1->AddText("PV D0bar Sin(2#phi): " + D0barSin2Phi_asymmetryStr4 + " #pm " + D0barSin2Phi_asymmetryErrStr4); 
    pt1->AddText("CPV Sin(2#phi): " + Sin2PhiCP_asymmetryStr4 + " #pm " + Sin2PhiCP_asymmetryErrStr4); 
    pt1->AddText("D0-D0bar : " + DiffSin2PhiStr4); 
    quad4HistAD->SetTitle("Quad 4: " + DecayTwo);
  quad4HistAD->Draw(); 
  quad4DSHistAD->Draw("same");
  line4->Draw("same");
  leg->Draw("same"); 
  pt1->Draw("same");  
  quad4CompareADCanvas->SaveAs("Quad4CompAD.png"); 
  leg->Clear("same"); 
  pt1->Clear("same");  
  std::cout<<"thisistext"<<std::endl; 
 
}



