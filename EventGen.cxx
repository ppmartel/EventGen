// Compton Scattering and Pi0P/PiPN Photoproduction Event Generator
// Designed for use with the MAMI A2 Geant4 Simulation
//
// Author - P. Martel
// Version - 18 January 2013

#ifndef __CINT__

#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THStack.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TString.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TStopwatch.h"
#include "TLine.h"
#include "TArrow.h"
#include "TPaveText.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include <fstream>
using namespace std;

// Some physical constants.
#include "physics.h"
// Base Particle class
#include "BasePart.h"
// Base Generator class
#include "BaseGen.h"
// Compton Generator class
#include "CompGen.h"
// Pi0P Generator class
#include "Pi0PGen.h"
// PiPN Generator class
#include "PiPNGen.h"

int main()
{
  cout << "--------------------------------------------------" << endl;
  cout << "Compton Scattering and Pi0P/PiPN Photoproduction Event Generator" << endl;
  cout << "Designed for use with the MAMI A2 Geant4 Simulation" << endl << endl;
  cout << "Author - P. Martel" << endl;
  cout << "Version - 18 January 2013" << endl;
  cout << "--------------------------------------------------" << endl;

  // Set seed for random generator, otherwise the program will produce the
  // same set of random numbers each time it's run, which is sometimes useful
  // for testing the system. To do so, simply comment out the line.
  gRandom->SetSeed(0);

  TStopwatch swTimer;

  // Parameter Files
  
  gSystem->cd("par");
  TSystemDirectory *sdDir = new TSystemDirectory("files","./");
  TList *lList = sdDir->GetListOfFiles();
  lList->Sort();
  TIter itFile(lList);
  TSystemFile *sfFile;
  TString sFile;
  Int_t iLength;

  // Default beam and target polarization values

  Float_t fBeamT, fBeamL, fBeamP, fTargT, fTargL, fTargP;

  Float_t fBeamLo, fBeamHi;
  Int_t iBeamLo, iBeamHi;

  // Values to keep track of process and target selections

  const Int_t iNProc = 3;
  const Int_t iNType = 3;

  TString sProc[iNProc] = {"Comp","Pi0P","PiPN"};
  TString sType[iNType] = {"Normal","Incoherent","Coherent"};
  TString sWght;

  TString sTarg[6] = {"p","n","3He","4He","12C","16O"};
  Float_t fTargMass[6] = {kMP_MEV,kMN_MEV,kM_HE3_MEV,kM_HE4_MEV,kM_C12_MEV,kM_O16_MEV};

  TString sReco[6] = {"p","n","3He","4He","12C","16O"};
  Float_t fRecoMass[6] = {kMP_MEV,kMN_MEV,kM_HE3_MEV,kM_HE4_MEV,kM_C12_MEV,kM_O16_MEV};
  Int_t iRecoG3id[6] = {14, 13, 49, 47, 67, 69};

  Int_t iProcN, iTypeN, iWghtN;
  Int_t iProcS, iTypeS, iWghtS, iTargS, iRecoS;
  Int_t iBPolS, iTPolS;
  Int_t iTimeS;

  Int_t iNEvn = 1000000;
  Int_t iEvnN = 0;
  Int_t iEvnR = 0;

  // Choose reaction process

  cout << "Choose process:" << endl;
  for(iProcN=0; iProcN<iNProc; iProcN++) cout << iProcN+1 << ") " << sProc[iProcN] << endl;

  cin >> iProcS;
  cout << "--------------------------------------------------" << endl;

  if(iProcS<1 || iProcS>iNProc){
    cout << "Invalid process" << endl;
    return 0;
  }
  iProcS--;

  // Choose type of reaction

  cout << "Choose type:" << endl;
  for(iTypeN=0; iTypeN<iNType; iTypeN++) cout << iTypeN+1 << ") " << sType[iTypeN] << endl;

  cin >> iTypeS;
  cout << "--------------------------------------------------" << endl;

  if(iTypeS<1 || iTypeS>iNType){
    cout << "Invalid type" << endl;
    return 0;
  }
  iTypeS--;

  if(iTypeS!=2){
    cout << "Choose weighting:" << endl;
    cout << "1) Isot" << endl;
    iWghtN = 1;
    
    while((sfFile=(TSystemFile*)itFile())){
      sFile = sfFile->GetName();
      iLength = sFile.Length();
      if((iLength >= 10) && sFile.BeginsWith(sProc[iProcS]) && sFile.EndsWith(".txt")){
	sFile.Remove(0,1+sFile.First("_"));
	sFile.Remove(sFile.First("_"));
	cout << iWghtN+1 << ") " << sFile << endl;
	iWghtN++;
      }
    }
    
    cin >> iWghtS;
    cout << "--------------------------------------------------" << endl;
    
    if(iWghtS<1 || iWghtS>iWghtN){
      cout << "Invalid weighting" << endl;
    return 0;
    }
    iWghtS--;

    if(iWghtS == 0) sWght = "Isot";
    else{
      iWghtN = 1;
      itFile.Reset();
      while((sfFile=(TSystemFile*)itFile())){
	sFile = sfFile->GetName();
	iLength = sFile.Length();
	if((iLength >= 10) && sFile.BeginsWith(sProc[iProcS]) && sFile.EndsWith(".txt")){
	  sFile.Remove(0,1+sFile.First("_"));
	  sFile.Remove(sFile.First("_"));
	  if(iWghtS == iWghtN) sWght = sFile;
	  iWghtN++;
	}
      }
    }
  }

  else if(iTypeS==2){
    cout << "Choose weighting:" << endl;
    cout << "1) Isot" << endl;
    cout << "2) Parameterization" << endl;
    cin >> iWghtS;
    cout << "--------------------------------------------------" << endl;
    

    if(iProcS==2 && iWghtS==2){
      cout << "Coherent parameterization not provided for PiPN" << endl;
      return 0;
    }

    if(iWghtS<1 || iWghtS>2){
      cout << "Invalid weighting" << endl;
    return 0;
    }
    iWghtS--;

    if(iWghtS == 0) sWght = "Isot";
    else sWght = "parameterized";
  }

  // Choose target

  // Normal
  if(iTypeS==0){
    iTargS = 0;
    if(iProcS!=2) iRecoS = 0;
    else iRecoS = 1;
  }
  // Incoherent
  else if(iTypeS==1){
    iTargS = 4;
    if(iProcS!=2) iRecoS = 0;
    else iRecoS = 1;
  }
  // Coherent
  else{
    cout << "Choose target:" << endl;
    cout << "1) 3He" << endl;
    cout << "2) 4He" << endl;
    cout << "3) 12C" << endl;
    cout << "4) 16O" << endl;

    cin >> iTargS;
    cout << "--------------------------------------------------" << endl;

    if(iTargS<1 || iTargS>4){
      cout << "Invalid target" << endl;
      return 0;
    }
    iTargS++;
    iRecoS = iTargS;
  }

  cout << sTarg[iTargS] << " target chosen" << endl;
  cout << sReco[iRecoS] << " recoil chosen" << endl;
  cout << "--------------------------------------------------" << endl;

  // Set description of generator

  TString sName = (sType[iTypeS]+" "+sProc[iProcS]+" on "+sTarg[iTargS]+" using "+sWght+" weighting");

  // Set the beam and target polarizations for normal reactions if desired
  // For incoherent or coherent reactions, set target polarization to zero

  if(iTypeS==0 && iWghtS!=0){
    cout << "Choose beam polarization:" << endl;
    cout << "1) Unpolarized" << endl;
    cout << "2) Linear" << endl;
    cout << "3) Circular" << endl;

    cin >> iBPolS;

    if(iBPolS<1 || iBPolS>3){
      cout << "Invalid polarization" << endl;
      return 0;
    }

    if(iBPolS==1) fBeamT = fBeamL = fBeamP = 0;
    else if(iBPolS==2){
      cout << "Set polarization magnitude = ";
      cin >> fBeamT;
      if(fBeamT<0 || fBeamT>1){
	cout << "Invalid polarization" << endl;
	return 0;
      }
      cout << "Set polarization phi (deg) = ";
      cin >> fBeamP;
      if((TMath::Abs(fBeamP))>180){
	cout << "Invalid polarization" << endl;
	return 0;
      }
      fBeamL = 0;
    }
    else if(iBPolS==3){
      cout << "Set polarization magnitude = ";
      cin >> fBeamL;
      if((TMath::Abs(fBeamL))>1){
	cout << "Invalid polarization" << endl;
	return 0;
      }
      fBeamT = 0;
      fBeamP = 0;
    }
    cout << "--------------------------------------------------" << endl;

    cout << "Choose target polarization:" << endl;
    cout << "1) Unpolarized" << endl;
    cout << "2) Transverse" << endl;
    cout << "3) Longitudinal" << endl;

    cin >> iTPolS;

    if(iTPolS<1 || iTPolS>3){
      cout << "Invalid polarization" << endl;
      return 0;
    }

    if(iTPolS==1) fTargT = fTargL = fTargP = 0;
    else if(iTPolS==2){
      cout << "Set polarization magnitude = ";
      cin >> fTargT;
      if((TMath::Abs(fTargT))>1){
	cout << "Invalid polarization" << endl;
	return 0;
      }
      cout << "Set polarization phi (deg) = ";
      cin >> fTargP;
      if((TMath::Abs(fTargP))>180){
	cout << "Invalid polarization" << endl;
	return 0;
      }
      fTargL = 0;
    }
    else if(iTPolS==3){
      cout << "Set polarization magnitude = ";
      cin >> fTargL;
      if((TMath::Abs(fTargL))>1){
	cout << "Invalid polarization" << endl;
	return 0;
      }
      fTargT = 0;
      fTargP = 0;
    }
    cout << "--------------------------------------------------" << endl;
  }
  else{
    fBeamT = 0;
    fBeamL = 0;
    fBeamP = 0;
    fTargT = 0;
    fTargL = 0;
    fTargP = 0;
  }

  // Select beam energy range and determine if parameter files are available

  cout << "Minimum tagged photon energy (MeV) = ";
  cin >> fBeamLo;
  cout << "--------------------------------------------------" << endl;

  cout << "Maximum tagged photon energy (MeV) = ";
  cin >> fBeamHi;
  cout << "--------------------------------------------------" << endl;

  if(fBeamLo > fBeamHi){
    cout << "Invalid energy range" << endl;
    return 0;
  }

  TString sBase = sProc[iProcS]+"_"+sWght;
  Int_t iFileE;
  iBeamLo = 0;
  iBeamHi = 2000;

  if(iTypeS!=2 && iWghtS!=0){
    if(iProcS == 0) sBase += "_lab_";
    else sBase += "_cm_";

    itFile.Reset();
    while((sfFile=(TSystemFile*)itFile())){
      sFile = sfFile->GetName();
      iLength = sFile.Length();
      if((iLength >= 10) && sFile.BeginsWith(sBase)){
	sFile.ReplaceAll(sBase,"");
	iFileE = sFile.Atoi();
	if((iFileE > iBeamLo) && (iFileE <= fBeamLo)) iBeamLo = iFileE;
	if((iFileE < iBeamHi) && (iFileE >= fBeamHi)) iBeamHi = iFileE;
      }
    }
    
    if((iBeamLo == 0) || (iBeamHi == 2000)){
      cout << "Out of range of parameter files" << endl;
      return 0;
    }
  }
  else{
    iBeamLo = TMath::FloorNint(fBeamLo);
    iBeamHi = TMath::CeilNint(fBeamHi);
  }

  delete sdDir;
  delete lList;

  gSystem->cd("..");

  // Make Bremsstrahlung distribution for event selection

  //TF1 *fBeam = new TF1("fBeam", "1/x", fBeamLo, fBeamHi);
  TF1 *fBeam = new TF1("fBeam", "1", fBeamLo, fBeamHi);
  Bool_t bBeam = kTRUE;
  Double_t dBeam = iBeamLo;
  Bool_t bTime = kTRUE;

  // If selecting one energy value, turn off Brem distribution and create
  // cross section table just below and just above this value

  if(fBeamLo == fBeamHi){
    bBeam = kFALSE;
    bTime = kFALSE;
  }
  if(iBeamLo == iBeamHi){
    iBeamLo -= 5;
    iBeamHi += 5;
  }

  if(iWghtS == 0) bTime = kFALSE;

  Float_t fEnrE = 450.0;
  Float_t fEnrP = 160.0;

  Double_t dTagC = 0;
  Double_t dTagS = 1;
  Double_t dTagM = 1;

  // Variables for conversion constant

  Double_t dUBCM = 1e-30;
  Double_t dEffD = 1.0;
  Double_t dEffT = 1.0;
  Double_t dDenT = 4e23;
  Double_t dTime = 1;

  if(bTime){
    
    cout << "Choose run limit:" << endl;
    cout << "1) Run time" << endl;
    cout << "2) Number of events" << endl;
    
    cin >> iTimeS;
    
    if(iTimeS==2) bTime = kFALSE;
    else if(iTimeS<1 || iTimeS>2){
      cout << "Invalid run limit" << endl;
      return 0;
    }
  }

  if(bTime){

    // Variables for timed running

    cout << "Enter electron beam energy (MeV)" << endl;
    cin >> fEnrE;
    cout << "--------------------------------------------------" << endl;

    cout << "Enter lowest tagged photon energy (MeV)" << endl;
    cout << "(below which the tagger is turned off)" << endl;
    cin >> fEnrP;
    cout << "--------------------------------------------------" << endl;

    if(fEnrP == fBeamLo) dTagM = 1e6*60;
    else if(fEnrP < fBeamLo){
      dTagM = 1e6*60*fEnrP/fBeamLo;
      fEnrP = fBeamLo;
    }
    else{
      cout << "Tagger is turned off for part of the desired energy range" << endl;
      return 0;
    }

    // Tagger scalers are deadtime inhibited, so the calculation
    // automatically takes the livetime into account
    /*
    cout << "Enter average livetime" << endl;
    cin >> dEffD;
    cout << "--------------------------------------------------" << endl;
    if((dEffD>1) || (dEffD<0)){
      cout << "Invalid livetime" << endl;
      return 0;
    }
    */

    cout << "Enter average tagging efficiency" << endl;
    cin >> dEffT;
    cout << "--------------------------------------------------" << endl;
    if((dEffT>1) || (dEffT<0)){
      cout << "Invalid tagging efficiency" << endl;
      return 0;
    }

    cout << "Enter area target density (4.2e23 LH2, 9.1e22 FST)" << endl;
    cin >> dDenT;
    cout << "--------------------------------------------------" << endl;
    if(dDenT<0){
      cout << "Invalid target density" << endl;
      return 0;
    }

    cout << "Enter running time (min)" << endl;
    cout << "(though will still stop at a million events)" << endl;
    cin >> dTime;
    cout << "--------------------------------------------------" << endl;
    if(dTime<0){
      cout << "Invalid running time" << endl;
      return 0;
    }
    dTagM = dTagM*dTime;

  }
  else{

    cout << "Enter number of events:" << endl;
    cin >> iNEvn;
    cout << "--------------------------------------------------" << endl << endl;
    if(iNEvn<0){
      cout << "Invalid number of events" << endl;
      return 0;
    }

  }

  // Determine width of highest counting tagger channel, roughly

  TF1 *fChan = new TF1("fChan", "pol4", 0.0, 1.0);
  fChan->SetParameters(0.001415,0.003859,0.001417,-0.007596,0.004085);
  Float_t fEnrL = fEnrP;
  Float_t fEnrH = (fEnrP+(fEnrE*(fChan->Eval((fEnrE-fEnrP)/fEnrE))));

  // Determine conversion constant

  Double_t dConv = 1;
  if(bTime){
    dConv = dUBCM*dEffD*dEffT*dDenT;
    if((iProcS == 0) && (iTypeS != 2)) dConv = dConv*0.001;
  }

  // Set upper limit of scaling factor such that it would
  // require 100 throws to fill the max couting tagger channel

  Double_t dScalM = (dTagM/100.0);
  Double_t dScale = 1.0;

  // End of setup

  cout << sName << endl << endl;
  cout << "--------------------------------------------------" << endl << endl;

  TString sFile1 = "out/hist.root";
  TString sFile2 = "out/ntpl.root";
  TString sFile3 = "out/tree.root";

  // Compton generator

  if(iProcS == 0){

    // Create and initialize generator

    CompGen *pgen = new CompGen(sName,sTarg[iTargS],sBase,fTargMass[iTargS],fRecoMass[iRecoS],iRecoG3id[iRecoS],iBeamLo,iBeamHi);
    pgen->SetPol(fBeamT, fBeamL, fBeamP, fTargT, fTargL, fTargP);
    pgen->SetOut(kTRUE,kTRUE,kTRUE);
    pgen->Init();
    dScale = pgen->SetConv(dConv);
    if(bTime && dScale>dScalM) dScale = dScalM;
    dTagS = dScale;
    dScale = pgen->SetConv(dConv*dScale);

    // Events loop

    swTimer.Start();
    if(bTime){
      while(iEvnN<iNEvn && dTagC<dTagM){
	dBeam = fBeam->GetRandom();
	if(dBeam>=fEnrL && dBeam<fEnrH) dTagC += dTagS;
	if(pgen->NewEvent(dBeam)) iEvnN++;
	else iEvnR++;
      }
    }
    else{
      while(iEvnN<iNEvn){
	if(bBeam) dBeam = fBeam->GetRandom();
	if(pgen->NewEvent(dBeam)) iEvnN++;
	else iEvnR++;
      }
    }
    swTimer.Stop();

    // Read out results
  
    cout << iEvnN << " events accepted" << endl;
    cout << iEvnR << " events rejected" << endl;
    if(bTime){
      cout << dTagC << " events tagged (max of " << dTagM << ")" << endl;
      cout << "Equates to " << dTime*dTagC/dTagM << " min run time" << endl;
    }
    cout << "Computed in " << swTimer.RealTime() << " sec." << endl << endl;
    cout << "--------------------------------------------------" << endl;
    cout << endl;

    pgen->SaveHists(sFile1);
    pgen->SaveNtuple(sFile2);
    pgen->SaveTree(sFile3);
    delete pgen;
  }

  // Pi0P Photoproduction generator

  else if(iProcS == 1){

    // Create and initialize generator

    Pi0PGen *pgen = new Pi0PGen(sName,sTarg[iTargS],sBase,fTargMass[iTargS],fRecoMass[iRecoS],iRecoG3id[iRecoS],iBeamLo,iBeamHi);
    pgen->SetPol(fBeamT, fBeamL, fBeamP, fTargT, fTargL, fTargP);
    pgen->SetOut(kTRUE,kTRUE,kTRUE);
    pgen->Init();
    dScale = pgen->SetConv(dConv);
    if(bTime && dScale>dScalM) dScale = dScalM;
    dTagS = dScale;
    dScale = pgen->SetConv(dConv*dScale);

    // Events loop

    swTimer.Start();
    if(bTime){
      while(iEvnN<iNEvn && dTagC<dTagM){
	dBeam = fBeam->GetRandom();
	if(dBeam>=fEnrL && dBeam<fEnrH) dTagC += dTagS;
	if(pgen->NewEvent(dBeam)) iEvnN++;
	else iEvnR++;
      }
    }
    else{
      while(iEvnN<iNEvn){
	if(bBeam) dBeam = fBeam->GetRandom();
	if(pgen->NewEvent(dBeam)) iEvnN++;
	else iEvnR++;
      }
    }
    swTimer.Stop();

    // Read out results
  
    cout << iEvnN << " events accepted" << endl;
    cout << iEvnR << " events rejected" << endl;
    if(bTime){
      cout << dTagC << " events tagged (max of " << dTagM << ")" << endl;
      cout << "Equates to " << dTime*dTagC/dTagM << " min run time" << endl;
    }
    cout << "Computed in " << swTimer.RealTime() << " sec." << endl << endl;
    cout << "--------------------------------------------------" << endl;
    cout << endl;

    pgen->SaveHists(sFile1);
    pgen->SaveNtuple(sFile2);
    pgen->SaveTree(sFile3);
    delete pgen;
  }
  
  // PiPN Photoproduction generator

  else if(iProcS == 2){

    // Create and initialize generator

    PiPNGen *pgen = new PiPNGen(sName,sTarg[iTargS],sBase,fTargMass[iTargS],fRecoMass[iRecoS],iRecoG3id[iRecoS],iBeamLo,iBeamHi);
    pgen->SetPol(fBeamT, fBeamL, fBeamP, fTargT, fTargL, fTargP);
    pgen->SetOut(kTRUE,kTRUE,kTRUE);
    pgen->Init();
    dScale = pgen->SetConv(dConv);
    if(bTime && dScale>dScalM) dScale = dScalM;
    dTagS = dScale;
    dScale = pgen->SetConv(dConv*dScale);

    // Events loop

    swTimer.Start();
    if(bTime){
      while(iEvnN<iNEvn && dTagC<dTagM){
	dBeam = fBeam->GetRandom();
	if(dBeam>=fEnrL && dBeam<fEnrH) dTagC += dTagS;
	if(pgen->NewEvent(dBeam)) iEvnN++;
	else iEvnR++;
      }
    }
    else{
      while(iEvnN<iNEvn){
	if(bBeam) dBeam = fBeam->GetRandom();
	if(pgen->NewEvent(dBeam)) iEvnN++;
	else iEvnR++;
      }
    }
    swTimer.Stop();

    // Read out results
  
    cout << iEvnN << " events accepted" << endl;
    cout << iEvnR << " events rejected" << endl;
    if(bTime){
      cout << dTagC << " events tagged (max of " << dTagM << ")" << endl;
      cout << "Equates to " << dTime*dTagC/dTagM << " min run time" << endl;
    }
    cout << "Computed in " << swTimer.RealTime() << " sec." << endl << endl;
    cout << "--------------------------------------------------" << endl;
    cout << endl;

    pgen->SaveHists(sFile1);
    pgen->SaveNtuple(sFile2);
    pgen->SaveTree(sFile3);

    delete pgen;
  }

  cout << "--------------------------------------------------" << endl << endl;
  cout << ("Finished "+sName) << endl << endl;
  cout << "--------------------------------------------------" << endl << endl;

  return 0;
}

#endif
