// Compton Scattering and Pi0P/PiPN Photoproduction Event Generator
// Designed for use with the MAMI A2 Geant4 Simulation
//
// Author - P. Martel
// Version - 23 October 2012

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
  cout << "Version - 17 October 2012" << endl;
  cout << "--------------------------------------------------" << endl;

  // Set seed for random generator, otherwise the program will produce the
  // same set of random numbers each time it's run, which is sometimes useful
  // for testing the system. To do so, simply comment out the line.
  gRandom->SetSeed(0);

  TStopwatch swTimer;

  Int_t i=0, j=0;

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
  Int_t iProcS, iTypeS, iWghtS, iTargS, iRecoS, iNEvn;
  Int_t iBPolS, iTPolS;
  
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

  if(iProcS==2 && iTypeS==2){
    cout << "Coherent process not allowed for PiPN" << endl;
    return 0;
  }

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

  TString sName = (sType[iTypeS]+" "+sProc[iProcS]+" on "+sTarg[iTargS]);
  if(iTypeS!=2) sName += (" using "+sWght+" weighting");

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

  cout << "Enter minimum tagged photon energy (integers of MeV)" << endl;
  cin >> iBeamLo;
  cout << "--------------------------------------------------" << endl;

  cout << "Enter maximum tagged photon energy (integers of MeV)" << endl;
  cin >> iBeamHi;
  cout << "--------------------------------------------------" << endl;

  if(iBeamLo > iBeamHi){
    cout << "Invalid energy range" << endl;
    return 0;
  }

  TString sBase = sProc[iProcS]+"_"+sWght;

  if(iTypeS!=2 && iWghtS!=0){
    if(iProcS == 0) sBase += "_lab_";
    else sBase += "_cm_";
    
    TString sFnLo = sBase;
    sFnLo += iBeamLo;
    TString sFnHi = sBase;
    sFnHi += iBeamHi;

    if(!(lList->Contains(sFnLo) && lList->Contains(sFnHi))){
      cout << "Invalid energy range" << endl;
      return 0;
    }
  }

  delete sdDir;
  delete lList;

  gSystem->cd("..");

  // Make Bremsstrahlung distribution for event selection

  TF1 *fBeam = new TF1("fBeam", "1/x", iBeamLo, iBeamHi);
  Bool_t bBeam = kTRUE;
  Double_t dBeam = iBeamLo;

  // If selecting one energy value, turn off Brem distribution and create
  // cross section table just below and just above this value

  if(iBeamLo == iBeamHi){
    bBeam = kFALSE;
    dBeam = iBeamLo;
    iBeamLo -= 5;
    iBeamHi += 5;
  }

  // Other initial setup parts

  cout << "Enter number of events:" << endl;
  cin >> iNEvn;
  cout << "--------------------------------------------------" << endl << endl;

  cout << sName << endl << endl;
  cout << "--------------------------------------------------" << endl << endl;

  TString sFile1 = "out/hist.root";
  TString sFile2 = "out/ntpl.root";
  TString sFile3 = "out/tree.root";

  // Compton generator

  if(iProcS == 0){

    // Create and initialize generator

    CompGen *cgen = new CompGen(sName,sTarg[iTargS],sBase,fTargMass[iTargS],fRecoMass[iRecoS],iRecoG3id[iRecoS],iBeamLo,iBeamHi);
    cgen->SetPol(fBeamT, fBeamL, fBeamP, fTargT, fTargL, fTargP);
    cgen->Init();

    // Events loop

    swTimer.Start();
    while(i<iNEvn){
      if(bBeam) dBeam = fBeam->GetRandom();
      if(cgen->NewEvent(dBeam)) i++;
      else{
	j++;
	continue;
      }
    }
    swTimer.Stop();

    // Read out results
  
    cout << i << " events accepted" << endl;
    cout << j << " events rejected" << endl;
    cout << "Computed in " << swTimer.RealTime() << " sec." << endl << endl;
    cout << "--------------------------------------------------" << endl;
    cout << endl;

    cgen->SaveHists(sFile1);
    cgen->SaveNtuple(sFile2);
    cgen->SaveTree(sFile3);
    delete cgen;
  }

  // Pi0P Photoproduction generator

  else if(iProcS == 1){

    // Create and initialize generator

    Pi0PGen *pgen = new Pi0PGen(sName,sTarg[iTargS],sBase,fTargMass[iTargS],fRecoMass[iRecoS],iRecoG3id[iRecoS],iBeamLo,iBeamHi);
    pgen->SetPol(fBeamT, fBeamL, fBeamP, fTargT, fTargL, fTargP);
    pgen->Init();

    // Events loop

    swTimer.Start();
    while(i<iNEvn){
      if(bBeam) dBeam = fBeam->GetRandom();
      if(pgen->NewEvent(dBeam)) i++;
      else{
	j++;
	continue;
      }
    }
    swTimer.Stop();

    // Read out results
  
    cout << i << " events accepted" << endl;
    cout << j << " events rejected" << endl;
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
    pgen->Init();

    // Events loop

    swTimer.Start();
    while(i<iNEvn){
      if(bBeam) dBeam = fBeam->GetRandom();
      if(pgen->NewEvent(dBeam)) i++;
      else{
	j++;
	continue;
      }
    }
    swTimer.Stop();

    // Read out results
  
    cout << i << " events accepted" << endl;
    cout << j << " events rejected" << endl;
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
