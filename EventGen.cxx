// Compton Scattering and Pi0 Photoproduction Event Generator
// Designed for use with the MAMI A2 Geant4 Simulation
//
// Author - P. Martel
// Version - 21 June 2011

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
#include "ComptonGen.h"
// Pi0 Generator class
#include "Pi0PhotGen.h"

int main()
{
  cout << "--------------------------------------------------" << endl;
  cout << "Compton Scattering and Pi0 Photoproduction Event Generator" << endl;
  cout << "Designed for use with the MAMI A2 Geant4 Simulation" << endl << endl;
  cout << "Author - P. Martel" << endl;
  cout << "Version - 21 June 2011" << endl;
  cout << "--------------------------------------------------" << endl;

  // Set seed for random generator, otherwise the program will produce the
  // same set of random numbers each time it's run, which is sometimes useful
  // for testing the system. To do so, simply comment out the line.
  gRandom->SetSeed(0);

  TStopwatch swTimer;

  Int_t i=0, j=0;

  // Default beam and target polarization values
  Float_t fPolG = 0.7, fPolT = 0.8;

  Int_t iBeamLo, iBeamHi;

  // Values to keep track of process and target selections

  TString sProc[4] = {"Compton","Pi0 (MAID)","Pi0 (SAID)","Pi0"};
  TString sFnlo[3] = {"Compton_lab_","Pi0_maid_cm_","Pi0_said_cm_"};
  TString sFnhi[3] = {"Compton_lab_","Pi0_maid_cm_","Pi0_said_cm_"};
  TString sType[3] = {"Normal","Incoherent","Coherent"};
  TString sTarg[6] = {"p","n","3He","4He","12C","16O"};
  TString sName;

  Float_t fRecoMass[6] = {kMP_MEV,kMN_MEV,kM_HE3_MEV,kM_HE4_MEV,kM_C12_MEV,kM_O16_MEV};
  Int_t iRecoG3id[6] = {14, 13, 49, 47, 67, 69};
  
  Int_t iProcSel, iTypeSel, iTargSel, iRecoSel, iNEvn;
  char cCngPol[128];
  
  // Choose reaction process

  cout << "Choose process:" << endl;
  cout << "1) Compton" << endl;
  cout << "2) Pi0" << endl;

  cin >> iProcSel;
  cout << "--------------------" << endl;

  if(iProcSel<1 || iProcSel>2){
    cout << "Invalid process" << endl;
    return 0;
  }
  iProcSel--;

  // Choose type of reaction

  cout << "Choose type:" << endl;
  cout << "1) Normal" << endl;
  cout << "2) Incoherent" << endl;
  cout << "3) Coherent" << endl;

  cin >> iTypeSel;
  cout << "--------------------" << endl;

  if(iTypeSel<1 || iTypeSel>3){
    cout << "Invalid type" << endl;
    return 0;
  }
  iTypeSel--;

  if(iProcSel==1){
    if(iTypeSel==0 || iTypeSel==1){
      cout << "Choose Pi0 database:" << endl;
      cout << "1) MAID" << endl;
      cout << "2) SAID" << endl;

      cin >> iProcSel;
      cout << "--------------------" << endl;

      if(iProcSel<1 || iProcSel>2){
	cout << "Invalid process" << endl;
	return 0;
      }
    }
    else iProcSel = 3;
  }

  // For incoherent or coherent reactions, choose target

  if(iTypeSel==0){
    cout << "p target chosen" << endl;
    cout << "--------------------" << endl;
    iTargSel = 0;
  }
  else if(iTypeSel==1){
    cout << "12C target chosen" << endl;
    cout << "--------------------" << endl;
    iTargSel = 4;
  }
  else{
    cout << "Choose target:" << endl;
    cout << "1) 3He" << endl;
    cout << "2) 4He" << endl;
    cout << "3) 12C" << endl;
    cout << "4) 16O" << endl;

    cin >> iTargSel;
    cout << "--------------------" << endl;

    if(iTargSel<1 || iTargSel>4){
      cout << "Invalid target" << endl;
      return 0;
    }
    iTargSel++;
  }

  // For incoherent reactions, choose recoil nucleon
  
  if(iTypeSel==1){
    /*
    cout << "Choose recoil:" << endl;
    cout << "1) p" << endl;
    cout << "2) n" << endl;

    cin >> iRecoSel;
    cout << "--------------------" << endl;

    if(iRecoSel<1 || iRecoSel>2){
      cout << "Invalid recoil" << endl;
      return 0;
    }
    iRecoSel--;    
    */

    // For now, recoil from incoherent reaction is always a proton
    // eventually a neutron recoil may be added in

    iRecoSel = 0;
  }
  else iRecoSel = iTargSel;

  // Set the beam and target polarizations for normal reactions if desired
  // For incoherent or coherent reactions, set target polarization to zero

  if(iTypeSel==0){
    cout << "Polarization settings:" << endl;
    cout << "Beam pol = " << fPolG << endl;
    cout << "Target pol = " << fPolT << endl;
    cout << "Change these?" << endl;

    cin.ignore();
    cin.getline(cCngPol,128);

    if(!strncmp(cCngPol,"y",4) || !strncmp(cCngPol,"yes",4) || !strncmp(cCngPol,"Y",4) || !strncmp(cCngPol,"Yes",4) || !strncmp(cCngPol,"YES",4)){
      cout << "Set Beam pol = ";
      cin >> fPolG;
      if((TMath::Abs(fPolG))>1){
	cout << "Invalid polarization" << endl;
	return 0;
      }
      cout << "Set Target pol = ";
      cin >> fPolT;
      if((TMath::Abs(fPolT))>1){
	cout << "Invalid polarization" << endl;
	return 0;
      }
    }
    cout << "--------------------------------------------------" << endl;
  }
  else fPolT = 0;

  // Select beam energy range and determine if parameter files are available

  cout << "Enter minimum tagged photon energy (integers of MeV)" << endl;
  cin >> iBeamLo;
  cout << "--------------------------------------------------" << endl;

  cout << "Enter maximum tagged photon energy (integers of MeV)" << endl;
  cin >> iBeamHi;
  cout << "--------------------------------------------------" << endl;

  if(iTypeSel!=2){
    sFnlo[iProcSel] += iBeamLo;
    sFnhi[iProcSel] += iBeamHi;

    TSystemDirectory *sdParDir = new TSystemDirectory("parameters","./par");
    TList *lParList = sdParDir->GetListOfFiles();
    lParList->Sort();
    
    if(!(lParList->Contains(sFnlo[iProcSel]) &&
	 lParList->Contains(sFnhi[iProcSel]) &&
	 iBeamHi>iBeamLo)){
      cout << "Invalid energy range" << endl;
      return 0;
    }
    
    delete sdParDir;
    delete lParList;
  }

  // Make Bremsstrahlung distribution for event selection
  TF1 *fBeam = new TF1("fBeam", "1/x", iBeamLo, iBeamHi);

  // Other initial setup parts

  cout << "Enter number of events:" << endl;
  cin >> iNEvn;
  cout << "--------------------------------------------------" << endl << endl;

  sName = (sType[iTypeSel]+" "+sProc[iProcSel]+" on "+sTarg[iTargSel]);

  cout << (sName+" selected") << endl << endl;
  cout << "--------------------------------------------------" << endl << endl;

  TString sFile1 = "out/hist.root";
  TString sFile2 = "out/ntpl.root";
  
  // Compton generator

  if(sProc[iProcSel].Contains("Compton")){

    // Create and initialize generator

    ComptonGen *cgen = new ComptonGen(sName,sTarg[iTargSel],fRecoMass[iRecoSel],iRecoG3id[iRecoSel],iBeamLo,iBeamHi);
    cgen->Init(fPolG, fPolT);

    // Events loop

    swTimer.Start();
    while(i<iNEvn){
      if(cgen->NewEvent(fBeam->GetRandom())) i++;
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
    delete cgen;
  }
  
  // Pi0 Photoproduction generator

  else if(sProc[iProcSel].Contains("Pi0")){

    // Create and initialize generator

    Pi0PhotGen *pgen = new Pi0PhotGen(sName,sTarg[iTargSel],fRecoMass[iRecoSel],iRecoG3id[iRecoSel],iBeamLo,iBeamHi);
    pgen->Init(fPolG, fPolT);

    // Events loop

    swTimer.Start();
    while(i<iNEvn){
      if(pgen->NewEvent(fBeam->GetRandom())) i++;
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
    delete pgen;
  }

  return 0;
}

#endif
