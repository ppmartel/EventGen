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
#include "TRandom.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TStopwatch.h"
#include <fstream>
using namespace std;

int TestGen();
//void SpecModel(Float_t*);
void SpecModel(TLorentzVector*);
TString GenNames(Int_t, Int_t*);
Float_t Sqr(Float_t);

int main()
{
  return TestGen();
}

#endif

// Some physical constants.
#include "physics.h"
// Some useful functions.
#include "functions.h"

int TestGen() 
{
  Int_t i, j, k, m, AngN, NAng, PhiN, NPhi, EnrN, NEnr, EvnN, NEvn, proc;
  Int_t qfDet = 0, kfDet = 0, p1Det = 0, p2Det=0;
  Int_t NComp[3][3] = {{0}};
  Int_t NPi0[3][3][3] = {{{0}}};
  Bool_t RandSel = kFALSE;
  TString gnames;
  Int_t npart, ptag[10];

  TStopwatch timer;

  cout << "Choose process:" << endl;
  cout << "1) p(g,g)p" << endl;
  cout << "2) p(g,pi0)p" << endl;
  cout << "3) A(g,g)pX" << endl;
  cout << "4) A(g,pi0)pX" << endl;
  cout << "5) Random" << endl;
  cout << "--------------------" << endl;

  cin >> proc;

  if(proc==1) cout << "p(g,g)p selected" << endl << endl;
  else if(proc==2) cout << "p(g,pi0)p selected" << endl << endl;
  else if(proc==3) cout << "A(g,g)pX selected" << endl << endl;
  else if(proc==4) cout << "A(g,pi0)pX selected" << endl << endl;
  else if(proc==5){
    cout << "Random selection" << endl << endl;
    RandSel = kTRUE;
  }
  else{
    cout << "Invalid selection" << endl;
    return 0;
  }

  if(proc==1 || proc==3){
    npart = 2;
    ptag[0] = 14;
    ptag[1] = 1;
  }
  if(proc==2 || proc==4){
    npart = 4;
    ptag[0] = 14;
    ptag[1] = 7;
    ptag[2] = 1;
    ptag[3] = 1;
  }
  gnames = GenNames(npart, ptag);

  AngN = 0;
  NAng = 181;
  PhiN = 0;
  NPhi = 181;
  EnrN = 0;
  NEnr = 3;
  EvnN = 0;
  NEvn = 100000;

  Float_t Co_Angle[181] = {0.0};
  Float_t Co_SolAng = 0;
  Float_t Co_SolAng_tot = 0;
  Float_t Co_U[21][181] = {{0.0}};
  Float_t Co_P[3][181] = {{0.0}};
  Float_t Co_M[3][181] = {{0.0}};
  Float_t Co_A[3][181] = {{0.0}};
  Float_t Co[21][181][181] = {{{0.0}}};
  Float_t Co_m[2][181][181] = {{{0.0}}};
  Float_t Co_b[2][181][181] = {{{0.0}}};
  Float_t Co_max[21] = {0.0};
  Float_t Co_max_m[2] = {0.0};
  Float_t Co_max_b[2] = {0.0};
  Float_t Co_U_m[2][181] = {{0.0}};
  Float_t Co_U_b[2][181] = {{0.0}};
  Float_t Co_U_max[21] = {0.0};
  Float_t Co_U_max_m[2] = {0.0};
  Float_t Co_U_max_b[2] = {0.0};
  Float_t Co_tot[21] = {0.0};
  Float_t Co_tot_m[2] = {0.0};
  Float_t Co_tot_b[2] = {0.0};

  Float_t Pi_Angle[181] = {0.0};
  Float_t Pi_SolAng = 0;
  Float_t Pi_SolAng_tot = 0;
  Float_t Pi_S[3][181] = {{0.0}};
  Float_t Pi_T[3][181] = {{0.0}};
  Float_t Pi_E[3][181] = {{0.0}};
  Float_t Pi_F[3][181] = {{0.0}};
  Float_t Pi_U[21][181] = {{0.0}};
  Float_t Pi[21][181][181] = {{{0.0}}};
  Float_t Pi_m[2][181][181] = {{{0.0}}};
  Float_t Pi_b[2][181][181] = {{{0.0}}};
  Float_t Pi_max[21] = {0.0};
  Float_t Pi_max_m[2] = {0.0};
  Float_t Pi_max_b[2] = {0.0};
  Float_t Pi_U_m[2][181] = {{0.0}};
  Float_t Pi_U_b[2][181] = {{0.0}};
  Float_t Pi_U_max[21] = {0.0};
  Float_t Pi_U_max_m[2] = {0.0};
  Float_t Pi_U_max_b[2] = {0.0};
  Float_t Pi_tot[21] = {0.0};
  Float_t Pi_tot_m[2] = {0.0};
  Float_t Pi_tot_b[2] = {0.0};

  Float_t Pol_g = 0.7;
  Float_t Pol_t = 0.8;

  char Co_Text[3][100] = {"par/Comp_lab_","par/Comp_lab_","par/Comp_lab_"};
  char Pi_Text[3][100] = {"par/Pi0_cm_","par/Pi0_cm_","par/Pi0_cm_"};

  Int_t BeamE = 280;
  Int_t Energy[3] = {BeamE-10, BeamE, BeamE+10};
  char EnergyText[3][4];

  const char *Co_File[3];
  const char *Pi_File[3];

  for(EnrN=0; EnrN<NEnr; EnrN++){
    sprintf(EnergyText[EnrN],"%d",Energy[EnrN]);

    strcat(Co_Text[EnrN], EnergyText[EnrN]);
    strcat(Pi_Text[EnrN], EnergyText[EnrN]);
    
    Co_File[EnrN] = Co_Text[EnrN];
    Pi_File[EnrN] = Pi_Text[EnrN];

    AngN = 0;
    
    ifstream fin0(Co_File[EnrN]);
    while(!fin0.eof()){
      fin0>>Co_Angle[180-AngN];
      fin0>>Co_U[EnrN*10][180-AngN];
      fin0>>Co_P[EnrN][180-AngN];
      fin0>>Co_M[EnrN][180-AngN];
      
      AngN++;
    }
    fin0.close();
    
    AngN = 0;
    
    ifstream fin1(Pi_File[EnrN]);
    while(!fin1.eof()){
      fin1>>Pi_Angle[AngN];
      fin1>>Pi_S[EnrN][AngN];
      fin1>>Pi_T[EnrN][AngN];
      fin1>>Pi_E[EnrN][AngN];
      fin1>>Pi_F[EnrN][AngN];
      
      AngN++;
    }
    fin1.close();
  }

  Float_t qi_m, ki_m, qf_m, kf_m, qi_e, ki_e, qf_e, kf_e;

  Float_t qi_p, ki_p, qf_p, kf_p, vtx_x, vtx_y, vtx_z;

  Float_t qth, qph, qth_cm, qph_cm, pth_pi, pph_pi;

  Float_t weight, cross, cross_max, proc_weight[4][21], proc_tot[21];
  Float_t targ_weight[2] = {(10./42.),(32./42.)};

  char proc_text[4][10] = {"Comp","Pi0","IncohComp","IncohPi0"};
  
  TLorentzVector qi, ki, qf, kf, qi_cm, ki_cm, qf_cm, kf_cm, ptot, ptot_cm;

  TLorentzVector p1, p2, p1_pi, p2_pi;

  TVector3 vtx, lab_to_cm, cm_to_lab, lab_to_pi, pi_to_lab;

  qi_m = 0;
  ki_m = kMP_MEV;
  qf_m = 0;
  kf_m = kMP_MEV;
  ki_e = ki_m;

  TNtuple *h1 = new TNtuple("h1", "TMCUserGenerator", gnames);

  // Array for filling ntuple
  Float_t var[100];

  TH1F *hqi[5][6], *hki[5][6], *hqf[5][6], *hkf[5][6], *hp1[5][6], *hp2[5][6];

  THStack *sqi[6], *ski[6], *sqf[6], *skf[6], *sp1[6], *sp2[6];

  char h_name[10], h_title[25];

  char qi_part[10] = "Photon_i";
  char ki_part[10] = "Proton_i";
  char qf_part[5][10] = {"Scatter_f","Photon_f","Pi0_f","Photon_f","Pi0_f"};
  char kf_part[10] = "Proton_f";
  char p1_part[10] = "Decay_1";
  char p2_part[10] = "Decay_2";
  char title_end[6][15] = {"E (MeV)","xMom (MeV)","yMom (MeV)","zMom (MeV)","Mass (MeV)","Theta (deg)"};

  Float_t hist_lim[6][2] = {{0,300},{-500,500},{-500,500},{-500,500},{0,1000},{0,180}};

  for(i=0; i<5; i++){
    for(j=0; j<6; j++){
      sprintf(h_name,"hqi_%d_%d",i,j);
      sprintf(h_title,"%s %s",qi_part,title_end[j]);
      hqi[i][j] = new TH1F(h_name,h_title,1001,hist_lim[j][0],hist_lim[j][1]);

      sprintf(h_name,"hki_%d_%d",i,j);
      sprintf(h_title,"%s %s",ki_part,title_end[j]);
      hki[i][j] = new TH1F(h_name,h_title,1001,hist_lim[j][0],hist_lim[j][1]);

      sprintf(h_name,"hqf_%d_%d",i,j);
      sprintf(h_title,"%s %s",qf_part[i],title_end[j]);
      hqf[i][j] = new TH1F(h_name,h_title,1001,hist_lim[j][0],hist_lim[j][1]);

      sprintf(h_name,"hkf_%d_%d",i,j);
      sprintf(h_title,"%s %s",kf_part,title_end[j]);
      hkf[i][j] = new TH1F(h_name,h_title,1001,hist_lim[j][0],hist_lim[j][1]);

      sprintf(h_name,"hp1_%d_%d",i,j);
      sprintf(h_title,"%s %s",p1_part,title_end[j]);
      hp1[i][j] = new TH1F(h_name,h_title,1001,hist_lim[j][0],hist_lim[j][1]);

      sprintf(h_name,"hp2_%d_%d",i,j);
      sprintf(h_title,"%s %s",p2_part,title_end[j]);
      hp2[i][j] = new TH1F(h_name,h_title,1001,hist_lim[j][0],hist_lim[j][1]);
    }
  }

  for(j=0; j<6; j++){
    sprintf(h_name,"sqi_%d",j);
    sprintf(h_title,"%s %s",qi_part,title_end[j]);
    sqi[j] = new THStack(h_name,h_title);
    
    sprintf(h_name,"ski_%d",j);
    sprintf(h_title,"%s %s",ki_part,title_end[j]);
    ski[j] = new THStack(h_name,h_title);
    
    sprintf(h_name,"sqf_%d",j);
    sprintf(h_title,"%s %s",qf_part[0],title_end[j]);
    sqf[j] = new THStack(h_name,h_title);
    
    sprintf(h_name,"skf_%d",j);
    sprintf(h_title,"%s %s",kf_part,title_end[j]);
    skf[j] = new THStack(h_name,h_title);
    
    sprintf(h_name,"sp1_%d",j);
    sprintf(h_title,"%s %s",p1_part,title_end[j]);
    sp1[j] = new THStack(h_name,h_title);
    
    sprintf(h_name,"sp2_%d",j);
    sprintf(h_title,"%s %s",p2_part,title_end[j]);
    sp2[j] = new THStack(h_name,h_title);
  }

  TCanvas *c1 = new TCanvas("c1", "Co_Sigma", 800, 1000);
  TCanvas *c2 = new TCanvas("c2", "Pi_Sigma", 800, 1000);
  TCanvas *c3 = new TCanvas("c3", "Incoh_Co_Sigma", 800, 1000);
  TCanvas *c4 = new TCanvas("c4", "Incoh_Pi_Sigma", 800, 1000);

  TH2F *hco[3];
  TH2F *hpi[3];

  TH1F *hco_U[3];
  TH1F *hpi_U[3];

  TH1F *hproc_i[21];
  TH1F *hproc_f[21];

  TH2F *hangle = new TH2F("hangle","Proton vs Scatter Angle",181,0,181,181,0,181);
  TH2F *hdeang = new TH2F("hdeang","Photon2 vs Photon1 Angle",181,0,181,181,0,181);
  TH2F *hmissm = new TH2F("hmissm","Missing Mass vs Proton K",301,0,301,301,0,301);

  char Co_Name[3][10] = {"hco_0", "hco_1", "hco_2"};
  char Pi_Name[3][10] = {"hpi_0", "hpi_1", "hpi_2"};

  char Co_U_Name[3][10] = {"hco_U_0", "hco_U_1", "hco_U_2"};
  char Pi_U_Name[3][10] = {"hpi_U_0", "hpi_U_1", "hpi_U_2"};

  char Proc_i_Name[21][15];
  char Proc_f_Name[21][15];

  for(i=0; i<21; i++){
    sprintf(Proc_i_Name[i],"hproc_i_%d",i);
    sprintf(Proc_f_Name[i],"hproc_f_%d",i);
  }

  for(EnrN=0; EnrN<3; EnrN++){
    qi_e = Energy[EnrN];
    qi.SetPxPyPzE(0,0,qi_e,qi_e);
    ki.SetPxPyPzE(0,0,0,ki_e);
    ptot = qi + ki;

    hco[EnrN] = new TH2F(Co_Name[EnrN],"Comp Cross Sections",181,0,181,181,0,181);
    hpi[EnrN] = new TH2F(Pi_Name[EnrN],"Pi0 Cross Sections",181,0,181,181,0,181);

    hco_U[EnrN] = new TH1F(Co_U_Name[EnrN],"Incoh Comp Cross Sections",181,0,181);
    hpi_U[EnrN] = new TH1F(Pi_U_Name[EnrN],"Incoh Pi0 Cross Sections",181,0,181);

    for(AngN=0; AngN<NAng; AngN++){
      Co_A[EnrN][AngN]=(Pol_g*Pol_t*((Co_P[EnrN][AngN]-Co_M[EnrN][AngN])/(Co_P[EnrN][AngN]+Co_M[EnrN][AngN])));

      for(PhiN=0; PhiN<NPhi; PhiN++){
	Co[EnrN*10][AngN][PhiN] = (Co_U[EnrN*10][AngN]*(1+Co_A[EnrN][AngN]*cos(PhiN*kD2R)));
	hco[EnrN]->Fill(AngN, PhiN, Co[EnrN*10][AngN][PhiN]);

	if(Co[EnrN*10][AngN][PhiN]>Co_max[EnrN*10]){
	  Co_max[EnrN*10] = Co[EnrN*10][AngN][PhiN];
	}
	
	Pi[EnrN*10][AngN][PhiN] = (Pi_S[EnrN][AngN]*(1+Pol_t*sin(-PhiN*kD2R)*Pi_T[EnrN][AngN]+Pol_g*(ptot.Beta()*ptot.Gamma()*Pi_E[EnrN][AngN]+Pol_t*cos(-PhiN*kD2R)*Pi_F[EnrN][AngN])));
	hpi[EnrN]->Fill(AngN, PhiN, Pi[EnrN*10][AngN][PhiN]);

	if(Pi[EnrN*10][AngN][PhiN]>Pi_max[EnrN*10]){
	  Pi_max[EnrN*10] = Pi[EnrN*10][AngN][PhiN];
	}	
      }

      hco_U[EnrN]->Fill(AngN, Co_U[EnrN*10][AngN]);
      if(Co_U[EnrN*10][AngN]>Co_U_max[EnrN*10]){
	Co_U_max[EnrN*10] = Co_U[EnrN*10][AngN];
      }
      
      Pi_U[EnrN*10][AngN]=(Pi_S[EnrN][AngN]*(1+Pol_g*(ptot.Beta()*ptot.Gamma()*Pi_E[EnrN][AngN])));
      hpi_U[EnrN]->Fill(AngN, Pi_U[EnrN*10][AngN]);
      if(Pi_U[EnrN*10][AngN]>Pi_U_max[EnrN*10]){
	Pi_U_max[EnrN*10] = Pi_U[EnrN*10][AngN];
      }
    }

    Co_SolAng_tot = 0;
    Pi_SolAng_tot = 0;

    for(AngN=0; AngN<NAng-1; AngN++){
      Co_SolAng = ((cos(Co_Angle[AngN]*kD2R)-cos(Co_Angle[AngN+1]*kD2R))*kD2R);
      Pi_SolAng = ((cos(Pi_Angle[AngN]*kD2R)-cos(Pi_Angle[AngN+1]*kD2R))*kD2R);

      for(PhiN=0; PhiN<NPhi-1; PhiN++){
	Co_tot[EnrN*10] += (2*Co[EnrN*10][AngN][PhiN]*Co_SolAng);
	Co_SolAng_tot += (2*Co_SolAng);

	Pi_tot[EnrN*10] += (2*Pi[EnrN*10][AngN][PhiN]*Pi_SolAng);
	Pi_SolAng_tot += (2*Pi_SolAng);
      }
    }
  }

  Co_max_m[0] = ((Co_max[10]-Co_max[0])/(Energy[1]-Energy[0]));
  Co_max_b[0] = (Co_max[10]-(Co_max_m[0]*Energy[1]));
  Co_max_m[1] = ((Co_max[20]-Co_max[10])/(Energy[2]-Energy[1]));
  Co_max_b[1] = (Co_max[10]-(Co_max_m[1]*Energy[1]));

  Co_U_max_m[0] = ((Co_U_max[10]-Co_U_max[0])/(Energy[1]-Energy[0]));
  Co_U_max_b[0] = (Co_U_max[10]-(Co_U_max_m[0]*Energy[1]));
  Co_U_max_m[1] = ((Co_U_max[20]-Co_U_max[10])/(Energy[2]-Energy[1]));
  Co_U_max_b[1] = (Co_U_max[10]-(Co_U_max_m[1]*Energy[1]));

  Co_tot_m[0] = ((Co_tot[10]-Co_tot[0])/(Energy[1]-Energy[0]));
  Co_tot_b[0] = (Co_tot[10]-(Co_tot_m[0]*Energy[1]));
  Co_tot_m[1] = ((Co_tot[20]-Co_tot[10])/(Energy[2]-Energy[1]));
  Co_tot_b[1] = (Co_tot[10]-(Co_tot_m[1]*Energy[1]));

  Pi_max_m[0] = ((Pi_max[10]-Pi_max[0])/(Energy[1]-Energy[0]));
  Pi_max_b[0] = (Pi_max[10]-(Pi_max_m[0]*Energy[1]));
  Pi_max_m[1] = ((Pi_max[20]-Pi_max[10])/(Energy[2]-Energy[1]));
  Pi_max_b[1] = (Pi_max[10]-(Pi_max_m[1]*Energy[1]));

  Pi_U_max_m[0] = ((Pi_U_max[10]-Pi_U_max[0])/(Energy[1]-Energy[0]));
  Pi_U_max_b[0] = (Pi_U_max[10]-(Pi_U_max_m[0]*Energy[1]));
  Pi_U_max_m[1] = ((Pi_U_max[20]-Pi_U_max[10])/(Energy[2]-Energy[1]));
  Pi_U_max_b[1] = (Pi_U_max[10]-(Pi_U_max_m[1]*Energy[1]));

  Pi_tot_m[0] = ((Pi_tot[10]-Pi_tot[0])/(Energy[1]-Energy[0]));
  Pi_tot_b[0] = (Pi_tot[10]-(Pi_tot_m[0]*Energy[1]));
  Pi_tot_m[1] = ((Pi_tot[20]-Pi_tot[10])/(Energy[2]-Energy[1]));
  Pi_tot_b[1] = (Pi_tot[10]-(Pi_tot_m[1]*Energy[1]));

  for(AngN=0; AngN<NAng; AngN++){
    for(PhiN=0; PhiN<NPhi; PhiN++){
      Co_m[0][AngN][PhiN] = ((Co[10][AngN][PhiN]-Co[0][AngN][PhiN])/(Energy[1]-Energy[0]));
      Co_b[0][AngN][PhiN] = (Co[10][AngN][PhiN]-(Co_m[0][AngN][PhiN]*Energy[1]));
      Co_m[1][AngN][PhiN] = ((Co[20][AngN][PhiN]-Co[10][AngN][PhiN])/(Energy[2]-Energy[1]));
      Co_b[1][AngN][PhiN] = (Co[10][AngN][PhiN]-(Co_m[1][AngN][PhiN]*Energy[1]));

      Pi_m[0][AngN][PhiN] = ((Pi[10][AngN][PhiN]-Pi[0][AngN][PhiN])/(Energy[1]-Energy[0]));
      Pi_b[0][AngN][PhiN] = (Pi[10][AngN][PhiN]-(Pi_m[0][AngN][PhiN]*Energy[1]));
      Pi_m[1][AngN][PhiN] = ((Pi[20][AngN][PhiN]-Pi[10][AngN][PhiN])/(Energy[2]-Energy[1]));
      Pi_b[1][AngN][PhiN] = (Pi[10][AngN][PhiN]-(Pi_m[1][AngN][PhiN]*Energy[1]));
    }

    Co_U_m[0][AngN] = ((Co_U[10][AngN]-Co_U[0][AngN])/(Energy[1]-Energy[0]));
    Co_U_b[0][AngN] = (Co_U[10][AngN]-(Co_U_m[0][AngN]*Energy[1]));
    Co_U_m[1][AngN] = ((Co_U[20][AngN]-Co_U[10][AngN])/(Energy[2]-Energy[1]));
    Co_U_b[1][AngN] = (Co_U[10][AngN]-(Co_U_m[1][AngN]*Energy[1]));

    Pi_U_m[0][AngN] = ((Pi_U[10][AngN]-Pi_U[0][AngN])/(Energy[1]-Energy[0]));
    Pi_U_b[0][AngN] = (Pi_U[10][AngN]-(Pi_U_m[0][AngN]*Energy[1]));
    Pi_U_m[1][AngN] = ((Pi_U[20][AngN]-Pi_U[10][AngN])/(Energy[2]-Energy[1]));
    Pi_U_b[1][AngN] = (Pi_U[10][AngN]-(Pi_U_m[1][AngN]*Energy[1]));

  }
  
  for(i=Energy[0]; i<=Energy[2]; i+=1){
    EnrN = (i-Energy[0]);
    if(i!=Energy[0] && i!=Energy[1] && i!=Energy[2]){
      for(AngN=0; AngN<NAng; AngN++){
	for(PhiN=0; PhiN<NPhi; PhiN++){
	  Co[EnrN][AngN][PhiN] = ((i*Co_m[0][AngN][PhiN])+Co_b[0][AngN][PhiN]);
	  Pi[EnrN][AngN][PhiN] = ((i*Pi_m[0][AngN][PhiN])+Pi_b[0][AngN][PhiN]);
	}
	Co_U[EnrN][AngN] = ((i*Co_U_m[0][AngN])+Co_U_b[0][AngN]);
	Pi_U[EnrN][AngN] = ((i*Pi_U_m[0][AngN])+Pi_U_b[0][AngN]);
      }
      Co_max[EnrN] = ((i*Co_max_m[0])+Co_max_b[0]);
      Co_U_max[EnrN] = ((i*Co_U_max_m[0])+Co_U_max_b[0]);
      Co_tot[EnrN] = ((i*Co_tot_m[0])+Co_tot_b[0]);

      Pi_max[EnrN] = ((i*Pi_max_m[0])+Pi_max_b[0]);
      Pi_U_max[EnrN] = ((i*Pi_U_max_m[0])+Pi_U_max_b[0]);
      Pi_tot[EnrN] = ((i*Pi_tot_m[0])+Pi_tot_b[0]);
    }

    proc_weight[0][EnrN] = Co_tot[EnrN]*targ_weight[0];
    proc_weight[1][EnrN] = 1000*Pi_tot[EnrN]*targ_weight[0];
    proc_weight[2][EnrN] = Co_tot[EnrN]*targ_weight[1];
    proc_weight[3][EnrN] = 1000*Pi_tot[EnrN]*targ_weight[1];
    
    proc_tot[EnrN] = (proc_weight[0][EnrN]+proc_weight[1][EnrN]+proc_weight[2][EnrN]+proc_weight[3][EnrN]);

    hproc_i[EnrN] = new TH1F(Proc_i_Name[EnrN],"Process weight initial",4,1,5);
    hproc_f[EnrN] = new TH1F(Proc_f_Name[EnrN],"Process weight final",4,1,5);

    for(j=0; j<4; j++){
      hproc_i[EnrN]->GetXaxis()->SetBinLabel((j+1),proc_text[j]);
      hproc_f[EnrN]->GetXaxis()->SetBinLabel((j+1),proc_text[j]);

      hproc_i[EnrN]->Fill((j+1), proc_weight[j][EnrN]);
    }
  }

  c1->cd();
  hco[1]->Draw("surf2");
  gPad->SetPhi(315);
  gPad->Modified();
  gPad->Update();
  c1->SaveAs("out/Comp_Sigma.eps");
  
  c2->cd();
  hpi[1]->Draw("surf2");
  gPad->SetPhi(315);
  gPad->Modified();
  gPad->Update();
  c2->SaveAs("out/Pi0_Sigma.eps");

  c3->cd();
  hco_U[1]->Draw();
  c3->SaveAs("out/IncohComp_Sigma.eps");

  c4->cd();
  hpi_U[1]->Draw();
  c4->SaveAs("out/IncohPi0_Sigma.eps");

  cout << endl << "Running..." << endl << endl;

  i=0;
  j=0;
  k=0;
  m=0;

  TF1 *f1 = new TF1("f1", "1/x", Energy[0], Energy[2]);
  
  timer.Start();

  while(i<NEvn){
    qi_e = f1->GetRandom();
    EnrN = (qi_e-Energy[0]);

    if(RandSel){
      proc = hproc_i[EnrN]->GetRandom();
    }

    hproc_f[EnrN]->Fill(proc);

    if(proc==1 || proc==3) qf_m = 0;
    else if(proc==2 || proc==4) qf_m = kMPI0_MEV;
    else{
      cout << "Invalid process" << endl;
      continue;
    }
    
    // Non-target window vertex (everthing else)
    vtx_z = 2.0*(-0.5 + gRandom->Rndm());
    
    // The while statement cuts off the gaussian xy values of vertex 
    // position so that they are inside the target.
    while(sqrt(Sqr(vtx_x = gRandom->Gaus(0,0.5))
	       +Sqr(vtx_y = gRandom->Gaus(0,0.5))) > 0.5);
    vtx.SetXYZ(vtx_x, vtx_y, vtx_z);
    
    qi.SetPxPyPzE(0,0,qi_e,qi_e);
    ki.SetPxPyPzE(0,0,0,ki_e);
    
    ptot = qi + ki;
    cm_to_lab = ptot.BoostVector();
    lab_to_cm = -ptot.BoostVector();

    qi_cm = qi;
    qi_cm.Boost(lab_to_cm);
    ki_cm = ki;
    ki_cm.Boost(lab_to_cm);
    ptot_cm = ptot;
    ptot_cm.Boost(lab_to_cm);
    
    qth_cm = acos(-1+2*gRandom->Rndm());
    qph_cm = kPI*(-1+2*gRandom->Rndm());
    
    qf_cm.SetPxPyPzE(1,1,1,0);
    qf_cm.SetRho(1);
    qf_cm.SetTheta(qth_cm);
    qf_cm.SetPhi(qph_cm);

    qf_e = ((ptot_cm.M2()+Sqr(qf_m)-Sqr(kf_m))/(2*ptot_cm.M()));
    qf_p = sqrt(Sqr(qf_e)-Sqr(qf_m));

    qf_cm = qf_cm*qf_p;
    qf_cm.SetE(qf_e);
    
    qf = qf_cm;
    qf.Boost(cm_to_lab);
    
    qth = qf.Theta();
    qph = qf.Phi();

    weight = gRandom->Rndm();
    
    if(proc==1){
      AngN = (qth*kR2D);
      PhiN = fabs(qph*kR2D);
      cross = Co[EnrN][AngN][PhiN];
      cross_max = Co_max[EnrN];
      /*
      if(qi_e<BeamE){
	cross = ((qi_e*Co_m[0][AngN][PhiN])+Co_b[0][AngN][PhiN]);
	cross_max = ((qi_e*Co_max_m[0])+Co_max_b[0]);
      }
      else{
	cross = ((qi_e*Co_m[1][AngN][PhiN])+Co_b[1][AngN][PhiN]);
	cross_max = ((qi_e*Co_max_m[1])+Co_max_b[1]);
      }
      */
      /*
      if(hco[1]->GetBinContent((AngN+1),(PhiN+181))<=(weight*Co_max[1])){
	j++;
	continue;
      }
      */
    }

    else if(proc==2){
      AngN = (qth_cm*kR2D);
      PhiN = fabs(qph_cm*kR2D);
      cross = Pi[EnrN][AngN][PhiN];
      cross_max = Pi_max[EnrN];
      /*
      if(qi_e<BeamE){
	cross = ((qi_e*Pi_m[0][AngN][PhiN])+Pi_b[0][AngN][PhiN]);
	cross_max = ((qi_e*Pi_max_m[0])+Pi_max_b[0]);
      }
      else{
	cross = ((qi_e*Pi_m[1][AngN][PhiN])+Pi_b[1][AngN][PhiN]);
	cross_max = ((qi_e*Pi_max_m[1])+Pi_max_b[1]);
      }
      */
      /*
      if(hpi[1]->GetBinContent((AngN+1),(PhiN+181))<=(weight*Pi_max[1])){
	j++;
	continue;
      }
      */
    }
    
    else if(proc==3){
      AngN = (qth*kR2D);
      PhiN = fabs(qph*kR2D);
      cross = Co_U[EnrN][AngN];
      cross_max = Co_U_max[EnrN];
      /*
      if(qi_e<BeamE){
	cross = ((qi_e*Co_U_m[0][AngN])+Co_U_b[0][AngN]);
	cross_max = ((qi_e*Co_U_max_m[0])+Co_U_max_b[0]);
      }
      else{
	cross = ((qi_e*Co_U_m[1][AngN])+Co_U_b[1][AngN]);
	cross_max = ((qi_e*Co_U_max_m[1])+Co_U_max_b[1]);
      }
      */
      /*
      if(hco[1]->GetBinContent((AngN+1),(PhiN+181))<=(weight*Co_max[1])){
	j++;
	continue;
      }
      */
    }

    else if(proc==4){
      AngN = (qth_cm*kR2D);
      PhiN = fabs(qph_cm*kR2D);
      cross = Pi_U[EnrN][AngN];
      cross_max = Pi_U_max[EnrN];
      /*
      if(qi_e<BeamE){
	cross = ((qi_e*Pi_U_m[0][AngN])+Pi_U_b[0][AngN]);
	cross_max = ((qi_e*Pi_U_max_m[0])+Pi_U_max_b[0]);
      }
      else{
	cross = ((qi_e*Pi_U_m[1][AngN])+Pi_U_b[1][AngN]);
	cross_max = ((qi_e*Pi_U_max_m[0])+Pi_U_max_b[0]);
      }
      */
      /*
      if(hpi[1]->GetBinContent((AngN+1),(PhiN+181))<=(weight*Pi_max[1])){
	j++;
	continue;
      }
      */
    }

    else{
      cout << "Invalid process" << endl;
      continue;
    }

    if(cross<=(weight*cross_max)){
      j++;
      continue;
    }
    
    if(proc==3 || proc==4){
      SpecModel(ki);
      
      ptot = qi + ki;
      cm_to_lab = ptot.BoostVector();
      lab_to_cm = -ptot.BoostVector();
      
      qi_cm = qi;
      qi_cm.Boost(lab_to_cm);
      ki_cm = ki;
      ki_cm.Boost(lab_to_cm);
      ptot_cm = ptot;
      ptot_cm.Boost(lab_to_cm);
      
      qf_cm = qf;
      qf_cm.Boost(lab_to_cm);
      
      qth_cm = qf_cm.Theta();
      qph_cm = qf_cm.Phi();
      
      qf_cm.SetPxPyPzE(1,1,1,0);
      qf_cm.SetRho(1);
      qf_cm.SetTheta(qth_cm);
      qf_cm.SetPhi(qph_cm);
      
      qf_e = ((ptot_cm.M2()+Sqr(qf_m)-Sqr(kf_m))/(2*ptot_cm.M()));
      qf_p = sqrt(Sqr(qf_e)-Sqr(qf_m));
      
      qf_cm = qf_cm*qf_p;
      qf_cm.SetE(qf_e);
    }
     
    if(qf_e<qf_m){
      m++;
      continue;
    }

    kf_e = (qi_cm.E()+ki_cm.E()-qf_cm.E());
    kf_p = qf_p;
    
    kf_cm = -qf_cm;
    kf_cm.SetE(kf_e);

    qf = qf_cm;
    qf.Boost(cm_to_lab);
    kf = kf_cm;
    kf.Boost(cm_to_lab);

    if(qf_cm.Px()+kf_cm.Px()!=0 || qf_cm.Px()+kf_cm.Px()!=0 ||
       qf_cm.Px()+kf_cm.Px()!=0){
      k++;
      continue;
    }

    if(proc==2 || proc==4){
      pi_to_lab = qf.BoostVector();
      lab_to_pi = -qf.BoostVector();
      
      pth_pi = acos(-1+2*gRandom->Rndm());
      pph_pi = kPI*(-1+2*gRandom->Rndm());
      
      p1_pi.SetPxPyPzE(1,1,1,0);
      p1_pi.SetRho(1);
      p1_pi.SetTheta(pth_pi);
      p1_pi.SetPhi(pph_pi);
      
      p1_pi = p1_pi*(qf_m/2);
      p2_pi = -p1_pi;
      
      p1_pi.SetE(qf_m/2);
      p2_pi.SetE(qf_m/2);
      
      p1 = p1_pi;
      p1.Boost(pi_to_lab);
      p2 = p2_pi;
      p2.Boost(pi_to_lab);
      
      if(p1_pi.Px()+p2_pi.Px()!=0 || p1_pi.Px()+p2_pi.Px()!=0 ||
	 p1_pi.Px()+p2_pi.Px()!=0){
	k++;
	continue;
      }

    }

    if(     (qf.Theta()*kR2D)<=20 ) qfDet = 0;
    else if((qf.Theta()*kR2D)<=160) qfDet = 1;
    else if((qf.Theta()*kR2D)<=180) qfDet = 2;
    else{
      cout << "Final angles invalid." << endl;
      continue;
    }

    if(     (kf.Theta()*kR2D)<=20 ) kfDet = 0;
    else if((kf.Theta()*kR2D)<=160) kfDet = 1;
    else if((kf.Theta()*kR2D)<=180) kfDet = 2;
    else{
      cout << "Final angles invalid." << endl;
      continue;
    }

    NComp[qfDet][kfDet]++;

    if(proc==2 || proc==4){
      if(     (p1.Theta()*kR2D)<=20 ) p1Det = 0;
      else if((p1.Theta()*kR2D)<=160) p1Det = 1;
      else if((p1.Theta()*kR2D)<=180) p1Det = 2;
      else{
	cout << "Final angles invalid." << endl;
	continue;
      }
      if(     (p2.Theta()*kR2D)<=20 ) p2Det = 0;
      else if((p2.Theta()*kR2D)<=160) p2Det = 1;
      else if((p2.Theta()*kR2D)<=180) p2Det = 2;
      else{
	cout << "Final angles invalid." << endl;
	continue;
      }

      NPi0[kfDet][p1Det][p2Det]++;
    }

    hqi[proc][0]->Fill(qi.E()-qi.M());
    hqi[proc][1]->Fill(qi.Px());
    hqi[proc][2]->Fill(qi.Py());
    hqi[proc][3]->Fill(qi.Pz());
    hqi[proc][4]->Fill(qi.M());
    hqi[proc][5]->Fill(qi.Theta()*kR2D);

    if(proc==3 || proc==4){
      hki[proc][0]->Fill(ki.E()-ki.M());
      hki[proc][1]->Fill(ki.Px());
      hki[proc][2]->Fill(ki.Py());
      hki[proc][3]->Fill(ki.Pz());
      hki[proc][4]->Fill(ki.M());
      hki[proc][5]->Fill(ki.Theta()*kR2D);
    }

    hqf[proc][0]->Fill(qf.E()-qf.M());
    hqf[proc][1]->Fill(qf.Px());
    hqf[proc][2]->Fill(qf.Py());
    hqf[proc][3]->Fill(qf.Pz());
    hqf[proc][4]->Fill(qf.M());
    hqf[proc][5]->Fill(qf.Theta()*kR2D);

    hkf[proc][0]->Fill(kf.E()-kf.M());
    hkf[proc][1]->Fill(kf.Px());
    hkf[proc][2]->Fill(kf.Py());
    hkf[proc][3]->Fill(kf.Pz());
    hkf[proc][4]->Fill(kf.M());
    hkf[proc][5]->Fill(kf.Theta()*kR2D);
 
    if(proc==2 || proc==4){
      hp1[proc][0]->Fill(p1.E()-p1.M());
      hp1[proc][1]->Fill(p1.Px());
      hp1[proc][2]->Fill(p1.Py());
      hp1[proc][3]->Fill(p1.Pz());
      hp1[proc][4]->Fill(p1.M());
      hp1[proc][5]->Fill(p1.Theta()*kR2D);
 
      hp2[proc][0]->Fill(p2.E()-p2.M());
      hp2[proc][1]->Fill(p2.Px());
      hp2[proc][2]->Fill(p2.Py());
      hp2[proc][3]->Fill(p2.Pz());
      hp2[proc][4]->Fill(p2.M());
      hp2[proc][5]->Fill(p2.Theta()*kR2D);

      hdeang->Fill((p1.Theta()*kR2D),(p2.Theta()*kR2D));
    }

    if((kf.E()-kf.M())>40) hangle->Fill((qf.Theta()*kR2D),(kf.Theta()*kR2D));
    hmissm->Fill((kf.E()-kf.M()),(qi.E()-(qf.E()+(kf.E()-kf.M()))));

    // Interaction vertex position
    var[0] = vtx.X();
    var[1] = vtx.Y();
    var[2] = vtx.Z();
    
    // Incident photon beam
    var[3] = qi.Px()/qi.P();
    var[4] = qi.Py()/qi.P();
    var[5] = qi.Pz()/qi.P();
    var[6] = qi.P()/1000;
    var[7] = qi.E()/1000;
    
    // Recoil nucleus
    var[8] = kf.Px()/kf.P();
    var[9] = kf.Py()/kf.P();
    var[10] = kf.Pz()/kf.P();
    var[11] = kf.P()/1000;
    var[12] = kf.E()/1000;
    
    // Scattered Particle
    var[13] = qf.Px()/qf.P();
    var[14] = qf.Py()/qf.P();
    var[15] = qf.Pz()/qf.P();
    var[16] = qf.P()/1000;
    var[17] = qf.E()/1000;

    if(proc==2 || proc==4){
      // Photon 1
      var[18] = p1.Px()/p1.P();
      var[19] = p1.Py()/p1.P();
      var[20] = p1.Pz()/p1.P();
      var[21] = p1.P()/1000;
      var[22] = p1.E()/1000;
      
      // Photon 2
      var[23] = p2.Px()/p2.P();
      var[24] = p2.Py()/p2.P();
      var[25] = p2.Pz()/p2.P();
      var[26] = p2.P()/1000;
      var[27] = p2.E()/1000;
    }

    // Fill ntuple
    h1->Fill(var);

    i++;
  }

  timer.Stop();
  cout << i << " events in " << timer.RealTime() << " sec." << endl << endl;

  cout << "Recoil particle:" << endl;
  cout << kf.Px() << "  \t" << kf.Py() << "  \t" << kf.Pz() << endl;
  cout << kf.Theta()*kR2D << "  \t" << kf.Phi()*kR2D << "  \t" << kf.Rho() << endl;

  cout << "Scattered particle:" << endl;
  cout << qf.Px() << "  \t" << qf.Py() << "  \t" << qf.Pz() << endl;
  cout << qf.Theta()*kR2D << "  \t" << qf.Phi()*kR2D << "  \t" << qf.Rho() << endl;

  if(proc==2 || proc==4){
    cout << "Decay Photon 1:" << endl;
    cout << p1.Px() << "  \t" << p1.Py() << "  \t" << p1.Pz() << endl;
    cout << p1.Theta()*kR2D << "  \t" << p1.Phi()*kR2D << "  \t" << p1.Rho() << endl;
    
    cout << "Decay Photon 2:" << endl;
    cout << p2.Px() << "  \t" << p2.Py() << "  \t" << p2.Pz() << endl;
    cout << p2.Theta()*kR2D << "  \t" << p2.Phi()*kR2D << "  \t" << p2.Rho() << endl;
  }

  cout << endl << "\t\t\tProton" << endl;
  cout << "\t\tTAPS\tCB\tOut" << endl;
  cout << "\tTAPS\t" << NComp[0][0] << "\t" << NComp[0][1] << "\t" << NComp[0][2] << endl;
  cout << "Scatter\tCB\t" << NComp[1][0] << "\t" << NComp[1][1] << "\t" << NComp[1][2] << endl;
  cout << "\tOut\t" << NComp[2][0] << "\t" << NComp[2][1] << "\t" << NComp[2][2] << endl;

  if(proc==2 || proc==4){
    cout << endl << "\t\t\t\t\t\tProton" << endl;
    cout << "\t\t\tTAPS\t\t\tCB\t\t\tOut" << endl;
    cout << "\t\t\tPhoton1\t\t\tPhoton1\t\t\tPhoton1" << endl;
    cout << "\t\tTAPS\tCB\tOut\tTAPS\tCB\tOut\tTAPS\tCB\tOut" << endl;
    cout << "\tTAPS\t" << NPi0[0][0][0] << "\t" << NPi0[0][0][1] << "\t" << NPi0[0][0][2] << "\t" << NPi0[1][0][0] << "\t" << NPi0[1][0][1] << "\t" << NPi0[1][0][2] << "\t" << NPi0[2][0][0] << "\t" << NPi0[2][0][1] << "\t" << NPi0[2][0][2] << endl;
    cout << "Photon2\tCB\t" << NPi0[0][1][0] << "\t" << NPi0[0][1][1] << "\t" << NPi0[0][1][2] << "\t" << NPi0[1][1][0] << "\t" << NPi0[1][1][1] << "\t" << NPi0[1][1][2] << "\t" << NPi0[2][1][0] << "\t" << NPi0[2][1][1] << "\t" << NPi0[2][1][2] << endl;
    cout << "\tOut\t" << NPi0[0][2][0] << "\t" << NPi0[0][2][1] << "\t" << NPi0[0][2][2] << "\t" << NPi0[1][2][0] << "\t" << NPi0[1][2][1] << "\t" << NPi0[1][2][2] << "\t" << NPi0[2][2][0] << "\t" << NPi0[2][2][1] << "\t" << NPi0[2][2][2] << endl;
  }

  cout << endl << i << " accepted, " << j << " discarded, " << k << " unequal in cm, " << m << " below threshold." << endl;

  TFile hfile("out/test.root", "RECREATE", "MC_Ntuple_File");

  if(RandSel){
    hproc_i[10]->Write();
    hproc_f[10]->Write();
    for(i=0; i<6; i++){
      for(j=1; j<5; j++){
	hqi[0][i]->Add(hqi[j][i]);
	hki[0][i]->Add(hki[j][i]);
	hqf[0][i]->Add(hqf[j][i]);
	hkf[0][i]->Add(hkf[j][i]);
	hp1[0][i]->Add(hp1[j][i]);
	hp2[0][i]->Add(hp2[j][i]);

	hqi[j][i]->SetFillColor(j+1);
	sqi[i]->Add(hqi[j][i]);
	hki[j][i]->SetFillColor(j+1);
	ski[i]->Add(hki[j][i]);
	hqf[j][i]->SetFillColor(j+1);
	sqf[i]->Add(hqf[j][i]);
	hkf[j][i]->SetFillColor(j+1);
	skf[i]->Add(hkf[j][i]);
	hp1[j][i]->SetFillColor(j+1);
	sp1[i]->Add(hp1[j][i]);
	hp2[j][i]->SetFillColor(j+1);
	sp2[i]->Add(hp2[j][i]);
      }
      hqi[0][i]->Write();
      hki[0][i]->Write();
      hqf[0][i]->Write();
      hkf[0][i]->Write();
      hp1[0][i]->Write();
      hp2[0][i]->Write();

      sqi[i]->Write();
      ski[i]->Write();
      sqf[i]->Write();
      skf[i]->Write();
      sp1[i]->Write();
      sp2[i]->Write();
    }
  }
  else{
    for(i=0; i<6; i++){
      if((hqi[proc][i]->GetMaximum())>0) hqi[proc][i]->Write();
      if((hki[proc][i]->GetMaximum())>0) hki[proc][i]->Write();
      if((hqf[proc][i]->GetMaximum())>0) hqf[proc][i]->Write();
      if((hkf[proc][i]->GetMaximum())>0) hkf[proc][i]->Write();
      if((hp1[proc][i]->GetMaximum())>0) hp1[proc][i]->Write();
      if((hp2[proc][i]->GetMaximum())>0) hp2[proc][i]->Write();
    }
  }

  hangle->Write();
  hdeang->Write();
  hmissm->Write();
  
  h1->Write();

  hfile.Close();
  
  return 0;
}
