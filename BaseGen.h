// Compton Scattering and Pi0P/PiPN Photoproduction Event Generator
// Designed for use with the MAMI A2 Geant4 Simulation
//
// Base class for event generators
//
// Author - P. Martel
// Version - 02 November 2012

class BaseGen {
 protected:
  TString sProcName;
  TString sTargName;
  TString sBaseName;
  Bool_t bIncoh;
  Bool_t bCoher;
  Bool_t bIsotW;
  Int_t iBeamLo;
  Int_t iBeamHi;
  Int_t iBeamSt;
  Int_t iAngSt;
  Int_t iNEnr;
  Int_t iNAng;
  Int_t iNPhi;
  Float_t fBeamT;
  Float_t fBeamL;
  Float_t fBeamP;
  Float_t fTargT;
  Float_t fTargL;
  Float_t fTargP;
  Float_t *fEnr;
  Float_t ***fPar;
  TH3F *hCrossSec;
  TH1F *hCrossMax;
  TH1F *hCrossTot;
  TLorentzVector ptot;
  TVector3 cm_to_lab, lab_to_cm;
  Float_t vtx_x, vtx_y, vtx_z;
 public:
  BaseGen(TString, TString, TString, Int_t, Int_t);
  ~BaseGen();
  void SetPol(Float_t, Float_t, Float_t, Float_t, Float_t, Float_t);
  void InitBase(const char*);
  void CrossGen();
  Float_t GetCross(Int_t, Int_t, Int_t);
  void Collision2B(BasePart&, BasePart&, BasePart&, BasePart&);
  void Decay2B(BasePart&, BasePart&, BasePart&);
  void SpecModel(BasePart&);
  void InitNtuple(Int_t, Int_t*);
  void SaveNtuple(TString);
  void SaveTree(TString);
  void NewVertex();
  Float_t Sqr(Float_t x){return(x*x);};
  Bool_t Reject(Float_t, Float_t, Float_t);
  TNtuple *h1;
  TTree *t1;
  Float_t var[100];
};

BaseGen::BaseGen(TString name, TString target, TString base, Int_t beamlo, Int_t beamhi){

  // Set reaction information

  sProcName = name;
  sTargName = target;
  sBaseName = base;

  if(sProcName.Contains("Incoherent")) bIncoh = kTRUE;
  else bIncoh = kFALSE;

  if(sProcName.Contains("Coherent")) bCoher = kTRUE;
  else bCoher = kFALSE;

  if(sProcName.Contains("Isot")) bIsotW = kTRUE;
  else bIsotW = kFALSE;

  iBeamLo = beamlo;
  iBeamHi= beamhi;
  iBeamSt = (beamhi-beamlo);

  // Begin with an unpolarized state

  fBeamT = 0;
  fBeamL = 0;
  fBeamP = 0;
  fTargT = 0;
  fTargL = 0;
  fTargP = 0;

  cout << "Constructing base generator" << endl;

  // Initialize cross section histograms, bin sizes will be adjusted later
  /*
  hCrossSec = new TH3F(sProcName,sProcName,iBeamSt+1,iBeamLo,(iBeamHi+1),181,0,181,181,0,181);
  hCrossMax = new TH1F(sProcName+"Max",(sProcName+" - max"),iBeamSt+1,iBeamLo,(iBeamHi+1));
  hCrossTot = new TH1F(sProcName+"Tot",(sProcName+" - tot"),iBeamSt+1,iBeamLo,(iBeamHi+1));
  */
  hCrossSec = new TH3F("hCrossSec",sProcName,iBeamSt+1,iBeamLo,(iBeamHi+1),181,0,181,181,0,181);
  hCrossMax = new TH1F("hCrossMax",(sProcName+" - max"),iBeamSt+1,iBeamLo,(iBeamHi+1));
  hCrossTot = new TH1F("hCrossTot",(sProcName+" - tot"),iBeamSt+1,iBeamLo,(iBeamHi+1));
  
  // Initialize output tree

  t1 = new TTree("OutTree","Generator kinematics tree");

};

BaseGen::~BaseGen(){

  cout << "Deleting base generator" << endl << endl;

  delete fEnr;
  delete fPar;
  delete hCrossSec;
  delete hCrossMax;
  delete hCrossTot;

};

void BaseGen::SetPol(Float_t fBTin, Float_t fBLin, Float_t fBPin, Float_t fTTin, Float_t fTLin, Float_t fTPin){

  // Set polarization values

  fBeamT = fBTin;
  fBeamL = fBLin;
  fBeamP = fBPin;
  fTargT = fTTin;
  fTargL = fTLin;
  fTargP = fTPin;

};

void BaseGen::InitBase(const char* cComFile){

  // Initialization for Normal and Incoherent reactions, reads in data files
  // from the 'par' directory to construct cross section tables

  const Int_t iNParC = 9;

  Int_t iBeamE = 0, iParN = 0, iEnrN = 0, iAngN = 0, iNPar = iNParC;

  Float_t *fParT;

  fParT = new Float_t[iNParC];
  for(iParN=0; iParN<iNPar; iParN++) fParT[iParN] = 0;

  Int_t iLength = strlen(cComFile);

  gSystem->cd("par");
    
  TSystemDirectory *sdParDir = new TSystemDirectory("parameters","./");
  TList *lParList = sdParDir->GetListOfFiles();
  lParList->Sort();
  TIter ParNext(lParList);
  TSystemFile *sfParFile;
  const char* cParName;

  // Loop through the available parameter files once to determine how many exist
  // and how many angles they each have
  
  while((sfParFile=(TSystemFile*)ParNext())){
    cParName = sfParFile->GetName();
    if(strncmp(cParName,cComFile,iLength) == 0){
      iBeamE = atoi(&cParName[iLength]);
      if(iBeamE >= iBeamLo && iBeamE <= iBeamHi){
	iAngN = 0;
	ifstream fin(cParName);
	while(!fin.eof()){
	  for(iParN=0; iParN<iNPar; iParN++) fin >> fParT[iParN];
	  iAngN++;
	}
	fin.close();
	iAngN--;
	if((iEnrN>0) && (iAngN!=iNAng)){
	  cout << "Files have different numbers of angles" << endl;
	  gSystem->Exit(0);
	}
	iNAng = iAngN;
	iEnrN++;
      }
    }
  }

  iNEnr = iEnrN;
  iNAng = iAngN;
  iNPhi = ((2*iAngN)-1);

  const Int_t iNEnrC = iEnrN;
  const Int_t iNAngC = iAngN;

  // Determine step sizes for energy and angle
  
  iBeamSt = ((iBeamHi-iBeamLo)/(iNEnr-1));
  iAngSt = (180/(iNAng-1));

  // Create arrays to hold input from parameter files

  fPar = new Float_t**[iNParC];
  for(iParN=0; iParN<iNPar; iParN++){
    fPar[iParN] = new Float_t*[iNEnrC];
    for(iEnrN=0; iEnrN<iNEnr; iEnrN++){
      fPar[iParN][iEnrN] = new Float_t[iNAngC];
      for(iAngN=0; iAngN<iNAng; iAngN++){
	fPar[iParN][iEnrN][iAngN] = 0;
      }
    }
  }

  fEnr = new Float_t[iNEnrC];
  for(iEnrN=0; iEnrN<iNEnr; iEnrN++){
    fEnr[iEnrN] = 0;
  }

  iEnrN = 0;

  ParNext.Reset();

  // Loop through the available parameter files again to fill the arrays
  
  while((sfParFile=(TSystemFile*)ParNext())){
    cParName = sfParFile->GetName();
    if(strncmp(cParName,cComFile,iLength) == 0){
      iBeamE = atoi(&cParName[iLength]);
      if(iBeamE >= iBeamLo && iBeamE <= iBeamHi){
	fEnr[iEnrN] = iBeamE;
	iAngN = 0;
	ifstream fin(cParName);
	while(!fin.eof()){
	  for(iParN=0; iParN<iNPar; iParN++) fin >> fPar[iParN][iEnrN][iAngN];
	  iAngN++;
	}
	fin.close();
	iEnrN++;
      }
    }
  }

  // Check to make sure the files match, and that the step sizes are uniform

  for(iEnrN=0; iEnrN<iNEnr; iEnrN++){
    for(iAngN=1; iAngN<iNAng; iAngN++){
      if((fPar[0][iEnrN][iAngN]-fPar[0][iEnrN][iAngN-1])!=iAngSt){
	cout << "Parameter files have uneven angle steps" << endl;
	gSystem->Exit(0);
      }
    }
  }
  for(iEnrN=1; iEnrN<iNEnr; iEnrN++){
    if((fEnr[iEnrN]-fEnr[iEnrN-1])!=iBeamSt){
      cout << "Parameter files have uneven energy steps" << endl;
      gSystem->Exit(0);
    }
    for(iAngN=0; iAngN<iNAng; iAngN++){
      if(fPar[0][iEnrN][iAngN]!=fPar[0][iEnrN-1][iAngN]){
	cout << "Parameter files have different angles" << endl;
	gSystem->Exit(0);
      }
    }
  }

  cout << endl << "Found " << iNEnr << " energies (parameter files) with " << iNAng << " angle entries" << endl << "Giving " << iBeamSt << " MeV and " << iAngSt << " deg steps" << endl << endl;

  // Reset bin settings for the cross section histograms

  hCrossSec->SetBins(iNEnr,iBeamLo-0.5*iBeamSt,iBeamHi+0.5*iBeamSt,iNAng,0-0.5*iAngSt,180+0.5*iAngSt,iNPhi,-180-0.5*iAngSt,180+0.5*iAngSt);
  hCrossMax->SetBins(iNEnr,iBeamLo-0.5*iBeamSt,iBeamHi+0.5*iBeamSt);
  hCrossTot->SetBins(iNEnr,iBeamLo-0.5*iBeamSt,iBeamHi+0.5*iBeamSt);
  
  delete sdParDir;
  delete lParList;
  delete sfParFile;
  
  gSystem->cd("..");

};

void BaseGen::CrossGen(){

  Int_t iEnrN = 0, iAngN = 0, iPhiN = 0;
  Float_t fPhi, fAngL, fAngH;

  Float_t fSolAng = 0, fSolAngTot = 0, fCVal = 0, fCMax = 0, fCTot = 0;
  Float_t fCVll, fCVlu, fCVul, fCVuu, fCVav;

  // Construct cross sections and fill the corresponding histogram

  for(iEnrN=0; iEnrN<iNEnr; iEnrN++){
    fSolAngTot = 0;
    fCMax = 0;
    fCTot = 0;
    for(iAngN=0; iAngN<iNAng; iAngN++){
      fAngL = fPar[0][iEnrN][iAngN];
      fAngH = fPar[0][iEnrN][iAngN+1];

      if(iAngN < (iNAng-1)) fSolAng = (cos(fAngL*kD2R)-cos(fAngH*kD2R));

      for(iPhiN=0; iPhiN<iNPhi; iPhiN++){
	fPhi = (-180+(iPhiN*iAngSt));

	fCVal = GetCross(iEnrN, iAngN, iPhiN);

	hCrossSec->Fill(fEnr[iEnrN], fAngL, fPhi, fCVal);
	
	// Check for maximum cross section value
	if(fCVal > fCMax) fCMax = fCVal;

	// Compute average cross section over angular bin
	if((iAngN < (iNAng-1)) && (iPhiN < (iNPhi-1))){
	  fCVll = GetCross(iEnrN, iAngN, iPhiN);
	  fCVlu = GetCross(iEnrN, iAngN, iPhiN+1);
	  fCVul = GetCross(iEnrN, iAngN+1, iPhiN);
	  fCVuu = GetCross(iEnrN, iAngN+1, iPhiN+1);
	  fCVav = ((fCVll+fCVlu+fCVul+fCVuu)/4.);	  

	  // Sum up total cross section and solid angle (which should be 4pi)
	  fCTot += (fCVal*iAngSt*kD2R*fSolAng);
	  fSolAngTot += (iAngSt*kD2R*fSolAng);
	}
      }
    }

    hCrossMax->Fill(fEnr[iEnrN], fCMax);
    hCrossTot->Fill(fEnr[iEnrN], fCTot);

    cout << fEnr[iEnrN] << "  \t\t" << fCMax << "   \t" << fCTot;
    if(!(TMath::AreEqualRel(fSolAngTot,(4*TMath::Pi()),0.0001))) cout << "\t\tError - Sol Ang = " << fSolAngTot << " sr";
    cout << endl;
    
  }

};

Float_t BaseGen::GetCross(Int_t iEnrN, Int_t iAngN, Int_t iPhiN){

  // Determine cross section from read-in parameters
  
  Float_t fCVal = 0;
  Float_t fPhi = (-180+(iPhiN*iAngSt));
  
  Float_t fParU = fPar[1][iEnrN][iAngN];
  Float_t fParS = fPar[2][iEnrN][iAngN];
  Float_t fParT = fPar[3][iEnrN][iAngN];
  Float_t fParP = fPar[4][iEnrN][iAngN];
  Float_t fParG = fPar[5][iEnrN][iAngN];
  Float_t fParH = fPar[6][iEnrN][iAngN];
  Float_t fParE = fPar[7][iEnrN][iAngN];
  Float_t fParF = fPar[8][iEnrN][iAngN];

  fCVal = fParU*(1-(fBeamT*fParS*cos(2*(fBeamP-fPhi)*kD2R)));
  
  fCVal += fParU*(((fBeamL*fParF)-(fBeamT*fParH*sin(2*(fBeamP-fPhi)*kD2R)))*
		  fTargT*cos((fTargP-fPhi)*kD2R));
  fCVal += fParU*((fParT-(fBeamT*fParP*cos(2*(fBeamP-fPhi)*kD2R)))*
		  fTargT*sin((fTargP-fPhi)*kD2R));
  fCVal += fParU*(((-fBeamL*fParE)+(fBeamT*fParG*sin(2*(fBeamP-fPhi)*kD2R)))*
		  fTargL);
  
  return fCVal;

};

void BaseGen::Collision2B(BasePart& qi,BasePart& ki,BasePart& qf,BasePart& kf){

  // Two body collision kinematics

  Float_t ener, mom;

  // Total initial state and boost vectors
  
  ptot = qi.P4 + ki.P4;
  cm_to_lab = ptot.BoostVector();
  lab_to_cm = -ptot.BoostVector();

  // Boost initial particles to CM frame
  
  qi.BoostCM(lab_to_cm);
  ki.BoostCM(lab_to_cm);

  // Determine energy and momentum of 'scattered' particle
  
  ener = ((ptot.M2()+Sqr(qf.Mass)-Sqr(kf.Mass))/(2*ptot.M()));
  mom = TMath::Sqrt(Sqr(ener)-Sqr(qf.Mass));

  // If the energy of the 'scattered' particle has NOT been previously set
  // (therefore still zero) set its energy and momentum (magnitude) but
  // let it have an isotropic distribution in CM

  if(qf.Ener==0){  
    qf.SetP4CM(ener,mom);
    qf.BoostLab(cm_to_lab);
  }

  // If the energy of the 'scattered' particle HAS been previously set
  // (through an incoherent process) keep the already chosen direction
  // but adjust the energy and momentum (magnitude)

  else{
    qf.BoostCM(lab_to_cm);
    qf.SetP4CM(ener,mom,qf.ThetaCM,qf.PhiCM);
    qf.BoostLab(cm_to_lab);
  }

  // Determine energy and momentum of recoil particle

  ener = (qi.EnerCM + ki.EnerCM - qf.EnerCM);

  kf.SetP4CM(ener,mom,(180-qf.ThetaCM),qf.PhiCM);
  kf.RotateCM(180);
  kf.BoostLab(cm_to_lab);
  
};

void BaseGen::Decay2B(BasePart& k,BasePart& p1,BasePart& p2){

  // Two body decay kinematics

  Float_t ener, mom;

  // Total initial state and boost vectors
  
  ptot = k.P4;
  cm_to_lab = ptot.BoostVector();
  lab_to_cm = -ptot.BoostVector();

  // Boost initial particle to its CM frame

  // k.BoostCM(lab_to_cm); // Unnecessary, and bad for original CM frame info  
  
  // Determine energy and momentum of first decay particle

  ener = ((ptot.M2()+Sqr(p1.Mass)-Sqr(p2.Mass))/(2*ptot.M()));
  mom = TMath::Sqrt(Sqr(ener)-Sqr(p1.Mass));

  p1.SetP4CM(ener,mom);
  p1.BoostLab(cm_to_lab);

  // Set second decay particle back-to-back with first decay particle in CM

  p2.SetP4CM(ener,mom,(180-p1.ThetaCM),p1.PhiCM);
  p2.RotateCM(180);
  p2.BoostLab(cm_to_lab);

};

void BaseGen::SpecModel(BasePart& ki){

  // Spectral Model for Incoherent reactions
  // Takes input nucleon and determines a proper recoil energy and momentum 

  Float_t kpi, prob, func, ener, N_s, N_p;
  Float_t k_max, A[2], coeff[2], sigma[2], Emin[2], Emax[2], sf_max[2];
  Int_t L;

  // Carbon-12 values

  if(sTargName=="12C"){
    N_s = 0.52; N_p = 1.65;
    k_max = 300.0;
    A[0] = 43.5; A[1] = 13500.0;
    coeff[0] = 0.0; coeff[1] = 2.0;
    sigma[0] = 90.0; sigma[1] = 75.0;
    Emin[0] = 26.0; Emin[1] = 16.0;
    Emax[0] = 50.0; Emax[1] = 26.0;
    sf_max[0] = 0.26; sf_max[1] = 0.93;
  }

  // Oxygen-16 values
  // Note that they're identical to Carbon-12 at the moment, to be updated

  if(sTargName=="16O"){
    N_s = 0.52; N_p = 1.65;
    k_max = 300.0;
    A[0] = 43.5; A[1] = 13500.0;
    coeff[0] = 0.0; coeff[1] = 2.0;
    sigma[0] = 90.0; sigma[1] = 75.0;
    Emin[0] = 26.0; Emin[1] = 16.0;
    Emax[0] = 50.0; Emax[1] = 26.0;
    sf_max[0] = 0.26; sf_max[1] = 0.93;
  }

  if((gRandom->Rndm()) < N_s/(N_s+N_p)) L = 0;
  else L = 1;
  
  prob = 1;
  func = 0;

  while(prob>func){
    kpi = k_max*gRandom->Rndm();
    prob = sf_max[L]*gRandom->Rndm();
    func = A[L]*(TMath::Power((kpi/1000.0),(coeff[L]+2)))*(TMath::Exp(-Sqr(kpi/sigma[L])/2));
  }

  ener = (kMP_MEV-(Emin[L]+(Emax[L]-Emin[L])*gRandom->Rndm()));

  // Return nucleon with chosen energy and momentum

  ki.SetP4Lab(ener, kpi);

};

Bool_t BaseGen::Reject(Float_t fBeamE, Float_t fAng, Float_t fPhi){

  if(bIsotW) return kFALSE;

  // Weighting calculation for event rejection
  
  Bool_t bCheck = kFALSE;

  Float_t fCVal = hCrossSec->Interpolate(fBeamE, fAng, fPhi);
  //Float_t fCMax = hCrossMax->Interpolate(fBeamE);
  Float_t fCMax = hCrossMax->GetMaximum();

  // If the cross section for a given energy and angle is less than, or equal
  // to, the maximum cross section (at that energy) at all energies times a 
  // random number between 0 and 1, then we reject that event.

  if(fCVal <= (fCMax*gRandom->Rndm())) bCheck = kTRUE;

  return bCheck;

};

void BaseGen::InitNtuple(Int_t npart, Int_t* ptag) {

  // Setup Ntuple for MAMI A2 Geant4 simulation
  
  Int_t i, j;
  
  TString pstr[] = {"Px", "Py", "Pz", "Pt", "En"};
  TString beam = "X_vtx:Y_vtx:Z_vtx:Px_bm:Py_bm:Pz_bm:Pt_bm:En_bm";
  TString particles;
  TString names;
  
  for ( i = 0; i < npart; i++) {
    for ( j = 0; j < 5; j++) {
      particles.Append( pstr[j]);
      if ( ( i == (npart-1)) && ( j == 4))
	particles.Append( Form( "_l%02d%02d", i+1, ptag[i]));
      else
	particles.Append( Form( "_l%02d%02d:", i+1, ptag[i]));
    }
  }
  
  names = beam + ":" + particles;

  h1 = new TNtuple("h1", "TMCUserGenerator", names);
  
}

void BaseGen::SaveNtuple(TString sFile){

  // Write Ntuple to file, obviously

  TFile f1(sFile, "RECREATE", "MC_Ntuple_File");

  h1->Write();

  f1.Close();
  
};

void BaseGen::SaveTree(TString sFile){

  // Write tree to file, obviously

  TFile f1(sFile, "RECREATE", "MC_Tree_File");

  t1->Write();

  f1.Close();
  
};

void BaseGen::NewVertex() {

  // Choose a new random vertex for this event, given the dimensions of the
  // frozen spin target
  
  vtx_x = 0.5;
  vtx_y = 0.5;

  while((TMath::Sqrt(Sqr(vtx_x)+Sqr(vtx_y))) > 0.5){
    vtx_x = gRandom->Gaus(0,0.5);
    vtx_y = gRandom->Gaus(0,0.5);
  }

  vtx_z = 2.0*(-0.5 + gRandom->Rndm());
  
  var[0] = vtx_x;
  var[1] = vtx_y;
  var[2] = vtx_z;

}
