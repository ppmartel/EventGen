// Compton Scattering and Pi0 Photoproduction Event Generator
// Designed for use with the MAMI A2 Geant4 Simulation
//
// Pi0 photoproduction event generator class
//
// Author - P. Martel
// Version - 21 June 2011

class Pi0PhotGen : public BaseGen {
 protected:
  BasePart pPhoton, pTarget, pPi0, pRecoil, pDecay1, pDecay2;
  TLorentzVector lvMiss;
  TH2F *hPvP1, *hPvP1Lo, *hPvP1Hi, *hPvP2, *hPvP2Lo, *hPvP2Hi, *hMiM;
 public:
  Pi0PhotGen(TString, TString, Float_t, Int_t, Int_t, Int_t);
  ~Pi0PhotGen();
  void Init(Float_t, Float_t);
  void InitCoher();
  void CrossPol(Float_t, Float_t);
  Bool_t NewEvent(Float_t);
  void Reset();
  void SaveHists(TString);
};

Pi0PhotGen::Pi0PhotGen(TString sName, TString sTarget, Float_t fRecoMass, Int_t iRecoG3id, Int_t beamlo, Int_t beamhi) : BaseGen(sName, sTarget, beamlo, beamhi), pPhoton("Photon",0), pTarget("Target",fRecoMass), pPi0("Pi0",kMPI0_MEV), pRecoil("Recoil",fRecoMass), pDecay1("Decay1",0), pDecay2("Decay2",0) {

  cout << ("Constructing generator for "+sProcName) << endl;

  // Setup Geant4 Ntuple

  Int_t ptag[4] = {iRecoG3id,7,1,1};

  InitNtuple(4,ptag);

  // Create some histograms for testing
  
  hPvP1 = new TH2F("hPvP1","Recoil vs Decay Photon1 Theta",180,0,180,180,0,180);
  hPvP1->GetXaxis()->SetTitle("Photon Theta (deg)");
  hPvP1->GetYaxis()->SetTitle("Recoil Theta (deg)");

  hPvP1Lo = new TH2F("hPvP1Lo","Recoil (E<50 MeV) vs Decay Photon1 Theta",180,0,180,180,0,180);
  hPvP1Lo->GetXaxis()->SetTitle("Photon Theta (deg)");
  hPvP1Lo->GetYaxis()->SetTitle("Recoil Theta (deg)");
  hPvP1Lo->SetMarkerColor(2);

  hPvP1Hi = new TH2F("hPvP1Hi","Recoil (E>=50 MeV) vs Decay Photon1 Theta",180,0,180,180,0,180);
  hPvP1Hi->GetXaxis()->SetTitle("Photon Theta (deg)");
  hPvP1Hi->GetYaxis()->SetTitle("Recoil Theta (deg)");

  hPvP2 = new TH2F("hPvP2","Recoil vs Decay Photon2 Theta",180,0,180,180,0,180);
  hPvP2->GetXaxis()->SetTitle("Photon Theta (deg)");
  hPvP2->GetYaxis()->SetTitle("Recoil Theta (deg)");

  hPvP2Lo = new TH2F("hPvP2Lo","Recoil (E<50 MeV) vs Decay Photon2 Theta",180,0,180,180,0,180);
  hPvP2Lo->GetXaxis()->SetTitle("Photon Theta (deg)");
  hPvP2Lo->GetYaxis()->SetTitle("Recoil Theta (deg)");
  hPvP2Lo->SetMarkerColor(2);

  hPvP2Hi = new TH2F("hPvP2Hi","Recoil (E>=50 MeV) vs Decay Photon2 Theta",180,0,180,180,0,180);
  hPvP2Hi->GetXaxis()->SetTitle("Photon Theta (deg)");
  hPvP2Hi->GetYaxis()->SetTitle("Recoil Theta (deg)");

  Int_t iMassCent = ((TMath::Nint(fRecoMass/50.))*50);

  hMiM = new TH2F("hMissM","Missing Mass vs Recoil K",100,0,500,60,iMassCent-150,iMassCent+150);
  hMiM->GetXaxis()->SetTitle("Recoil Kinetic Energy (MeV)");
  hMiM->GetYaxis()->SetTitle("Missing Mass (MeV)");

};

Pi0PhotGen::~Pi0PhotGen(){

  cout << ("Deleting generator for "+sProcName) << endl;

  delete hPvP1;
  delete hPvP1Lo;
  delete hPvP1Hi;
  delete hPvP2;
  delete hPvP2Lo;
  delete hPvP2Hi;
  delete hMiM;

};

void Pi0PhotGen::Init(Float_t fPolG, Float_t fPolT){

  // Initialization for Pi0 Photoproduction process

  if(bCoher) InitCoher();

  else{
    cout << endl;
    cout << "--------------------------------------------------" << endl << endl;
    cout << "Loading Pi0 Photoproduction data files" << endl;
    
    if(sProcName.Contains("MAID")) InitBase("Pi0_maid_cm_",5);
    else if(sProcName.Contains("SAID")) InitBase("Pi0_said_cm_",5);
    else{
      cout << "Unknown Pi0 files" << endl;
      gSystem->Exit(0);
    }

    cout << "E (MeV)\t\tMax CS (ub)\tTot CS (ub)" << endl;

    if(bIncoh || fPolT==0) CrossUnp();
    else CrossPol(fPolG, fPolT);
  }

  cout << endl;
  cout << "--------------------------------------------------" << endl << endl;
  cout << "Running" << endl << endl;

};

void Pi0PhotGen::InitCoher(){

  // Initialization for Coherent Pi0 Photoproduction process

  cout << endl;
  cout << "--------------------------------------------------" << endl << endl;
  cout << "Constructing Coherent cross sections table" << endl;

  Int_t iEnrN = 0, iAngN = 0;

  Float_t fCPeak, fKpi, fFFSlope, fQ2, fFq2, fCMax = 0, fScale, fRatio = 1;

  if(((iBeamHi-iBeamLo)%5) == 0) iNEnr = (1+((iBeamHi-iBeamLo)/5));
  else if(((iBeamHi-iBeamLo)%4) == 0) iNEnr = (1+((iBeamHi-iBeamLo)/4));
  else if(((iBeamHi-iBeamLo)%3) == 0) iNEnr = (1+((iBeamHi-iBeamLo)/3));
  else if(((iBeamHi-iBeamLo)%2) == 0) iNEnr = (1+((iBeamHi-iBeamLo)/2));
  else iNEnr = (1+(iBeamHi-iBeamLo));
  iNAng = 37;
  iNPhi = iNAng;

  const Int_t iNEnrC = iNEnr;
  const Int_t iNAngC = iNAng;

  // Determine step sizes for energy and angle
  
  iBeamSt = ((iBeamHi-iBeamLo)/(iNEnr-1));
  iAngSt = (180/(iNAng-1));

  // Create arrays to hold input

  fPar = new Float_t**[2];
  fPar[0] = new Float_t*[iNEnrC];
  fPar[1] = new Float_t*[iNEnrC];

  fEnr = new Float_t[iNEnrC];

  for(iEnrN=0; iEnrN<iNEnr; iEnrN++){
    fEnr[iEnrN] = (iBeamLo+(iEnrN*iBeamSt));

    fPar[0][iEnrN] = new Float_t[iNAngC];
    fPar[1][iEnrN] = new Float_t[iNAngC];

    fCPeak = (70.*(290.-fEnr[iEnrN])*(350.-fEnr[iEnrN])/(290.-200.)/(350.-200.)+180.*(fEnr[iEnrN]-200.)*(350.-fEnr[iEnrN])/(290.-200.)/(350.-290.)+200.*(fEnr[iEnrN]-200.)*(fEnr[iEnrN]-290.)/(350.-200.)/(350.-290.));
    fKpi = sqrt(Sqr(fEnr[iEnrN])-Sqr(kMPI0_MEV));
    fFFSlope = (1.05*(268.-fEnr[iEnrN])/(268.-223.)+1.2*(fEnr[iEnrN]-223.)/(268.-223.));
    fCMax = 0;

    for(iAngN=0; iAngN<iNAng; iAngN++){
      fPar[0][iEnrN][iAngN] = (iAngN*iAngSt);
      fQ2 = ((Sqr(fEnr[iEnrN])+Sqr(fKpi)-(2*fEnr[iEnrN]*fKpi*cos(fPar[0][iEnrN][iAngN]*kD2R)))/Sqr(197.));
      fFq2 = exp(-fQ2*fFFSlope);
      fPar[1][iEnrN][iAngN] = Sqr(fFq2*sin(fPar[0][iEnrN][iAngN]*kD2R));

      if(fPar[1][iEnrN][iAngN] > fCMax) fCMax = fPar[1][iEnrN][iAngN];
    }

    fScale = (fCPeak/fCMax);

    if(sTargName.Contains("3He")) fRatio = (3./12.);
    else if(sTargName.Contains("4He")) fRatio = (4./12.);
    else if(sTargName.Contains("16O")) fRatio = (16./12.);

    fFFSlope = ((1.05*(268.-fEnr[iEnrN])/(268.-223.)+1.2*(fEnr[iEnrN]-223.)/(268.-223.))*TMath::Power(fRatio,(1./3.)));

    for(iAngN=0; iAngN<iNAng; iAngN++){
      fQ2 = ((Sqr(fEnr[iEnrN])+Sqr(fKpi)-(2*fEnr[iEnrN]*fKpi*cos(fPar[0][iEnrN][iAngN]*kD2R)))/Sqr(197.));
      fFq2 = exp(-fQ2*fFFSlope);
      fPar[1][iEnrN][iAngN] = (Sqr(fFq2*sin(fPar[0][iEnrN][iAngN]*kD2R)*fRatio)*fScale);
    }

  }

  cout << endl << "Produced " << iNEnr << " energies with " << iNAng << " angle entries" << endl << "Giving " << iBeamSt << " MeV and " << iAngSt << " deg steps" << endl << endl;

  // Reset bin settings for the cross section histograms

  hCrossSec->SetBins(iNEnr,iBeamLo-0.5*iBeamSt,iBeamHi+0.5*iBeamSt,iNAng,-0.5*iAngSt,180+0.5*iAngSt,iNAng,-0.5*iAngSt,180+0.5*iAngSt);
  hCrossMax->SetBins(iNEnr,iBeamLo-0.5*iBeamSt,iBeamHi+0.5*iBeamSt);
  hCrossTot->SetBins(iNEnr,iBeamLo-0.5*iBeamSt,iBeamHi+0.5*iBeamSt);

  cout << "E (MeV)\t\tMax CS (ub)\tTot CS (ub)" << endl;

  CrossUnp();

};

void Pi0PhotGen::CrossPol(Float_t fPolG, Float_t fPolT){

  Int_t iEnrN = 0, iAngN = 0, iPhiN = 0;

  Float_t fSolAng = 0, fSolAngTot = 0, fCVal = 0, fCMax = 0, fCTot = 0;
  Float_t fCVll, fCVlu, fCVul, fCVuu, fCVav;

  // Construct cross sections and fill the corresponding histogram

  for(iEnrN=0; iEnrN<iNEnr; iEnrN++){
    fSolAngTot = 0;
    fCMax = 0;
    fCTot = 0;
    for(iAngN=0; iAngN<iNAng; iAngN++){
      if(iAngN < (iNAng-1)) fSolAng = ((cos(fPar[0][iEnrN][iAngN]*kD2R)-cos(fPar[0][iEnrN][iAngN+1]*kD2R))*kD2R);

      for(iPhiN=0; iPhiN<iNPhi; iPhiN++){
	fCVal = (fPar[1][iEnrN][iAngN]*(1+fPolT*fPar[2][iEnrN][iAngN]*sin(-iPhiN*iAngSt*kD2R)+fPolG*fPolT*fPar[4][iEnrN][iAngN]*cos(-iPhiN*iAngSt*kD2R)));
	hCrossSec->Fill(fEnr[iEnrN], fPar[0][iEnrN][iAngN], (iPhiN*iAngSt), fCVal);
	
	// Check for maximum cross section value
	if(fCVal > fCMax) fCMax = fCVal;

	// Compute average cross section over angular bin
	if((iAngN < (iNAng-1)) && (iPhiN < (iNPhi-1))){
	  fCVll = (fPar[1][iEnrN][iAngN]*(1+fPolT*fPar[2][iEnrN][iAngN]*sin(-iPhiN*iAngSt*kD2R)+fPolG*fPolT*fPar[4][iEnrN][iAngN]*cos(-iPhiN*iAngSt*kD2R)));
	  fCVlu = (fPar[1][iEnrN][iAngN]*(1+fPolT*fPar[2][iEnrN][iAngN]*sin(-(iPhiN+1)*iAngSt*kD2R)+fPolG*fPolT*fPar[4][iEnrN][iAngN]*cos(-(iPhiN+1)*iAngSt*kD2R)));
	  fCVul = (fPar[1][iEnrN][iAngN+1]*(1+fPolT*fPar[2][iEnrN][iAngN+1]*sin(-iPhiN*iAngSt*kD2R)+fPolG*fPolT*fPar[4][iEnrN][iAngN+1]*cos(-iPhiN*iAngSt*kD2R)));
	  fCVuu = (fPar[1][iEnrN][iAngN+1]*(1+fPolT*fPar[2][iEnrN][iAngN+1]*sin(-(iPhiN+1)*iAngSt*kD2R)+fPolG*fPolT*fPar[4][iEnrN][iAngN+1]*cos(-(iPhiN+1)*iAngSt*kD2R)));
	  fCVav = ((fCVll+fCVlu+fCVul+fCVuu)/4.);	  

	  // Sum up total cross section and solid angle (which should be 4pi)
	  fCTot += (2*fCVal*iAngSt*fSolAng);
	  fSolAngTot += (2*iAngSt*fSolAng);
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

Bool_t Pi0PhotGen::NewEvent(Float_t fBeamE){

  // Construct new event

  Reset();
  NewVertex();

  // Set initial state particles

  pPhoton.SetP4Lab(fBeamE,fBeamE,0,0);
  pTarget.SetP4Lab(pTarget.Mass,0,0,0);

  // Introduce a collision to determine energy, momentum, and direction
  // of final states particles

  Collision2B(pPhoton, pTarget, pPi0, pRecoil);

  // Check whether to accept the event or not

  if(Weight(fBeamE,pPi0.ThetaCM,pPi0.PhiCM)){
    return kFALSE;
  }

  // For an incoherent process, use the previously determined directions of the
  // final state particles but determine the proper energies and momenta
  // through a spectral model
  
  if(bIncoh){
    SpecModel(pTarget);
    Collision2B(pPhoton, pTarget, pPi0, pRecoil);    
  }

  // Make sure the event is above threshold

  if(pPi0.Ener<pPi0.Mass){
    return kFALSE;
  }

  // Decay the pi0 into two photons

  Decay2B(pPi0, pDecay1, pDecay2);

  // Fill the particle histograms and construct the Geant4 Ntuple

  pPhoton.HistLab();
  pTarget.HistLab();
  pPi0.HistLab();
  pRecoil.HistLab();
  pDecay1.HistLab();
  pDecay2.HistLab();

  pPhoton.FillNtuple(var,3);
  pRecoil.FillNtuple(var,8);
  pPi0.FillNtuple(var,13);
  pDecay1.FillNtuple(var,18);
  pDecay2.FillNtuple(var,23);
  
  // Fill ntuple

  h1->Fill(var);

  // Fill other test histograms

  hPvP1->Fill(pDecay1.Theta,pRecoil.Theta);
  hPvP2->Fill(pDecay2.Theta,pRecoil.Theta);
  if(pRecoil.KEner<50){
    hPvP1Lo->Fill(pDecay1.Theta,pRecoil.Theta);
    hPvP2Lo->Fill(pDecay2.Theta,pRecoil.Theta);
  }
  else{
    hPvP1Hi->Fill(pDecay1.Theta,pRecoil.Theta);
    hPvP2Hi->Fill(pDecay2.Theta,pRecoil.Theta);
  }

  lvMiss = (pPhoton.P4+pTarget.P4-pDecay1.P4-pDecay2.P4);

  hMiM->Fill(pRecoil.KEner,lvMiss.M());

  return kTRUE;
  
};

void Pi0PhotGen::Reset(){

  // Reset particle information

  pPhoton.Reset();
  pTarget.Reset();
  pPi0.Reset();
  pRecoil.Reset();
  pDecay1.Reset();
  pDecay2.Reset();
  
};

void Pi0PhotGen::SaveHists(TString sFile){

  // Write out test histograms

  cout << ("Saving "+sProcName+" histograms") << endl;

  Int_t i;

  TFile f1(sFile, "RECREATE", "MC_Hists_File");

  pPhoton.WriteHists();
  pTarget.WriteHists();
  pPi0.WriteHists();
  pRecoil.WriteHists();
  pDecay1.WriteHists();
  pDecay2.WriteHists();

  TCanvas *cPvP1 = new TCanvas("cPvP1", "Proton vs Photon Angle", 800, 1000);

  TLine *xlo = new TLine(20,0,20,180);
  TLine *xhi = new TLine(160,0,160,180);
  TLine *ylo = new TLine(0,20,180,20);
  TLine *yhi = new TLine(0,160,180,160);

  Int_t iArCoo[9][4] = {{30,30,10,10},{90,30,90,10},{150,30,170,10},{30,90,10,90},{90,90,90,90},{150,90,170,90},{30,150,10,170},{90,150,90,170},{150,150,170,170}};
  Int_t iPtCoo[9][4] = {{30,30,50,40},{80,30,100,40},{130,30,150,40},{30,85,50,95},{80,85,100,95},{130,85,150,95},{30,140,50,150},{80,140,100,150},{130,140,150,150}};
  TString sPtTxt[9] = {"TAPS,TAPS","CB,TAPS","Out,TAPS","TAPS,CB","CB,CB","Out,CB","TAPS,Out","CB,Out","Out,Out"};

  TArrow *ar[9];
  TPaveText *pt[9];

  for(i=0; i<9; i++){
    ar[i] = new TArrow(iArCoo[i][0],iArCoo[i][1],iArCoo[i][2],iArCoo[i][3],0.01,"|>");
    ar[i]->SetLineWidth(2);
    ar[i]->SetLineColor(4);
    ar[i]->SetFillColor(4);

    pt[i] = new TPaveText(iPtCoo[i][0],iPtCoo[i][1],iPtCoo[i][2],iPtCoo[i][3]);
    pt[i]->AddText(sPtTxt[i]);
  }

  hPvP1Lo->SetStats(kFALSE);
  hPvP1Lo->Draw();
  hPvP1Hi->Draw("same");
  xlo->Draw("same");
  xhi->Draw("same");
  ylo->Draw("same");
  yhi->Draw("same");
  for(i=0; i<9; i++){
    if(i!=4) ar[i]->Draw();
    pt[i]->Draw();
  }

  TCanvas *cPvP2 = new TCanvas("cPvP2", "Proton vs Photon Angle", 800, 1000);

  hPvP2Lo->SetStats(kFALSE);
  hPvP2Lo->Draw();
  hPvP2Hi->Draw("same");
  xlo->Draw("same");
  xhi->Draw("same");
  ylo->Draw("same");
  yhi->Draw("same");
  for(i=0; i<9; i++){
    if(i!=4) ar[i]->Draw();
    pt[i]->Draw();
  }

  hPvP1Lo->Write();
  hPvP1Hi->Write();

  hPvP2Lo->Write();
  hPvP2Hi->Write();

  cPvP1->Write();
  cPvP2->Write();

  hPvP1->Write();
  hPvP2->Write();
  hMiM->Write();

  f1.Close();

};
