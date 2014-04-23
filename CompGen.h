// Compton Scattering and Pi0P/PiPN Photoproduction Event Generator
// Designed for use with the MAMI A2 Geant4 Simulation
//
// Compton scattering event generator class
//
// Author - P. Martel

class CompGen : public BaseGen {
 protected:
  BasePart pPhoton, pTarget, pScatter, pRecoil;
  Float_t fPhotEk, fPhotEkCM;
  TLorentzVector lvMiss, lvScat, lvScatCM, lvReco, lvRecoCM;
  TH2F *hPvP, *hPvPLo, *hPvPHi, *hMiM;

 public:
  CompGen(TString, TString, TString, Float_t, Float_t, Int_t, Int_t, Int_t);
  ~CompGen();
  void Init();
  void InitCoher();
  Bool_t NewEvent(Float_t);
  void Reset();
  void SaveHists(TString);
};

CompGen::CompGen(TString sName, TString sTarget, TString sBase, Float_t fTargMass, Float_t fRecoMass, Int_t iRecoG3id, Int_t beamlo, Int_t beamhi) : BaseGen(sName, sTarget, sBase, beamlo, beamhi), pPhoton("Photon",0), pTarget("Target",fTargMass), pScatter("Scatter",0), pRecoil("Recoil",fRecoMass) {

  cout << "Constructing generator" << endl;

  // Setup Geant4 Ntuple

  Int_t ptag[2] = {iRecoG3id,1};

  InitNtuple(2,ptag);

  // Create some histograms for testing
  
  hPvP = new TH2F("hPvP","Recoil vs Photon Theta",180,0,180,180,0,180);
  hPvP->GetXaxis()->SetTitle("Photon Theta (deg)");
  hPvP->GetYaxis()->SetTitle("Recoil Theta (deg)");

  hPvPLo = new TH2F("hPvPLo","Recoil (E<50 MeV) vs Photon Theta",180,0,180,180,0,180);
  hPvPLo->GetXaxis()->SetTitle("Photon Theta (deg)");
  hPvPLo->GetYaxis()->SetTitle("Recoil Theta (deg)");
  hPvPLo->SetMarkerColor(2);

  hPvPHi = new TH2F("hPvPHi","Recoil (E>=50 MeV) vs Photon Theta",180,0,180,180,0,180);
  hPvPHi->GetXaxis()->SetTitle("Photon Theta (deg)");
  hPvPHi->GetYaxis()->SetTitle("Recoil Theta (deg)");

  Int_t iMassCent = ((TMath::Nint(fRecoMass/50.))*50);

  hMiM = new TH2F("hMissM","Missing Mass vs Recoil K",100,0,500,60,iMassCent-150,iMassCent+150);
  hMiM->GetXaxis()->SetTitle("Recoil Kinetic Energy (MeV)");
  hMiM->GetYaxis()->SetTitle("Missing Mass (MeV)");
};

CompGen::~CompGen(){

  cout << "Deleting generator" << endl;

  delete hPvP;
  delete hPvPLo;
  delete hPvPHi;
  delete hMiM;

};

void CompGen::Init(){

  // Initialization for Compton Scattering process

  if(bCoher && !bIsotW) InitCoher();

  else if(bIsotW){
    cout << endl;
    cout << "--------------------------------------------------" << endl << endl;
    cout << "Using isotropic Comtpon distributions" << endl;
  }

  else{
    cout << endl;
    cout << "--------------------------------------------------" << endl << endl;
    cout << "Loading Compton data files" << endl;

    InitBase(sBaseName);

    cout << "E (MeV)\t\tMax CS (nb)\tTot CS (nb)" << endl;

    CrossGen();
  }

  cout << endl;
  cout << "--------------------------------------------------" << endl << endl;
  cout << "Running" << endl << endl;

  if(bSaveT){
    t1->Branch("Phot",&fPhotEk);
    t1->Branch("PhotCM",&fPhotEkCM);
    t1->Branch("Scat",&lvScat);
    t1->Branch("ScatCM",&lvScatCM);
    t1->Branch("Reco",&lvReco);
    t1->Branch("RecoCM",&lvRecoCM);
  }

};

void CompGen::InitCoher(){

  // Initialization for Coherent Compton Scattering process

  cout << endl;
  cout << "--------------------------------------------------" << endl << endl;
  cout << "Constructing Coherent cross sections table" << endl;

  Int_t iParN = 0, iEnrN = 0, iAngN = 0, iNPar = 9;

  Float_t fCPeak, fSigTh, fHel40, fCar40, fRatio = 1;

  if(((iBeamHi-iBeamLo)%5) == 0) iNEnr = (1+((iBeamHi-iBeamLo)/5));
  else if(((iBeamHi-iBeamLo)%4) == 0) iNEnr = (1+((iBeamHi-iBeamLo)/4));
  else if(((iBeamHi-iBeamLo)%3) == 0) iNEnr = (1+((iBeamHi-iBeamLo)/3));
  else if(((iBeamHi-iBeamLo)%2) == 0) iNEnr = (1+((iBeamHi-iBeamLo)/2));
  else iNEnr = (1+(iBeamHi-iBeamLo));
  iNAng = 37;
  iNPhi = 73;
  
  const Int_t iNEnrC = iNEnr;
  const Int_t iNAngC = iNAng;
  const Int_t iNParC = iNPar;

  // Determine step sizes for energy and angle
  
  iBeamSt = ((iBeamHi-iBeamLo)/(iNEnr-1));
  iAngSt = (180/(iNAng-1));

  // Create arrays to hold input

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
    fEnr[iEnrN] = (iBeamLo+(iEnrN*iBeamSt));

    fPar[0][iEnrN] = new Float_t[iNAngC];
    fPar[1][iEnrN] = new Float_t[iNAngC];

    fCPeak = ((800.*(253.-fEnr[iEnrN])*(310.-fEnr[iEnrN])/(253.-206.)/(310.-206.)+1500.*(fEnr[iEnrN]-206.)*(310.-fEnr[iEnrN])/(253.-206.)/(310.-253.)+3000.*(fEnr[iEnrN]-206.)*(fEnr[iEnrN]-253.)/(310.-206.)/(310.-253.))/1000.);
    fSigTh = ((75.-(fEnr[iEnrN]-253.)*20./(310.-206.))/2.15);
    fHel40 = (fCPeak*(TMath::Exp(-(Sqr(40.))/(2*Sqr(fSigTh))))+0.02);
    fCar40 = (1.8*(250.-fEnr[iEnrN])*(300.-fEnr[iEnrN])/(250.-200.)/(300.-200.)+3.25*(fEnr[iEnrN]-200.)*(300.-fEnr[iEnrN])/(250.-200.)/(300.-250.)+3.25*(fEnr[iEnrN]-200.)*(fEnr[iEnrN]-250.)/(300.-200.)/(300.-250.));
    
    if(sTargName.Contains("3He")) fRatio = (3./4.);
    else if(sTargName.Contains("12C")) fRatio = (fCar40/fHel40);
    else if(sTargName.Contains("16O")) fRatio = ((16.*fCar40)/(12.*fHel40));

    for(iAngN=0; iAngN<iNAng; iAngN++){
      fPar[0][iEnrN][iAngN] = (iAngN*iAngSt);
      fPar[1][iEnrN][iAngN] = (((fCPeak*(TMath::Exp(-Sqr(fPar[0][iEnrN][iAngN])/(2*Sqr(fSigTh)))))+0.02)*fRatio);
    }
  }

  cout << endl << "Produced " << iNEnr << " energies with " << iNAng << " angle entries" << endl << "Giving " << iBeamSt << " MeV and " << iAngSt << " deg steps" << endl << endl;

  // Reset bin settings for the cross section histograms

  hCrossSec->SetBins(iNEnr,iBeamLo-0.5*iBeamSt,iBeamHi+0.5*iBeamSt,iNAng,0-0.5*iAngSt,180+0.5*iAngSt,iNPhi,-180-0.5*iAngSt,180+0.5*iAngSt);
  hCrossMax->SetBins(iNEnr,iBeamLo-0.5*iBeamSt,iBeamHi+0.5*iBeamSt);
  hCrossTot->SetBins(iNEnr,iBeamLo-0.5*iBeamSt,iBeamHi+0.5*iBeamSt);

  cout << "E (MeV)\t\tMax CS (ub)\tTot CS (ub)" << endl;

  CrossGen();
  
};

Bool_t CompGen::NewEvent(Float_t fBeamE){

  Bool_t bCheck = kTRUE;

  Double_t dCTot = hCrossTot->Interpolate(fBeamE);
  if(!bIsotW && ((dCTot*dConv) <= (gRandom->Rndm()))) return kFALSE;

  while(bCheck){

    // Construct new event

    Reset();
    NewVertex();

    // Set initial state particles

    pPhoton.SetP4Lab(fBeamE,fBeamE,0,0);
    pTarget.SetP4Lab(pTarget.Mass,0,0,0);

    // Introduce a collision to determine energy, momentum, and direction
    // of final states particles

    Collision2B(pPhoton, pTarget, pScatter, pRecoil);

    // Check whether to reject the event based on the selected weighting

    bCheck = Reject(fBeamE,pScatter.Theta,pScatter.Phi);
  }

  // For an incoherent process, use the previously determined directions of the
  // final state particles but determine the proper energies and momenta
  // through a spectral model
  
  if(bIncoh){
    SpecModel(pTarget);
    Collision2B(pPhoton, pTarget, pScatter, pRecoil);    
  }

  // Fill the particle histograms and construct the Geant4 Ntuple

  pPhoton.HistLab();
  pTarget.HistLab();
  pScatter.HistLab();
  pRecoil.HistLab();

  pPhoton.FillNtuple(var,3);
  pRecoil.FillNtuple(var,8);
  pScatter.FillNtuple(var,13);
  
  // Fill ntuple

  if(bSaveN) h1->Fill(var);
  
  // Fill tree

  fPhotEk = pPhoton.KEner;
  fPhotEkCM = pPhoton.KEnerCM;
  lvScat = pScatter.P4;
  lvScatCM = pScatter.P4CM;
  lvReco = pRecoil.P4;
  lvRecoCM = pRecoil.P4CM;

  if(bSaveT) t1->Fill();

  // Fill other test histograms

  hPvP->Fill(pScatter.Theta,pRecoil.Theta);
  if(pRecoil.KEner<50) hPvPLo->Fill(pScatter.Theta,pRecoil.Theta);
  else hPvPHi->Fill(pScatter.Theta,pRecoil.Theta);

  lvMiss = (pPhoton.P4+pTarget.P4-pScatter.P4);

  hMiM->Fill(pRecoil.KEner,lvMiss.M());

  return kTRUE;
  
};

void CompGen::Reset(){

  // Reset particle information

  pPhoton.Reset();
  pTarget.Reset();
  pScatter.Reset();
  pRecoil.Reset();
  
};

void CompGen::SaveHists(TString sFile){

  if(bSaveH){
    // Write out test histograms
    
    cout << "Saving histograms" << endl;
    
    Int_t i;
    
    TFile f1(sFile, "RECREATE", "MC_Hists_File");
    
    pPhoton.WriteHists();
    pTarget.WriteHists();
    pScatter.WriteHists();
    pRecoil.WriteHists();
    
    hPvPLo->Write();
    hPvPHi->Write();
    
    TCanvas *cPvP = new TCanvas("cPvP", "Proton vs Photon Angle", 800, 1000);
    
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
    
    hPvPLo->SetStats(kFALSE);
    hPvPLo->Draw();
    hPvPHi->Draw("same");
    xlo->Draw("same");
    xhi->Draw("same");
    ylo->Draw("same");
    yhi->Draw("same");
    for(i=0; i<9; i++){
      if(i!=4) ar[i]->Draw();
      pt[i]->Draw();
    }
    
    cPvP->Write();
    
    hPvP->Write();
    hMiM->Write();
    
    hCrossMax->Write();
    hCrossTot->Write();
    
    f1.Close();
  }

};
