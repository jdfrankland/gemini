#include "CRun.h"

/**
 * this constructor is everything at the moment.
/param iZcn - CN proton number
/param iAcn is the mass number
/param fEx is the excitation energy in MeV
/param l0 is the critical spin , specifying the max spin in fusion (unit hbar)
/param d0 is the diffuseness of the spin distribution (units hbar)
/param lmax is the maximum spin value considered. (<l0 for just er info)
/param plb is pi-lambda-bar squared in mb
/param numTot is number of Monte Carlo simulations
/param title0 is name of output root file (without ".root" extension)

 */
CRun::CRun(int iZcn, int iAcn, float fEx, float l0, float d0, int lmax, float plb,
	   int numTot,string title0,float vcm/*=0.*/, float thetaDetMin/*=0.*/,
           float thetaDetMax/*=360*/)
{

  cout << title0 << endl;
  string title = title0 + ".root";
  bool residue;
  bool residueDet;

  // calculate the CN spin distribution
  float prob[lmax+1];
  float sum = 0.;
  for (int l=0;l<=lmax;l++)
    {
      prob[l] =  (float)(2*l+1);

      //select either a fermi distribution of the erfc distribution
      //if (d0 > 0.) prob[l] /= (1.+exp(((float)l-l0)/d0));
      if (d0 > 0.) prob[l] *= erfc(((float)l-l0)/d0/4.*sqrt(3.14159))/2.;
      else if ( l > l0) prob[l] = 0.;
      sum += prob[l];
    }
  //normalise to unity
  for (int l=0;l<=lmax;l++)
    {
      prob[l] /= sum;
      if (l > 0) prob[l] += prob[l-1];
    }

  //root file for output
  TFile *f = new TFile(title.c_str(),"RECREATE");

  //crerate compound nucleus
  CNucleus CN(iZcn,iAcn);

  CNucleus.SetEvapMode(1); // turn on Hauser Feshbach
  //CNucleus::setUserGDR();
  //set GEMINI++ parameters
  //CNucleus::setSolution(1);
  //CNucleus::setFissionScaleFactor(7.38);
  //CNucleus::setAddToFisBarrier(7.);
  //CNucleus::setNoIMF();
  //CNucleus::setAddToFisBarrier(4.);
  //CNucleus::setLestone();
  //CLevelDensity::setAfAn(1.036);
  //CLevelDensity::setAimfAn(1.0);
  //CNucleus::setTimeTransient(1.);
  //CTlBarDist::setBarWidth(1.);
  //CTlBarDist::setBarWidth(0.);
  //CYrast::forceSierk();


  CN.setVelocityCartesian((float)0.,(float)0.,(float)0.);
  CAngle spin((float)0.,(float)0.);
  CN.setSpinAxis(spin);

  CN.printParameters();
  
  float asy[20]={0.};
  float asyMultPre[20] = {0.};
  float asyMultPost[20] = {0.};
  float asyMultTot[20] = {0.};

  float Nres = 0.;
  float NresDet = 0.;
  float sumAres = 0.;
  float sumAresDet = 0.;
  float Nfiss = 0.;
  float NpreSad = 0.;
  float NpreScis = 0.;
  float Npost = 0.;
  float Nalpha = 0.;
  float Nproton = 0.;  
  float Nneutron = 0.;
  float NLi6 = 0.;
  float NLi7 = 0.;
  float NBe7 = 0.;
  float Mfis = 0;
  float M2fis = 0.;
  float M0fis = 0.;
  double numberA = 0.;
  double averageA = 0.;

  TH1F histEgammaTot("EgammaTot","",100,0,50);
  histEgammaTot.GetXaxis()->SetTitle("E_{#gamma} [MeV]");
  histEgammaTot.GetYaxis()->SetTitle("#sigma(E_{#gamma}) [mb/MeV]");
  histEgammaTot.SetTitle("distribution of total gamma energy for all events");


  TH1F histEgamma("Egamma","",100,0,50);
  histEgamma.GetXaxis()->SetTitle("E_{#gamma} [MeV]");
  histEgamma.GetYaxis()->SetTitle("#sigma(E_{#gamma}) [mb/MeV]");
  histEgamma.SetTitle("distribution of gamma energies for all events");

  TH1F histEgammaER("EgammaER","",100,0,50);
  histEgammaER.GetXaxis()->SetTitle("E_{#gamma} [MeV] for residues");
  histEgammaER.GetYaxis()->SetTitle("#sigma(E_{#gamma}) [mb/MeV]");
  histEgammaER.SetTitle("distribution of gamma energies for evap. residues");

  TH1F histER("histER","",91,-0.5,90.5);
  histER.GetXaxis()->SetTitle("J_{CN} [hbar]");
  histER.GetYaxis()->SetTitle("#sigma(J) [mb]");
  histER.SetTitle("spin distribution of events that form residues");

  TH1F histERxn("histERxn","",91,-0.5,90.5);
  histERxn.GetXaxis()->SetTitle("J_{CN} [hbar]");
  histERxn.GetYaxis()->SetTitle("#sigma(J) [mb]");
  histERxn.SetTitle("spin distribution of residue that decay by neutrons only");
  TH1F histFis("histFis","",91,-0.5,90.5);
  histFis.GetXaxis()->SetTitle("J_{CN} [hbar]");
  histFis.GetYaxis()->SetTitle("#sigma(J) [mb]");
  histFis.SetTitle("spin distribution of events that fission");

  TH1F histFus("histFus","",91,-0.5,90.5);
  histFus.GetXaxis()->SetTitle("J [hbar]");
  histFus.GetYaxis()->SetTitle("#sigma(J) [mb]");
  histFus.SetTitle("spin distribution of all fusion events");

  TH1F histA("histA","",231,-0.5,230.5);
  histA.GetXaxis()->SetTitle("A");
  histA.GetYaxis()->SetTitle("#sigma(A) [mb]");
  histA.SetTitle(" inclusive mass distribution");


  TH1F histZ("histZ","",93,-0.5,92.5);
  histZ.GetXaxis()->SetTitle("Z");
  histZ.GetYaxis()->SetTitle("#sigma(Z) [mb]");
  histZ.SetTitle(" inclusive charge distribution");


  TH1F histZ_fis("histZ_fis","",93,-0.5,92.5);
  histZ_fis.GetXaxis()->SetTitle("Z");
  histZ_fis.GetYaxis()->SetTitle("#sigma(Z) [mb]");
  histZ_fis.SetTitle("charge distribution for fission events");

  TH1F histZ_nofis("histZ_nofis","",93,-0.5,92.5);
  histZ_nofis.GetXaxis()->SetTitle("Z");
  histZ_nofis.GetYaxis()->SetTitle("#sigma(Z) [mb]");
  histZ_nofis.SetTitle("charge distribution for non-fission events");

  TH1F histN("histN","",153,-0.5,152.5);
  histN.GetXaxis()->SetTitle("N");
  histN.GetYaxis()->SetTitle("#sigma(N) [mb]");
  histN.SetTitle(" inclusive neutron-number distribution");

  //the following are differential multiplicity 
  TH1F keNeutron("keNeutron","",50,0,50);
  keNeutron.GetXaxis()->SetTitle("E_{k} [MeV]");
  keNeutron.GetYaxis()->SetTitle("dm/dE [MeV^{-1}]");
  keNeutron.SetTitle("neutron energy spectra in coincidence with residues");

  TH1F keAlpha("keAlpha","",50,0,50);
  keAlpha.GetXaxis()->SetTitle("E_{k} [MeV]");
  keAlpha.GetYaxis()->SetTitle("dm/dE [MeV^{-1}]");
  keAlpha.SetTitle("#alpha particle energy spectra in coincidence with residues");


  TH1F keProton("keProton","",50,0,50);
  keProton.GetXaxis()->SetTitle("E_{k} [MeV]");
  keProton.GetYaxis()->SetTitle("dm/dE [MeV^{-1}]");
  keProton.SetTitle("proton energy spectra in coincidence with residues");


  TH1F keLi6("keLi6","",60,0,60);
  keLi6.GetXaxis()->SetTitle("E_{k} [MeV]");
  keLi6.GetYaxis()->SetTitle("dm/dE [MeV^{_1}]");
  keLi6.SetTitle("^{6}Li energy spectra in coincidence with residues");


  TH1F keLi7("keLi7","",60,0,60);
  keLi7.GetXaxis()->SetTitle("E_{k} [MeV]");
  keLi7.GetYaxis()->SetTitle("dm/dE [MeV^{-1}]");
  keLi7.SetTitle("^{7}Li energy spectra in coincidence with residues");


  TH1F keBe7("keBe7","",60,0,60);
  keBe7.GetXaxis()->SetTitle("E_{k} [MeV]");
  keBe7.GetYaxis()->SetTitle("dm/dE [MeV^{-1}]");
  keBe7.SetTitle("^{7}Be energy spectra in coincidence with residues");



  //the following are differential cross sections 
  TH1F sigNeutron("sigNeutron","",50,0,50);
  sigNeutron.GetXaxis()->SetTitle("E_{k} [MeV]");
  sigNeutron.GetYaxis()->SetTitle("d#sigma/dE [mb/MeV]");
  sigNeutron.SetTitle("inclusive neutron energy spectra");

  TH1F sigAlpha("sigAlpha","",50,0,50);
  sigAlpha.GetXaxis()->SetTitle("E_{k} [MeV]");
  sigAlpha.GetYaxis()->SetTitle("d#sigma/dE [mb/MeV]");
  sigAlpha.SetTitle("inclusive #alpha-particle energy spectra");


  TH1F sigProton("sigProton","",50,0,50);
  sigProton.GetXaxis()->SetTitle("E_{k} [MeV]");
  sigProton.GetYaxis()->SetTitle("d#sigma/dE [mb/MeV]");
  sigProton.SetTitle("inclusive proton energy spectra");

  TH2F histAL("histAL","",251,-0.5,250.5,100,0,100);
  histAL.GetXaxis()->SetTitle("A");
  histAL.GetYaxis()->SetTitle("J_{CN} [hbar]");
  histAL.SetTitle("inclusive mass and CN spin distribution");


  TH2F histZN("histZN","",152,-0.5,151.5,152,-0.5,151.5);
  histZN.GetXaxis()->SetTitle("N");
  histZN.GetYaxis()->SetTitle("Z");
  histZN.SetTitle("inclusive joint N_Z distributions of all fragments");



  TH1F keFF("keFF","",150,0,150);
  keFF.GetXaxis()->SetTitle("E_{k} [MeV]");
  keFF.GetYaxis()->SetTitle("d#sigma/dE [mb/MeV]");
  keFF.SetTitle("Fission Fragment kinetic-energy spectrum");


  TH1F velFF("velFF","",100,0,4.);
  velFF.GetXaxis()->SetTitle("v [cm/ns]");
  velFF.GetYaxis()->SetTitle("d#sigma/dv [mb/cm/ns]");
  velFF.SetTitle("Fission Fragment velocity spectrum");



  TH1F kePreSad("kePreSad","",100,0,30);
  kePreSad.GetXaxis()->SetTitle("E_{k} [MeV]");
  kePreSad.GetYaxis()->SetTitle("dm/dE [MeV^{-1}]");
  kePreSad.SetTitle("pre-saddle neutron multiplicity spectrum");

  TH1F kePreSS("keSS","",100,0,30);
  kePreSS.GetXaxis()->SetTitle("E_{k} [MeV]");
  kePreSS.GetYaxis()->SetTitle("dm/dE [MeV^{-1}]");
  kePreSS.SetTitle("saddle-to-scission neutron multiplicity spectrum");


  TH1F kePreSc("kePreSc","",100,0,30);
  kePreSc.GetXaxis()->SetTitle("E_{k} [MeV]");
  kePreSc.GetYaxis()->SetTitle("dm/dE [MeV^{-1}]");
  kePreSc.SetTitle("pre-scission neutron multiplicity spectrum");

  TH1F kePost("kePost","",100,0,30);
  kePost.GetXaxis()->SetTitle("E_{k} [MeV]");
  kePost.GetYaxis()->SetTitle("dm/dE [MeV^{-1}]");
  kePost.SetTitle("post scission neutron spectrum");


  bool f14=1;
  bool f12=1;
  bool f34=1;
  for (int i=0;i<numTot;i++)
    {
      float weight = 1.;
      //cout <<"event= " <<  i << endl;
      if (i > numTot*.25 && f14)
	{
	  cout << " 25%" << flush;
          f14 = 0;
	} 
      if (i > numTot*.5 && f12)
	{
	  cout << '\xd' << " 50%" << flush;
          f12 = 0;
	} 
      if (i > numTot*.75 && f34)
	{
	  cout << '\xd'<< " 75%" << endl;
          f34 = 0;
	} 
      float ran = CN.ran->Rndm();
      int l = 0;
      for (;;)
	{
	  if (ran < prob[l]) break;
          l++;
	}
      // l = 43; //rjc
      //fEx = 77.83;

      CN.setCompoundNucleus(fEx,(float)l);      

	{
	  //if (i%100==0)cout << "l= "<< l  << " i= " << i << endl;
	 CN.setWeightIMF(); //turn on enhanced IMF emission
         CN.decay();

        if (CN.abortEvent)
          {
	    CN.reset();
            continue;
          }

	histEgammaTot.Fill(CN.getSumGammaEnergy());
        for(int i =0; i<CN.getnGammaRays(); i++)
	  {
	    histEgamma.Fill(CN.getGammaRayEnergy(i));
	    if (CN.isResidue())histEgammaER.Fill(CN.getGammaRayEnergy(i));
	  }

         int iStable = CN.getNumberOfProducts();

         CNucleus *productER = CN.getProducts(iStable-1);
	 weight *= productER->getWeightFactor();

	 if(productER->iZ == iZcn)
	   {
	     averageA += (double)productER->iA*(double)weight;
             numberA += (double)weight;
	   }

         int iZres = productER->iZ;
         //int iAres = productER->iA;
         int multTot = 0;
         int iZ, iA;

         CNucleus * product = CN.getProducts(0);
         histFus.Fill(l,weight);

         if (CN.isResidue()) 
	   {
             float * vv;
             vv = productER->getVelocityVector();
             vv[2] += vcm;
             float vvv = sqrt(pow(vv[0],2)+pow(vv[1],2)+pow(vv[2],2));
             float AngleDeg = acos(vv[2]/vvv)*180./3.14159;
             if (AngleDeg > thetaDetMin && AngleDeg < thetaDetMax) 
	       residueDet = 1;
             else residueDet = 0;
	     residue = 1;
             histER.Fill(l,weight);
	     if (iZres == iZcn)
	         {
                    histERxn.Fill(l,weight);
	          }

             Nres += weight;
             sumAres += weight*(float)productER->iA;
             if (residueDet) 
	       {
               NresDet += weight;
               sumAresDet += weight*(float)productER->iA;
	       }
	   }
	 else
	   {
	     residue = 0;
             residueDet = 0;
	   }

	 if (CN.isSymmetricFission())
	   {
	     Nfiss += weight;
             histFis.Fill(l,weight);
	   }


         int iAmax = 0;
         int iAnext = 0;
         float vmax = 0.;
         float vnext = 0.;
         float emax = 0.;
         float enext = 0.;
         for (int j=0;j<iStable;j++)
	   {
             iZ = product->iZ;
             iA = product->iA;

             if (iA > iAmax) 
	       {
                 iAnext = iAmax;
                 vnext = vmax;
                 enext = emax;
                 iAmax = iA;
		 emax = product->getKE();
		 vmax = product->getVelocity();
	       }
	     else if (iA > iAnext)
	       {
                 iAnext = iA;
                 enext = product->getKE();
                 vnext = product->getVelocity();
	       }
             //cout << iZ << " " << iA << endl;
             if (product->getTime() < 0.) 
	       {
		 cout << "negative time" << endl;
                 cout << iZ << " " << iA << " " << 
		   product->getParent()->iZ << " " << 
		   product->getParent()->iA << " " 
		      << product->getParent()->fEx << " " 
		      << product->getParent()->fJ << endl;
	       }

             histZ.Fill(iZ,weight);

             if (CN.isSymmetricFission())histZ_fis.Fill(iZ,weight);
             else histZ_nofis.Fill(iZ,weight);
             histA.Fill(iA,weight);

             histAL.Fill(iA,l,weight);
             histN.Fill(iA-iZ,weight);
             histZN.Fill(iA-iZ,iZ,weight);
             if (iZ == 0 && iA == 1) 
	       {
		 if (residueDet) //iARes >= Ares)
		   {
   		   keNeutron.Fill(product->getKE(),weight);
                   Nneutron += weight;
		   }
   		   sigNeutron.Fill(product->getKE(),weight);
                 multTot++;
                 if (CN.isSymmetricFission())
		   {
                    if (product->origin == 0)
		      { 
		      kePreSad.Fill(product->getKE(),weight);
		      NpreSad += weight;
		      }
                    if (product->origin == 1) 
                      kePreSS.Fill(product->getKE(),weight);
                    if (product->origin <= 1) 
		      {
			NpreScis += weight;
                      kePreSc.Fill(product->getKE(),weight);
		      }
		    if (product->origin > 1) 
		      {
			Npost += weight;
                      kePost.Fill(product->getKE(),weight);
		      }
		   }
	       }
	     else if (iZ == 1 && iA == 1)
	       {
	        if(residueDet) 
	          {
		   keProton.Fill(product->getKE(),weight);
                   Nproton += weight;
	          }
		   sigProton.Fill(product->getKE(),weight);
	       }
	     else if (iZ == 2 && iA == 4)
	       {
		 if(residueDet) 
		   {
	             keAlpha.Fill(product->getKE(),weight);
                     Nalpha += weight;
		   }
	          sigAlpha.Fill(product->getKE(),weight);
		 }
	     else if (iZ == 3 && iA == 6)
	       {
		 if(residueDet) 
		   {

	             keLi6.Fill(product->getKE(),weight);
                     NLi6 += weight;
		   }
		 }
	     else if (iZ == 3 && iA == 7)
	       {
		 if(residueDet) 
		   {

	             keLi7.Fill(product->getKE(),weight);
                     NLi7 += weight;
		   }
		 }
	     else if (iZ == 4 && iA == 7)
	       {
		 if(residueDet) 
		   {

	             keBe7.Fill(product->getKE(),weight);
                     NBe7 += weight;
		   }
		 }

             if (iZ > 5 && CN.isSymmetricFission()) 
	       {
                 keFF.Fill(product->getKE(),weight);
                 velFF.Fill(product->getVelocity(),weight);
                 Mfis += (float)product->iA;
                 M2fis += pow((float)product->iA,2);
                 M0fis += 1.;
	       }
	     product=CN.getProducts();
	    }


         float Amax = (float)iAmax/(float)(iAmax+iAnext)*162.;
         float Anext = (float)iAnext/(float)(iAmax+iAnext)*162.;
	 int iasy = (int)(Amax/10);
         asy[iasy] += weight;
         asyMultPre[iasy] += weight*(float)CN.getMultPre();
         asyMultPost[iasy] += weight*(float)CN.getMultPost();
         asyMultTot[iasy] += weight*(float)multTot;
	 iasy = (int)(Anext/10);
         asy[iasy] += weight;
         asyMultPre[iasy] += weight*(float)CN.getMultPre();
         asyMultPost[iasy] += weight*(float)CN.getMultPost();
         asyMultTot[iasy] += weight*(float)multTot;
   

        CN.reset();
	}
      }

  title = title0+"M.dat";
  ofstream ofFile(title.c_str());
  for (int i=0;i<20;i++)
    {
      if (asy[i] == 0) continue;
      ofFile << i*10 + 5 << " " << asyMultPre[i]/asy[i] << " " <<
	asyMultPost[i]/asy[i] << " " << asyMultTot[i]/asy[i] << endl;
    }



  histA.Scale(plb/(float)numTot*sum);
  histZ.Scale(plb/(float)numTot*sum);
  histZ_fis.Scale(plb/(float)numTot*sum);
  histZ_nofis.Scale(plb/(float)numTot*sum);
  histN.Scale(plb/(float)numTot*sum);
  histZN.Scale(plb/(float)numTot*sum);
  velFF.Scale(plb/(float)numTot*sum/velFF.GetBinWidth(1));
  keFF.Scale(plb/(float)numTot*sum/keFF.GetBinWidth(1));



  //for multiplicity spectrra in coincidence with fission

  kePreSad.Scale(1./Nfiss/kePreSad.GetBinWidth(1));
  kePreSS.Scale(1./Nfiss/kePreSS.GetBinWidth(1));
  kePreSc.Scale(1./Nfiss/kePreSc.GetBinWidth(1));
  kePost.Scale(1./Nfiss/kePost.GetBinWidth(1));




  //for multiplicity spectra in coincidence with residues
  keAlpha.Scale(1./NresDet/keAlpha.GetBinWidth(1));
  keProton.Scale(1./NresDet/keProton.GetBinWidth(1));
  keNeutron.Scale(1./NresDet/keNeutron.GetBinWidth(1));
  keLi6.Scale(1./NresDet/keLi6.GetBinWidth(1));
  keLi7.Scale(1./NresDet/keLi7.GetBinWidth(1));
  keBe7.Scale(1./NresDet/keBe7.GetBinWidth(1));
  sigNeutron.Scale(plb/(float)numTot*sum/sigNeutron.GetBinWidth(1));
  sigProton.Scale(plb/(float)numTot*sum/sigProton.GetBinWidth(1));
  sigAlpha.Scale(plb/(float)numTot*sum/sigAlpha.GetBinWidth(1));

  histEgamma.Scale(plb/(float)numTot*sum/histEgamma.GetBinWidth(1));  
  histEgammaER.Scale(plb/(float)numTot*sum/histEgammaER.GetBinWidth(1));  
  histEgammaTot.Scale(plb/(float)numTot*sum/histEgammaTot.GetBinWidth(1));  

  f->Write();
  cout << "NresDet= " << NresDet << " Nneut= " << Nneutron << " NProt= " <<
    Nproton << " Nalpha= " << Nalpha << " NLi6= " << NLi6 << " NLi7= " << NLi7
       << " NBe7= " << NBe7 << endl;
  cout << "Li6 mult = " << NLi6/NresDet << endl;
  cout << "Li7 mult = " << NLi7/NresDet << endl;
  cout << "Be7 mult = " << NBe7/NresDet << endl;
  cout << "neutron mult= " << Nneutron/NresDet << endl;
  cout << "proton mult= " << Nproton/NresDet << endl;
  cout << "alpha mult= " << Nalpha/NresDet << endl;

  cout << " mean ER A = " << sumAres/Nres << endl;
  cout << " for det res = " << sumAresDet/NresDet << endl;

  float xER = Nres/(float)numTot*sum*plb;
  float xFiss2 = Nfiss/(float)numTot*sum*plb;
  cout << "sigmaER = " << xER << " mb " << endl;

  float xFus = 0.;
  for (int l=0;l<200;l++)
    {
      float xx =  (float)(2*l+1);
      if (d0 > 0.) xx /=(1.+exp(((float)l-l0)/d0));
      else if (l > l0) break;
      xFus += xx;
    }
  xFus *= plb;
  float xFis = xFus - xER;

  cout << "fusion xsec= " << xFus << " mb" << endl;
  cout << "fission xsec= " << xFis << " mb " << xFiss2 <<  " mb " <<   endl;

  if (Nfiss > 0.)
    {
     cout << "preSaddle neut mult = " << NpreSad/Nfiss << endl;
     cout << "preScis neut mult = " << NpreScis/Nfiss << endl;
     cout << "post neut mult = " << Npost/Nfiss << endl;

     float Mav = Mfis/M0fis;
     cout << "mean fission mass = " << Mav << endl;
     float sigma2 = M2fis/(M0fis-1) - M0fis/(M0fis-1)*pow(Mav,2);
     //float sigma = sqrt(sigma2);
     cout << "sigma2M= " << sigma2 << endl;
    }
  if (numberA > 0) cout << "average x for xn products is " << (float)iAcn-
   averageA/numberA << endl;



}
