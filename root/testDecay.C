{
  gROOT->Reset();
  gSystem->Load("./libs/libGemini.so");

  CNucleus CN(49,110); //constructor
  float fEx = 100;
  float fJ = 30.;
  CN.setCompoundNucleus(fEx,fJ); //specify the excitation energy and spin

  CN.setVelocityCartesian(); // set initial CN velocity to zero
  CAngle spin(CNucleus::pi/2,(float)0.);
  CN.setSpinAxis(spin); //set the direction of the CN spin vector

  for (int i=0;i<2;i++)
    {
     cout << "event = " << i << endl;
     CN.decay(); //decay the compound nucleus

     if (CN.abortEvent)
       {
	 cout << "gemini was not able to decay this event" << endl;
         CN.reset();
         continue;
       }

     // print of number of stable 
     cout << "number of products= " <<CN.getNumberOfProducts() << endl; 
                                               // fragments produced in decay

     CNucleus * products = CN.getProducts(0); //set pointer to firt
                                                    //stable product
     for(;;)
       {
         cout << products->iZ << " " << products->iA << endl;
         products = CN.getProducts();  // go to the next product
         if (products == NULL) break;
       }

   CN.reset();
    }

  CN.setNewIsotope(70,160,50.,10.); // totally different nucleus
  CNucleus CN2(40,75,30.,5.); //define another nucleus

  CNucleus::setTimeTransient((float)5.); //set fission delay
  CLevelDensity::setLittleA(7.3);
  CN.decay();   // let them both decay
  CN2.decay();

  cout << "write out producst from both nuclei " << endl;
  // write out stable products from both decays
  CNucleus *products = CN.getProducts(0);
  for(;;)
     {
       cout << products->iZ << " " << products->iA << endl;
       products = CN.getProducts();  // go to the next product
       if (products == NULL) break;
     }

}
