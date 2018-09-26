// -*- mode: c++ -*- 
//
#ifndef nuclide_
#define nuclide_

#include <string>
#include "CMass.h"
#include "CRandom.h"
#include <sstream>
//*****ROOT**********
//#include "TObject.h"

/**
 *!\brief info on each nucleus
 *
 * basic class CNuclide - stores the basic properties of the nucleus
 */

class CNuclide //: public TObject
{
protected:
    string strChemName; //!< gives isotopes and chemical name, e.g. 208Pb
    string strName; //!< identifation name 
    static const char * name[111]; //!< array containing name of all elements

    
    
public:
    // mod-TU static CMass mass; //!< mass excess class
    CMass *mass; //!< mass excess class   
    CRandom *ran;   //!< random number generator class 
    int iZ; //!< proton number
    int iN; //!< neutron number
    int iA; //!< mass number
    float fJ; //!< spin [hbar]
    float fExpMass; //!< mass excess [MeV]
    float fEx; //!< excitation energy [MeV]
    CNuclide(int iZ ,int iA);  
    CNuclide(int,int,string);
    void init(int,int);
    // mod-TU CNuclide(){};
    CNuclide(); // mod-TU
    float getExcessMass(); 
    const char* getSymbol();
    string getName();
    
    //*****ROOT*******
    //ClassDef(CNuclide,1) //Gemini Nuclide
};

#endif
