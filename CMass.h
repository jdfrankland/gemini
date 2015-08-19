#ifndef mass_
#define mass_
#include "CChart.h"
#include <iomanip>

using namespace std;

/**
 *!\brief mass excesses, pairing energy, etc
 *
 * Class associated with returning quanties associated with 
 * the mass formula
 */

class CMass
{
protected:
    CMass();                  //!< constructor
    static  CMass* fInstance; //!< instance member to make tis a singleton

    float * fExpMass;  //!<experimental mass array
    float * fCalMass;  //!<experimental mass array
    float * fFRM;      //!<finite range mass array
    float * fPair; //!< pairing correction
    //float * fBeta2;    //!<deformation array
    float * fShell;    //!< shell correction
    float * fShell2;    //!< shell correction
    
    void ReadFRDMFile();
    void ReadThomasFermiFile();
    
public:
    // mod-TU CMass();
    ~CMass();
    CChart *chart; //!< contains the considered region of the chart of nuclides
    static CMass* instance(); //!< instance member to make this a singleton

    float getExpMass(int iZ,int iA);
    float getCalMass(int iZ,int iA);
    float getShellCorrection(int iZ, int iA);
    float getShellCorrection2(int iZ, int iA);
    float getFiniteRangeMass(int iZ,int iA);
    float getFiniteRangeMass(float fZ,float fA);
    float getLiquidDropMass(int iZ ,int iA);
    float getPairing(int iZ,int iA);
    float getPairing2(int iZ,int iA);
    
};
#endif
