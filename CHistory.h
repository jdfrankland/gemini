#ifndef chistory_h
#define chistory_h

#include "CNucleus.h"
#include <map>
#include <limits>
#include <iostream>
#include <string>
#include <algorithm>
#include <sstream>

/**
 *!\brief figures out the history of a fragment
 *
 * Figures out the de-excitation history of a fragment, and condenses this
 * information as either a string of characters or a string of digits (i.e. an
 * integer). Each character/digit represents one or more steps in the
 * de-excitation process. Some of the possible processes may never be realised
 * by a particular de-excitation model, but the encoding is meant to be allow
 * comparisons among different de-excitation models. GEMINI++, for instance,
 * will never produce fragments through process 3 ('m'), i.e.
 * multifragmentation.
 *
 * The currently defined possible character values and their meanings are the
 * following:
 *
 * digit   char   meaning
 *   1      e     evaporation product
 *   2      E     evaporation residue
 *   3      m     multifragmentation
 *   4      a     light partner in asymmetric fission or IMF emission
 *   5      A     heavy partner in asymmetric fission or IMF emission
 *   6      f     light partner in fission
 *   7      F     heavy partner in fission
 *   8      s     saddle-to-scission emission
 *   9      n     non-statistical emission (decay)
 *
 * Consecutive de-excitation steps are sometimes encoded as a single step, to
 * save some space in the string. The rules are the following:
 * 1. if last step == 'E' and new step == 'e', condense the two as 'e';
 * 2. if last step == 'E' and new step == 'E', condense the two as 'E';
 * 3. if last step == 's' and new step == 's', condense the two as 's'.
 *
 * Usage:
 *
 * CNucleus *CN = new CNucleus(iZ, iA, fEx, fJ);
 * CN->decay();
 * CHistory theHistoryClass(CN);
 * int n = CN->GetNumberOfProducts();
 * for(int i=0; i<n; ++i) {
 *   CNucleus *product = CN->getProducts(i);
 *   std::string historyString = theHistoryClass.getHistoryAsString(product);
 *   // or
 *   int32_t historyInt = theHistoryClass.getHistory(product);
 * }
 */

class CHistory {

  typedef std::map<CNucleus *, std::string> HistoryMap;

  public:

  /**
   * destructor
   */
  ~CHistory() {};

  /**
   * returns the history of a fragment as an integer (a sequence of digits)
   * \param p is a pointer to the fragment
   */
  int32_t getHistory(CNucleus *p) {
    std::string hist = getHistoryAsString(p).substr(0, maxInt32Len);
    std::transform(hist.begin(),
        hist.end(),
        hist.begin(),
        theTranslator
        );
    std::stringstream ss;
    ss << hist;
    int32_t iHist;
    ss >> iHist;
    return iHist;
  }

  /**
   * returns the history of a fragment as a string (a sequence of chars)
   * \param p is a pointer to the fragment
   */
  std::string getHistoryAsString(CNucleus *p) {
    HistoryMap::const_iterator iter = theMap.find(p);
    if(iter != theMap.end())
      return iter->second;
    else {
      std::cerr << "Unknown CNucleus pointer in HistoryMap" << std::endl;
      return 0;
    }
  }

  /**
   * constructor
       \param np pointer to compound nucleus
   */
  CHistory(CNucleus *np) {
    evap = CEvap::instance();
    maxEvapZ = evap->maxZ;
    theMap[np]="";
    tagDaughters(np, "");
  }

  static const char Evaporation;
  static const char EvaporationResidue;
  static const char Multifragmentation;
  static const char AsymmetricFissionLight;
  static const char AsymmetricFissionHeavy;
  static const char SymmetricFissionLight;
  static const char SymmetricFissionHeavy;
  static const char SaddleToScission;
  static const char NonStatistical;

  private:
  void tagDaughters(CNucleus *n, std::string const &parentHistory);
  std::string addToHistory(char steptype, std::string const &prevhist);

  static class HistoryStringToDigitsTranslator {
    typedef std::map<char, char> StepMap;
    public:
      HistoryStringToDigitsTranslator() {
        theStepMap[Evaporation] = '1';
        theStepMap[EvaporationResidue] = '2';
        theStepMap[Multifragmentation] = '3';
        theStepMap[AsymmetricFissionLight] = '4';
        theStepMap[AsymmetricFissionHeavy] = '5';
        theStepMap[SymmetricFissionLight] = '6';
        theStepMap[SymmetricFissionHeavy] = '7';
        theStepMap[SaddleToScission] = '8';
        theStepMap[NonStatistical] = '9';
      }
      char operator()(char const c);
    private:
      StepMap theStepMap;
  } theTranslator;

  static const int maxInt32Len;
  HistoryMap theMap;
  int maxEvapZ; //!< maximum Z for evaporation
  static CEvap *evap;
};
#endif // chistory_h
