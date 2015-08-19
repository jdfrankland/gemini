#include "CHistory.h"

CEvap *CHistory::evap;

CHistory::HistoryStringToDigitsTranslator CHistory::theTranslator;

const int CHistory::maxInt32Len = numeric_limits<int32_t>::digits10;

const char CHistory::Evaporation = 'e';
const char CHistory::EvaporationResidue = 'E';
const char CHistory::Multifragmentation = 'm';
const char CHistory::AsymmetricFissionLight = 'a';
const char CHistory::AsymmetricFissionHeavy = 'A';
const char CHistory::SymmetricFissionLight = 'f';
const char CHistory::SymmetricFissionHeavy = 'F';
const char CHistory::SaddleToScission = 's';
const char CHistory::NonStatistical = 'n';

void CHistory::tagDaughters(CNucleus *n, std::string const &parentHistory) {
  CNucleus *daughterLight = n->getLightDaughter();
  CNucleus *daughterHeavy = n->getHeavyDaughter();

  if(!daughterLight && !daughterHeavy) {
    return;
  }

  // Saddle-to-scission transition or gamma decay
  if(!daughterLight) {
    // Gamma decay
    theMap[daughterHeavy] = parentHistory;
    tagDaughters(daughterHeavy, parentHistory);
    return;
  }

  std::string lightHistory, heavyHistory;

  if(n->isSaddleToScission()) {

    if(daughterLight->iZ <= maxEvapZ) {
      lightHistory = addToHistory(Evaporation, addToHistory(SaddleToScission, parentHistory));
      heavyHistory = addToHistory(SaddleToScission, parentHistory);
    } else {
      lightHistory = addToHistory(SymmetricFissionLight, parentHistory);
      heavyHistory = addToHistory(SymmetricFissionHeavy, parentHistory);
    }
    theMap[daughterLight] = lightHistory;
    theMap[daughterHeavy] = heavyHistory;
  } else {
    if(n->isNotStatistical()) {
      lightHistory = addToHistory(NonStatistical, parentHistory);
      heavyHistory = addToHistory(NonStatistical, parentHistory);
    } else if(daughterLight->iZ <= maxEvapZ) {
      lightHistory = addToHistory(Evaporation, parentHistory);
      heavyHistory = addToHistory(EvaporationResidue, parentHistory);
    } else {
      lightHistory = addToHistory(AsymmetricFissionLight, parentHistory);
      heavyHistory = addToHistory(AsymmetricFissionHeavy, parentHistory);
    }
    theMap[daughterLight] = lightHistory;
    theMap[daughterHeavy] = heavyHistory;
  }
  tagDaughters(daughterLight, lightHistory);
  tagDaughters(daughterHeavy, heavyHistory);
}

//************************************************************
/**
 * add a step to the history variable
 */
std::string CHistory::addToHistory(char steptype, std::string const &prevhist)
{
  std::string history;

  size_t histLenMinus1 = prevhist.length();
  char last;
  if(histLenMinus1>0) {
    histLenMinus1--;
    last = prevhist.at(histLenMinus1);
  } else {
    last = '\0';
  }

  if( last == EvaporationResidue && steptype == Evaporation )
    history = prevhist.substr(0,histLenMinus1) + Evaporation;
  else if( last == EvaporationResidue && steptype == EvaporationResidue )
    history = prevhist;
  else if( last == SaddleToScission && steptype == SaddleToScission )
    history = prevhist;
  else
    history = prevhist + steptype;

  return history;
}

char CHistory::HistoryStringToDigitsTranslator::operator()(char const c) {
  std::map<char,char>::const_iterator iter = theStepMap.find(c);
  if(iter!=theStepMap.end())
    return iter->second;
  else
    return '\0';
}

