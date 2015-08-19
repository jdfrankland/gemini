/**
 *!\brief weighted Monte Carlo
 *
 * Class CWeight is a base class that deals with a weighted 
 * Monte Carlo scheme. It is used to enhance the probabilty of IMF emission.
 * To compensate for this, each event is given a weight.
 * This weight should be used when histogramming events.
 *
 */

class CWeight
{
 protected:

  float fact;    //!< weighting factor
  int iWeight;  //!< ==0, no weighting
  float runningWeight; //!< running weight of event
  void findFactor(float Glight, float Gimf, float Gfission, float Ggamma);
 public:

  int  chooseChannel(float Glight, float Gimf, float Gfission, float Ggamma, 
                                    float xran);
  void setWeightIMF();
  float getWeightFactor();
};
