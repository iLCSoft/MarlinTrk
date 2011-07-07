#ifndef IMarlinTrkFitter_h
#define IMarlinTrkFitter_h

#include <exception>

namespace MarlinTrk{
  class IMarlinTrack ;
}


namespace MarlinTrk{

  
  
  class IMarlinTrkFitter {
    
  public:
    
    virtual ~IMarlinTrkFitter() {};
    
    // initialise track fitter system
    virtual void init() = 0 ;
    
    // instantiate its implementation of the IMarlinTrack 
    virtual MarlinTrk::IMarlinTrack* createTrack() = 0 ;

    // take multiple scattering into account during the fit
    virtual void includeMultipleScattering( bool msOn ) = 0 ;

    // take energy loss into account during the fit
    virtual void includeEnergyLoss( bool energyLossOn ) = 0 ;

    
        
  protected:
    
  private:
    
    IMarlinTrkFitter& operator=( const IMarlinTrkFitter&) ; // disallow assignment operator 
      
  } ;
  
} // end of MarlinTrk namespace 

#endif

