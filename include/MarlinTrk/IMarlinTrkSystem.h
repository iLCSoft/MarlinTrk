#ifndef IMarlinTrkSystem_h
#define IMarlinTrkSystem_h

#include "MarlinTrkDiagnostics.h"

#include <exception>
#include "ConfigFlags.h"

namespace MarlinTrk{
  class IMarlinTrack ;
}


namespace MarlinTrk {
  
  /** Exception thrown in IMarlinTrk namespace (implemetations of IMarlinTrkSystem 
   *  and IMarlinTrack).
   */
  class Exception : public std::exception {
    
  protected:
    std::string message ;
    
    Exception(){  /*no_op*/ ; }
    
  public:
    virtual ~Exception() throw() { /*no_op*/; }
    
    Exception( const std::string& text ){
      message = "MarlinTrk::Exception: " + text ;
    }
    
    virtual const char* what() const  throw() { return  message.c_str() ; }
    
  };
  
  //----------------------------------------------------------------------------------------------------
  
  /** Base class for tracking system implementations in MarlinTrk.
   *
   * @version $Id$
   * @author S.Aplin, F. Gaede DESY
   */
  
  class IMarlinTrkSystem {
    
  public:
    
    /** 'Enums' for configuration options to be used with setOption().
     */
    struct CFG {
      /** Use multiple scattering in the track fits. */
      static const unsigned  useQMS   = 1 ;
      /** Use multiple scattering in the track fits. */
      static const unsigned  usedEdx  = 2 ;
      /** Use smoothing when calling fit( bool fitDirection ) */
      static const unsigned  useSmoothing = 3 ;
      //---
      static const unsigned  size     = 4 ;
      
    } ;
    
    
    /** D'tor - cleans up any allocated resources.*/
    virtual ~IMarlinTrkSystem() {};
    
    
    /** Sets the specified option ( one of the constants defined in IMarlinTrkSystem::CFG ) 
     *  to the given value.
     */
    void setOption(unsigned CFGOption, bool val) ;  
    
    /** Return the option's current value - false if option not defined.
     */
    bool getOption( unsigned CFGOption) ;
    
    /** String with all configuration options and their current values. 
     */ 
    std::string getOptions() ;
    
    
    /** Initialise tracking system - to be called after configuration with setOption() -
     *  IMarlinTrkSystem cannot be used before a call to init().
     */
    virtual void init() = 0 ;
    
    
    /** Return an instance of IMarlinTrack corresponding to the current implementation.
     */
    virtual MarlinTrk::IMarlinTrack* createTrack() = 0 ;
    
    
#ifdef MARLINTRK_DIAGNOSTICS_ON
    
    /** Return the pointer to the Diagnositics Object. Forseen for internal diagnostics, only available when complied with MARLINTRK_DIAGNOSTICS_ON defined. 
     */
    virtual void * getDiagnositicsPointer() = 0 ;

#endif
    
  protected:
    
    ConfigFlags _cfg ;
    
    /** Register the possible configuration options
     */ 
    void registerOptions() ;
    
    
    
  private:
    
    IMarlinTrkSystem& operator=( const IMarlinTrkSystem&) ; // disallow assignment operator 
    
    
  } ;
  
} // end of MarlinTrk namespace 

#endif

