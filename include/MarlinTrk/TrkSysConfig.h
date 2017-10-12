#ifndef TrkSysConfig_h
#define TrkSysConfig_h


namespace MarlinTrk{

  class IMarlinTrkSystem ;

  /** Helper class for temporarilly setting a configuration option of the tracking system for
   *  the current scope. The old seeting will be restored, when the class
   *  gets out of scope.  
   *
   * @author F. Gaede DESY
   * @date 10/2017
   */


  template <unsigned CFG>
  class TrkSysConfig{
    
  public:
    
    TrkSysConfig() = delete ;
    
    TrkSysConfig(IMarlinTrkSystem* trksys, bool value) : _trksys( trksys ){

      _value = _trksys->getOption(CFG) ;
      _trksys->setOption(CFG,value) ;
    }

    ~TrkSysConfig(){
      _trksys->setOption(CFG,_value) ;
    }

  protected:
    IMarlinTrkSystem* _trksys = 0;    
    bool _value = false ;
  } ;
 

}
#endif
