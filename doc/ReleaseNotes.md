# v02-09

* 2021-09-29 Bohdan Dudar ([PR#22](https://github.com/iLCSoft/MarlinTrk/pull/22))
  - Now checking for the intersection with the endcaps even if intersection with the barrel has been found. The closest track state to the last tracker hit is saved

* 2021-09-27 Thomas Madlener ([PR#23](https://github.com/iLCSoft/MarlinTrk/pull/23))
  - Migrate CI from travis to github actions.

# v02-08

* 2019-08-23 Andre Sailer ([PR#16](https://github.com/iLCSoft/MarlinTrk/pull/16))
  - CMake: drop export of Lib dependencies, which doesn't work in cmake 3.12
  - CMake: fix CLHEP discovery, first find CLHEP then DD4hep (and implicitly geant4)

# v02-07

* 2018-03-13 Marko Petric ([PR#11](https://github.com/iLCSoft/MarlinTrk/pull/11))
  - Fix for iLCSoft/LCIO#35

* 2018-03-23 Frank Gaede ([PR#12](https://github.com/iLCSoft/MarlinTrk/pull/12))
  - use the origin as reference point in MarlinTrk::createPrefit() for initial helix
      - this improves the fit probability for Si-tracks (in ILD)
      - similar issues reported by CLICdp

* 2018-03-28 Marko Petric ([PR#13](https://github.com/iLCSoft/MarlinTrk/pull/13))
  - Fix for the removal of DDSurfaces which have been merged into DDRec 
    -  includes from `DDSurfaces` -> `DDRec`
    - namespace `DDSurfaces` -> `dd4hep::rec`

* 2017-11-30 Andre Sailer ([PR#10](https://github.com/iLCSoft/MarlinTrk/pull/10))
  - Performance optimisation, avoiding dynamic_cast, see also  AIDAsoft/aidaTT#19
  - Fix compiler warnings for uninitialised members

# v02-06

* 2017-10-12 Frank Gaede ([PR#9](https://github.com/ilcsoft/MarlinTrk/pull/9))
  - fix the re-setting of eLoss and QMS configuration parameters in TrkSysConfig

# v02-05

* 2017-10-12 Shaojun Lu ([PR#8](https://github.com/iLCSoft/MarlinTrk/pull/8))
  - Fix the segmentation fault which is caused by assigning the wrong number to _indexMap.
      - Only when trajectory added measurement successfully, then to increase pointLabel by one, and assign it to _indexMap.
      - Both DDKalTest and AidaTT(GBL) may be used by  "RefitProcessor" and "FullLDCTracking_MarlinTrk" for track fitting.

* 2017-10-12 Frank Gaede ([PR#7](https://github.com/iLCSoft/MarlinTrk/pull/7))
  - add new helper class TrkSysConfig for setting the correct system configuration per event (scope)
  - reduce initial error on tanLambda in `MarlinDDKalTestTrack::initialise()`

* 2017-09-11 Frank Gaede ([PR#6](https://github.com/iLCSoft/MarlinTrk/pull/6))
  - reduce initial uncertainty for kappa/omega for track fits in
     `MarlinDDKalTestTracks::initialise(bool fitDirection)` 
        - fixes a problem in TPC track fits, when first three points define a very different curvature
        - should have no effect on Si-Tracking

# v02-04

* 2017-06-20 Andre Sailer ([PR#4](https://github.com/iLCSoft/MarlinTrk/pull/4))
  - Adapt to namespace changes in DD4hep

* 2017-06-27 Andre Sailer ([PR#5](https://github.com/iLCSoft/MarlinTrk/pull/5))
  - Clean up of tkrSystem and DDKalDetectors at the end of lifetime

* 2017-05-17 Frank Gaede ([PR#2](https://github.com/iLCSoft/MarlinTrk/pull/2))
  - replace gear::Vector3D with DDSurfaces::Vector3D in IMarlinTrack interface
        - provide wrapper functions with old signature for backwards compatibility
        - remove all other usages of Gear
  - remove obsolete MarlinKalTest and MarlinKalTestTrack
         - from now on use only MarlinDDKalTest
  - fix all warnings seen with llvm (on mac)

# v02-03

# v02-02
F. Gaede
* added setMass()/getMass to IMarlinTrack interface
* implemented for KalTest, DDKalTest and aidaTT tracks
* made compatible with c++11
* removed -ansi -pedantic -Wno-long-long
* fixed narrowing in initializer lists

v02-01
F.Gaede
* MarlinTrkUtils.cc
* fix condition for setting the z0 to 0 for calo TS for curling tracks 
* generalized the subdetectorHitNumbers to work for any detector (assuming that the tracker system ids are smaller than the Ecal system id !)
* fix bug with track state at cale endcap for old tracking ( use lcio::ILDDetID::ECAL also for endcap ) 
* MarlinAidaTTTrack: added a prefit to the initialize: use up to ~25 hits for an initial fit (w/o material) before doing the actual fit with all hits 
* added IMarlinTrkSystem::name() in order to allow tracking utilities to adapt to special features of fitters 
* allow to fall back to IMarlinTrkTrack implementations initialise() in createFinalisedLCIOTrack() -> for aidaTT see procedure above
* set precision of 2nd direction v to 0. for 1D hits
* ignore virtual surface with no material (e.g. inside the beam pipe)
* fixed issue w/ result labels 
* fixed point labels if QMS is turned off
* treat 1D measurements correctly
* finalized implementation of IMarlinTrkTrack interface
* can be used in RefittProcessor 
* added first implementation if IMarlinTrk for aidaTT
* needs debugging ...
* MarlinDDKalTest.cc
* fixed track state at calo creation 
* consider only DDCylinderMeasLayers as IP layers 

# v02-00-01
* made createTrackStateAtCaloFace() work for DD4hep (and old Mokka/Gear ) based tracking

# v02-00
* added implementation of MarlinDDKalTest: KalTest implementation of MarlinTrk working with DD4hep surfaces and DDKalTest
* changed Marlin(DD)KalTestTrack::initialise( const EVENT::TrackState& ts, double bfield_z, bool fitDirection ) to ignore the bfield_z parameter and get it from the first hit
* fixes issue in MarlinTrkUtils/finaliseLCIOTrack() where the method above is used  but the bfield is not available
* should consider changing IMarlinTrkTrack interface ...
* use new LCDD detector lists to create DDKalDetectors only for "tracker" and "passive" sub detectors, i.e. calorimeters are ignored
* made MarlinTrk::Factory a singleton
* sub sequent calls to init() can only change the options (MSQ,Eloss)
* Factory caches one IMarlinTrkSystem of every type 
* USERS SHOULD NO LONGER DELETE THE IMarlinTrkSystem POINTER IN THEIR CODE (Marlin processor)
 
* fixed MarlinTrkUtils::createFinalisedLCIOTrack for fitting insode out (forward ):
* smooth back to last constrained hit 
* add innermost hits with a proper Kalman filter step to get track state at the IP
* added fit_direction to finaliseLCIOTrack() et al
* made treatment of fit_direction more consistent

* added an IP layer to MarlinDDKalTest ( DDCylinderMeasLayer w/ smallest radius )
* debug printout in MarlinTrkUtils

* ignore missing track state at calo face (for now) -> just print warning
* added debug printout

# v01-11
 * Improved readability of the MarlinTrkUtils functions
* Added LCTPC-specific modifications:
* geometry building is now based on the detector name from GEAR
* Magnetic field may be zero -> avoid the value of infinity for kappa
* Interaction Point (IP) meas. layer may not exist for LCTPC
* Do not throw exception if the layer is not found, return "no_intersection" code instead. This avoids the problem if some of the state is not or cannot be available, e.g. the IP and CaloFace states are not (yet?) defined for the LCTPC geometry

# v01-10-01
* Added some more diagnostics, no change to algorithms.

# v01-10
* General 
* Made Debug output more consistent

* MarlinKalTestTrack
* Changed smooth() to do _kaltrack->SmoothAll(), previously _hitIndexAtPositiveNDF + 1 
* Corrected orientation regarding transporting inwards or outwards in propagate.
* Fixed problem where initial covariance term kappa,tanL wrongly set as kappa,z0. Minimal impact excepted. (Tino Calancha) 
* Ensure that hits which are rejected for reasons other that Chi2 cut are added to the list of outliers. 

# v01-09
* General 
* Made Debug output more consistent

* IMarlinTrack
* Moved constant definitions of return codes outside of the header file. 
* Added additional error return, for the case where no site are filtered when calling fit.

* MarlinKalTestTrack
* Prefer translation over rotation of the trackstate early in the fit, when using simple helix initialisation.

* MarlinTrkUtils
* Use EVENT::TrackState in place of IMPL::TrackStateImpl where appropriate. 
* The TrackState for the initialisation can be at any reference point.

# v01-08
Improvements in Diagnostics.

* MarlinKalTestTrack
* Print trkhit id in add hit

* TrkAnalysisTree
* Initialise the index values with -1
* Added run and event number, added reciprocal weights for relations, corrected array sizes in places.

# v01-07
* IMarlinTrack
* Added two new methods to IMarlinTrack interface: int getNDF( int& ndf ) and int getTrackerHitAtPositiveNDF( EVENT::TrackerHit*& trkhit ). 
* MarlinKalTestTrack
* smooth() no longer smooths back to 4th hit, but goes back the the hit at which the fit was constrained, meaning that it can adapt to 1D hits properly.
* Fixed the issue where the upper part of the triangular symmetric cov matrix was not copied into the square TMatrixD cov(5,5).

* HelixFit
* Correct pass by reference in argument list

* MarlinTrkUtils
* placed definition of finaliseLCIOTrack in header making it publicly available. 
* Corrected track state at first and last hit.
* finaliseLCIOTrack now has the option to take the state at the last hit and calo face as input, for the case where these have been calculated separately, e.g. in an in-out out-in strategy.
* Due to occasional problems in the covariance matrix when smoothing back to sites added before the fit was fully constrained, the last hit and the calo face use the smoothed values from the site at which the fit was fully constrained. The fit is propagated to the point of closest approach to the last hit, and to the calo face. 	 
* Improved Diagnostics.

# v01-06 
* Added convenient utility methods for fitting and producing LCIO tracks.
* Added SET.
* Moved initial pivot from the origin to the point where the initial helix crossed the layer containing the first hit to be filtered. 
* Replaced NULL C macro with 0 for pointers value.
* Corrected phi range for moveRefPoint in HelixTrack.
* Corrected the direction of the helix fit in HelixTrack
* Added diagnostic functionality. This needs to be enabled via the #define statement MARLINTRK_DIAGNOSTICS_ON in MarlinTrkDiagnostics.h. The recording of diagnostics is controlled via DiagnosticsController, with recording off by default. To enable recording use getDiagnositicsPointer on the IMarlinTrkSystem and cast to DiagnosticsController, then call init(std::string root_file_name, std::string root_tree_name, bool recording_on) with recording_on set to true.

# v01-05
* MarlinKalTestTrack:
* Protect against initializing before any hits have been added.
* Allow starting from strip hit in debug mode.

* HelixTrack:
* Take account of direction of helix for fit.

* LCIOTrackPropagators:
* Corrected phi range for propagation.

# v01-04 
* Changed fit method to optionally respect a maxDeltaChi2 in IMarlinTrack interface: fit( double maxChi2Increment=DBL_MAX ) ;
* Added new constructor to HelixTrack that generates helix from three point.
* Corrected sign of d0 when initialising from an LCIO track state in MarlinKalTestTrack.
* Fill GEAR MeasurementSurfaceStore using ILDMeasurementSurfaceStoreFiller from KalDet.
* Cleaned up depedencies

# v01-03
* Added cylindrical SIT for LOI data.
* improve use of MarlinTrk namespace.
* Fixed use of covariance matrix in initialise. 
* Changed getSiteFromLCIOHit from taking an iterator to taking a pointer to a site. 
* Corrected smoothed state when used.

# v01-02
* To correctly determine the fits at the first hit and the last hit methods have be added to IMarlinTrk to get the list of pointers to EVENT::TrackerHit and the chi2 increment for both hits included in the fit as well as those rejected as outliers. 
* To test the chi2 increment which would result from adding a hit to the fit the method testChi2Increment has been added to IMarlinTrk. This method will not alter the fit, it simply provides a method to test the inclusion of a hit. 
* MarlinKalTestTrack has been adapted to use the pointer to LCIO TrackerHits in ILDVTrackHit.


# v01-01
* Support for KalDet Multilayer added
* Added calo face to support detector
* Made use of new method in ILDVMeasLayer "getIntersectionAndCellID"
* Change GetMLName to GetName in line with KalTest
* Fixed phi range of helixtrack
* MarlinKalTestTrack fixed definition of phi0 in initialise(TS& ts) corrected conversion from LCIO to KalTest track parameters in initialise

# v01-00
* IMarlinTrack: Interface for generic tracks in MarlinTrk. The interface should provide the functionality to perform track finding and fitting. 
* IMarlinTrkSystem: Base class for tracking system implementations in MarlinTrk.
* MarlinKalTestTrack: Tracking System used to create a KaltTest Kalman fitter* instantiates and holds the detector geometry. 
* LCIOTrackPropagators: Functions to perform the geometric propagation of LCIO tracks, for closed solutions.
* HelixFit: C++ rewrite of the aleph Fortran routine TFITHL, which provides a fast helix fit

	
