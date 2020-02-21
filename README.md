# ROLO_LunarCalibration
IDL programs for estimating lunar irradiance with the ROLO model at a given date. These programs were used in:
- Kouyama et al., Development of an application scheme for the SELENE/SP lunar reflectance model for radiometric calibration of hyperspectral and multispectral sensors, Planetary and Space Science, 2016.
- Kouyama et al., Lunar calibration for ASTER VNIR and TIR with observations of the Moon in 2003 and 2017, Remote Sensing, accepted.

This repository will be composed of:
- Program for converting TLE to SPK kernel
- Program for calculating obsrevation geometry at a given date using SPICE toolkit
- Program for estimating lunar irradiance at the given geometry with the ROLO model

(Not completed yet.)

Usage example:
```
  obs_date = '2019-12-12T00:00:00'
  obs_sat_name = 'Terra'
  target_body_name = 'Moon'
  target_frame_name = 'Moon_ME'

  obs_geo = get_geometry_with_spice_rolo_for_LEO(obs_date $
               ,obs_sat_name=obs_sat_name $
               ,target_body_name=target_body_name $
               ,target_frame_name=target_frame_name)

  sub_obs_lat_deg = obs_geo.ssc_latitude
  sub_obs_lon_deg = obs_geo.ssc_longitude
  sub_solar_lon = obs_geo.ssl_longitude * !dpi/180.
  sub_solar_lat = obs_geo.ssl_latitude * !dpi/180.

  ;; Distance parameters
  AU = 149597870700d0/1000d0 ;; 1 AU (km)
  Moon_D = 384400d ;; km

  ;; Normalizing distance
  sl_distance = obs_geo.sl_distance / AU
  lo_distance = obs_geo.lo_distance / Moon_D
```
