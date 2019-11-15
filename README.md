# ROLO_LunarCalibration
IDL programs for estimating lunar irradiance with ROLO model at given date.

Usage example
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
