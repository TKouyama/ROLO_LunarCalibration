;;
;; ROLOに必要な観測ジオメトリを返すIDL関数
;; 必要なKernelはあらかじめ読み込んであることを仮定
;;
;; lunar_pos.tmの記述を参照
;; CSPICE_FURNSH, 'lunar_pos.tm'
;; で読み込める. 各衛星の軌道SPKは別途個別用に用意
;;

;; Author: Toru Kouyama (t.kouyama@aist.go.jp)

;;
;; 2015.08.12: phase angleを出力に追加
;; 2019.11.15: LEO用に清書.
;;
function get_geometry_with_spice_ROLO_for_LEO $
           ,Obs_date $                  ;;  観測時間
           ,obs_sat_name=obs_sat_name $ ;;  spkに登録した衛星名
           ,target_body_name=target_body_name $ ;; 観測対象 基本的に月, i.e. "Moon"
           ,target_frame_name=target_frame_name ;; 観測対象の持つ座標系, i.e. "Moon_ME"

  ;; 1AU
  AU = 149597870700d0/1000d0 ;; 1 AU (km)

  ;; Date
  if n_elements(Obs_date) eq 0 then begin
    print,"Date is not defined."
    stop
  endif

  utctim = Obs_date; ex '2003-04-14T22:10:00'
  ;; UTC => Et
  CSPICE_STR2ET, utctim, et

  ;;
  ;; Lunar position seen from a satellite in ITRF93 frame (Earth frame)
  ;; with light-time and path corrections at et
  ;;
  CSPICE_SPKPOS, target_body_name , et, 'ITRF93', 'LT+S',obs_sat_name, lunar_pos, ltime_LT
  LO_distance = sqrt(total(lunar_pos^2.))

  ;;
  ;; Solar position seen from Moon_ME frame with light-time and path corrections
  ;; at "et-ltime_LE"
  CSPICE_SPKPOS, 'Sun' , (et-ltime_LT), target_frame_name, 'LT+S',target_body_name, Sun_pos_lunar, ltime_SL
  SL_distance = sqrt(total(sun_pos_lunar^2.))
  SL_distance_au = SL_distance/AU

  ;;
  ;; Sub solar point at the obervation
  ;;
  CSPICE_SUBSLR, "Intercept: ellipsoid",target_body_name,et,target_frame_name, "LT+S", $
                 obs_sat_name, subpoint_solar, subepc, subsrfvec_solar
  ;; Convert to longitude, latitude
  CSPICE_RECLAT, subpoint_solar, radii, sub_solar_lon, sub_solar_lat

  ;;
  ;; Sub observer point on the Moon
  ;;
  CSPICE_SUBPNT, "Intercept: ellipsoid",target_body_name,et,target_frame_name, "LT+S",obs_sat_name,$
                  subpoint_terra, subepc, subsrfvec_terra
  ;; Convert to longitude, latitude
  CSPICE_RECLAT, subpoint_terra, t_radii, sub_terra_lon, sub_terra_lat

  ;;
  ;; Location of observer in Earth frame
  ;;
  ;; 地球フレーム上でのObserver位置 ;; 近いのと光時補正の影響が良くわからないので "NONE"
  CSPICE_SPKPOS, obs_sat_name , et, 'ITRF93', 'NONE','EARTH', Terra_pos, ltime
  ;; Sub Terra lon lat on Earth
  CSPICE_RECLAT, Terra_pos, ET_ditance, sub_terra_lon_e, sub_terra_lat_e

  ;;
  ;; Phase angle
  ;;
  u_subpoint_solar = subpoint_solar/sqrt(total(subpoint_solar^2.))
  u_subpoint_terra = subpoint_terra/sqrt(total(subpoint_terra^2.))
  phase_angle_deg = acos(total(u_subpoint_terra*u_subpoint_solar))*180./!dpi


  ;; 結果をまとめた構造体を返す. 神山IDL用
  tags = ['SL_distance','LO_distance','ssl_longitude','ssl_latitude' $
          ,'ssc_longitude', 'ssc_latitude','phase']
  obs_geo = create_struct(tags $
    ,SL_distance $
    ,LO_distance $
    ,sub_solar_lon*180./!dpi $
    ,sub_solar_lat*180./!dpi $
    ,sub_terra_lon*180./!dpi $
    ,sub_terra_lat*180./!dpi $
    ,phase_angle_deg $
    )

  return, obs_geo
end
