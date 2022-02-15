;;
;; ROLOの月照度を求める
;; SP_vs_ROLO_phase_angle
;; 観測者座標を自分で入れることが可能
;; 2021.07.16
;;

forward_function get_geometry_with_spice_Rolo_location
forward_function get_params_rolo
forward_function calculate_rolo_reflectance
forward_function generate_rolo_spec


;;
;; 設定ファイルを読み込む
;;
function read_cl_file_for_AIST_ROLO,ifname_cl_file $
                                    ,ifldname_kernel $
                                    ,metakr $
                                    ,ifname_ROLO_csv $
                                    ,ifname_soil $
                                    ,ifname_rock $
                                    ,solfname_rolo $
                                    ,obs_sat_spk $
                                    ,ifname_tle $
                                    ,obs_date $
                                    ,ofldname $
                                    ,ofname_SPK $
                                    ,log_output=log_output
  err = 0

  openr,1,ifname_cl_file, ERROR = err
  if err ne 0 then begin
    close,1
    return,err
  endif

  tmp_line=''
  while(not eof(1)) do begin
    readf,1,tmp_line
    ;print,tmp_line
    if strmid(tmp_line,0,1) eq '#' then begin
      ;; #はコメントアウト
    endif else begin
      tmp_line = strcompress(tmp_line,/remove)
      result = strpos(tmp_line,'=')
      if result ne -1 then begin
        tmp_str = strmid(tmp_line,0,result)
        tmp_len = strlen(tmp_line)
        tmp_value = strmid(tmp_line,result+1,tmp_len-result)

        case tmp_str of
          'KERNEL_file_directory': ifldname_kernel = tmp_value
          'METAKR_file_name': METAKR = tmp_value
          'ROLO_param_path': ifname_ROLO_csv = tmp_value
          'SOIL_path': ifname_soil = tmp_value
          'ROCK_path': ifname_rock = tmp_value
          'SOL_model_path': solfname_rolo = tmp_value
          'SPK_path': obs_sat_spk = tmp_value
          'TLE_file_name': ifname_tle = tmp_value
          'OBS_date': obs_date = tmp_value
          'OUTPUT_directory': ofldname = tmp_value
          else: print,"!!",tmp_value
        endcase

        if n_elements(log_output) then AIST_ROLO_log_print,tmp_line

      endif
    endelse

  endwhile
  close,1

  if n_elements(log_output) then AIST_ROLO_log_print,""

  ;; ファイル名からオブジェクト名作成 ;;
  tle_object = file_basename(ifname_tle, '.tle')

  ;; Outputするファイル名 ;;
  ofldname_spk =file_dirname(ifname_tle)+'\'
  ofname_SPK = ofldname_spk + tle_object + '.bsp'

  if n_elements(fldname_kernel) eq 0 or n_elements(METAKR) eq 0 then err = -1
 
  
  return,err
end

;;
;; main
;;
pro estimate_ROLO_brignthess_band,ifname_cl_file, log_output=log_output

  if n_elements(ifname_cl_file) eq 0 then begin
    ifname_cl_file = 'C:\work\ENVI_IDL\LunaCal\AIST_ROLO\AIST_ROLO_cal_file.txt'
  endif
  
  wrkdname = file_dirname(ifname_cl_file,/MARK_DIRECTORY)
  cd, wrkdname,current=current_fld
  
  result= read_cl_file_for_AIST_ROLO(ifname_cl_file $
                                    ,ifldname_kernel $
                                    ,metakr $
                                    ,ifname_ROLO_csv $
                                    ,ifname_soil $
                                    ,ifname_rock $
                                    ,solfname_rolo $
                                    ,obs_sat_spk $
                                    ,ifname_tle $
                                    ,obs_date $
                                    ,ofldname $
                                    ,ofname_SPK $
                                    ,log_output=log_output)


  ;; 日付の整形, "-"などを取り除く
  str_obs_date = strmid(obs_date,0,4)+strmid(obs_date,5,2)+strmid(obs_date,8,2) $
                 +'T'+strmid(obs_date,11,2)+strmid(obs_date,14,2)+strmid(obs_date,17,2)

  ;;
  ;; Output file name ;;
  ;;

  ofname = ofldname + path_sep() + 'ROLO_output_reflectance_band.csv'
  ofname_irradiance = ofldname + path_sep() + 'ROLO_output_irradiance_band.csv'
  ofname_irradiance_spec = ofldname + path_sep() + 'ROLO_output_irradiance_spec.csv'

  ;; 日付を含めたいとき ;;
  ;ofname = ofldname + path_sep() + 'ROLO_output_reflectance_'+str_obs_date+'.csv'
  ;ofname_irradiance = ofldname + path_sep() + 'ROLO_output_irradiance_'+str_obs_date+'.csv'

  ;;
  ;; Read ROLO parameters ;;
  ;;
  result = read_csv(ifname_ROLO_csv,count=csv_count)
  wav = result.field01
  wav_pos = where(wav gt 300 and wav lt 2500)
  wav = wav[wav_pos]
  ;print,long(wav)

  ;; Solar spectrum for ROLO center wav ;;
  get_radiance_rolo,sol_rad_rolo,solfname=solfname_rolo,wav=wav_rolo
  sol_rad_rolo /= 1000.

  ;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Spice Kernel読み込み  ;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; Genelar kernel
  cd, ifldname_kernel, current = old_dir
  CSPICE_FURNSH, METAKR
  cd, old_dir

  ;; Terra spk (position kernel)
  ;cd, obsvrkr_dir
  CSPICE_FURNSH, obs_sat_spk
  ;cd, old_dir

  ;; ID をspk fileから取り出す ;;
  winsiz = 2000l
  maxobj = 1000l
  cover = cspice_celld( WINSIZ )
  ids   = cspice_celli( MAXOBJ )
  cspice_spkobj, obs_sat_spk, ids
  obj = ids.base[ ids.data + 0 ]
  
  ;obs_sat_name = "EARTH" ;; 必要なカーネルがあればTerra以外でも計算可能
  ;obs_sat_name = '-999' ;; 必要なカーネルがあればTerra以外でも計算可能
  obs_sat_name = strtrim(string(obj),2)
  target_body_name = "MOON"
  target_frame_name = "MOON_ME"

  ;; 1 AU
  AU = 149597870700d0/1000d0 ;; 1 AU (km)

  ;; 月視直径 ;;
  Omega_M = 6.4177e-5

  ;;
  ;; Write header strings ;;
  ;;
  openw,1,ofname
  printf,1,'Date,Sun-Moon_distance,Moon-Obs_distance,Phase_angle,Sub solar longitude, Sub observer longitude, Sub observer latitude, 350,355,405,412,414,441,465,475,486,544,549,553,665, 693,703,745,763,774,865,872,882,928,939,942,1059,1243,1538,1633,1981,2126,2250,2383'
  close,1

  ;; Write header strings ;;
  openw,2,ofname_irradiance
  printf,2,'Date,Sun-Moon_distance,Moon-Obs_distance,Phase_angle,Sub solar longitude, Sub observer longitude, Sub observer latitude, 350,355,405,412,414,441,465,475,486,544,549,553,665, 693,703,745,763,774,865,872,882,928,939,942,1059,1243,1538,1633,1981,2126,2250,2383'
  close,2

  ;; Write header strings ;;
  openw,3,ofname_irradiance_spec
  printf,3,'Date,Sun-Moon_distance,Moon-Obs_distance,Phase_angle,Sub solar longitude, Sub observer longitude, Sub observer latitude'
  close,3


  for i=0, n_elements(obs_date)-1,1 do begin
    tmp_obs_date = obs_date[i]

    obs_geo = get_geometry_with_spice_rolo_location(tmp_obs_date $
      ,obs_sat_name=obs_sat_name $
      ,target_body_name=target_body_name $
      ;; 観測者位置を入れる ;;
      ,obs_lat_earth=obs_lat_earth $
      ,obs_lon_earth= obs_lon_earth $
      ,obs_alt_earth=obs_alt_earth $
      ,target_frame_name=target_frame_name)

    sub_obs_lat_deg = obs_geo.ssc_latitude
    sub_obs_lon_deg = obs_geo.ssc_longitude
    sub_solar_lon = obs_geo.ssl_longitude * !dpi/180.
    sub_solar_lat = obs_geo.ssl_latitude * !dpi/180.
    SL_distance = obs_geo.ST_distance /AU
    LO_distance = obs_geo.TO_distance
    ;stop


    ;; 距離補正項 ;;
    f_d = (SL_distance/1.)^2. * (LO_distance/384400.)^2.


    ;;
    ;; 位相角の計算
    ;; 符号は満ちの時負、欠けの時正
    ;;
    obs_direction = [cos(sub_obs_lat_deg*!dpi/180.) * cos(sub_obs_lon_deg*!dpi/180.) $
                    ,cos(sub_obs_lat_deg*!dpi/180.) * sin(sub_obs_lon_deg*!dpi/180.) $
                    ,sin(sub_obs_lat_deg*!dpi/180.)]

    sun_direction = [cos(sub_solar_lat) * cos(sub_solar_lon) $
                    ,cos(sub_solar_lat) * sin(sub_solar_lon) $
                    ,sin(sub_solar_lat)]

    phase_angle_deg = abs(acos(total(obs_direction*sun_direction))*180./!dpi)


    print,tmp_obs_date
    print,"Phase angle: ",phase_angle_deg, sub_solar_lon*180d/!dpi


    ;; ROLO 係数の読み取り ;;
    disk_ref = dblarr(n_elements(wav))
    disk_irradiance = dblarr(n_elements(wav))

    for wav_i=0, n_elements(wav)-1, 1 do begin
      sel_wav = wav[wav_i]
      error = get_params_rolo(sel_wav,ifname_ROLO_csv,param_a,param_b,param_d,coeff_c,coeff_p)


      ;; Rolo brightness ;;
      Disk_reflectance_rolo = calculate_rolo_reflectance(param_a,param_b,param_d,coeff_c,coeff_p $
        ,phase_angle_deg $
        ,sub_solar_lon $
        ,sub_obs_lat_deg $
        ,sub_obs_lon_deg)

      disk_ref[wav_i] = disk_reflectance_rolo
      disk_irradiance[wav_i] = (disk_ref[wav_i] * Omega_M * sol_rad_rolo[wav_i]/!dpi) / f_d * 1e6 ;; uW/m2/um

    endfor

    ;;
    ;; ROLOの連続スペクトルを求めるTest 2021.06.09 ;;
    ;;
    wav_n = 2300l
    srf_wav = dindgen(wav_n) + 300.
    ;ROLO_ref = generate_rolo_spec(phase_angle_deg, sub_solar_lon, sub_obs_lat_deg, sub_obs_lon_deg,srf_wav)
    ROLO_ref = generate_rolo_spec(ifname_ROLO_csv, ifname_soil, ifname_rock $
                                  ,phase_angle_deg, sub_solar_lon, sub_obs_lat_deg, sub_obs_lon_deg,srf_wav)

    ;; Solar brightness ;;
    ;; sol_rad_rolo = W/m2/sr/um;;
    get_radiance_rolo,sol_rad_rolo_spec,solfname=solfname_rolo, wav=srf_wav

    ;; ROLO irradiance with distance correction
    ROLO_irradiance_spec = ROLO_ref * Omega_M * sol_rad_rolo_spec / !dpi / f_d

    wset,0
    circlesym,thick=2
    plot,wav,disk_ref,psym=8,thick=1,xr=[300,2500],yr=[0.05,0.6],xs=1,ys=1,/ylog
    oplot,srf_wav,ROLO_ref
    wset,0

    wset,1
    if i eq 0 then begin
      plot,srf_wav,ROLO_irradiance_spec,yr=[0,0.01],ys=1
    endif else begin
      plot,srf_wav,ROLO_irradiance_spec,yr=[0,0.01],ys=1,/noerase
    endelse
    wset,0

    ;; 出力 ;;
    openw,1,ofname,/append
    printf,1 $
      ,tmp_obs_date $
      ,SL_distance $
      ,LO_distance $
      ,phase_angle_deg $
      ,sub_solar_lon*180./!dpi $
      ,sub_obs_lon_deg $
      ,sub_obs_lat_deg $
      ,disk_ref $
      ,format = '(1(A,","), 6(f,","), 32(f,","))'
    close,1

    ;; 出力 ;;
    openw,2,ofname_irradiance,/append
    printf,2 $
      ,tmp_obs_date $
      ,SL_distance $
      ,LO_distance $
      ,phase_angle_deg $
      ,sub_solar_lon*180./!dpi $
      ,sub_obs_lon_deg $
      ,sub_obs_lat_deg $
      ,disk_irradiance $
      ,format = '(1(A,","), 6(f,","), 32(f,","))'
    close,2

    ;;
    openw,3,ofname_irradiance_spec,/append

    printf,3 $
      ,tmp_obs_date $
      ,SL_distance $
      ,LO_distance $
      ,phase_angle_deg $
      ,sub_solar_lon*180./!dpi $
      ,sub_obs_lon_deg $
      ,sub_obs_lat_deg $
      ,format = '(1(A,","), 6(f,","), 32(f,","))'

    printf,3,'Wavelength, irradiance'
    for wav_i=0, n_elements(ROLO_irradiance_spec)-1, 1 do begin
      printf,3, srf_wav[wav_i], ',', ROLO_irradiance_spec[wav_i]
    endfor
    close,3

    wait,0.01
    ;stop


  endfor


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; SPICE kernel のクリア ;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;
  CSPICE_KCLEAR

  return
end

;;
;; 指定観測時の観測者ジオメトリを得る。
;; 地球上任意の地点を指定(緯度経度高度)できるよう拡張版(2017.11.24)
;;
function get_geometry_with_spice_Rolo_location $
  ,obs_date, obs_sat_name=obs_sat_name $
  ,obs_lat_earth=obs_lat_earth, obs_lon_earth= obs_lon_earth, obs_alt_earth=obs_alt_earth $
  ,target_body_name=target_body_name,target_frame_name=target_frame_name

  ;; Input kernel pool
  STRLEN = 50
  ;; 撮像条件と衛星姿勢から決まるものについて ;;
  ;; North vector :: fkが存在するならSpice からもらってくる
  NA_deg = 90. ;+31.07-0.

  ;; Date
  if n_elements(obs_date) eq 0 then begin
    print,"Date is not defined."
    stop
  endif
  utctim = obs_date; ex '2003-04-14T22:10:00'

  ;; UTC => Et
  CSPICE_STR2ET, utctim, et  ;; lsk\naif0010.tls
  ;; Lunar 3 radii
  cspice_bodvrd, target_body_name, 'RADII', 3, radii

  ;; Lunar position seen from a satellite in IAU_EARTH frame with light-time and path corrections at et
  CSPICE_SPKPOS, target_body_name , et, 'ITRF93', 'LT+S',obs_sat_name, lunar_pos, ltime_LT
  tmp_L_distance = sqrt(total(lunar_pos^2.))



  ;; Solar position seen from IAU_Moon frame with light-time and path corrections
  ;; at "et-ltime_LE"
  CSPICE_SPKPOS, 'Sun' , (et-ltime_LT), target_frame_name, 'LT+S',target_body_name, Sun_pos_lunar, ltime_SL
  SL_distance = sqrt(total(sun_pos_lunar^2.))
  AU = 149597870700d0/1000d0 ;; 1 AU (km)
  SL_distance_au = SL_distance/AU

  ;;!! 太陽光入射 ;;
  CSPICE_SUBSLR, "Intercept: ellipsoid",target_body_name,et,target_frame_name, "LT+S", $
    "Earth", subpoint_solar, subepc, subsrfvec_solar
  CSPICE_RECLAT, subpoint_solar, radii, sub_solar_lon, sub_solar_lat

  if n_elements(obs_lon_earth) eq 0 then begin
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; 地球フレーム上での観測者位置 ;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    ;; 地球フレーム上でのTerra位置 ;; 近いのと光時補正の影響が良くわからないので "NONE"
    ;CSPICE_SPKPOS, obs_sat_name , et, 'IAU_EARTH', 'NONE','EARTH', Terra_pos, ltime
    CSPICE_SPKPOS, obs_sat_name , et, 'ITRF93', 'NONE','EARTH', Terra_pos, ltime
    ;; Sub Terra lon lat on Earth
    CSPICE_RECLAT, Terra_pos, ET_ditance, sub_terra_lon_e, sub_terra_lat_e

    ;; 月-衛星間距離
    ;; 時刻 et にTerraに届く光を放つ月の位置が欲しい => "XLT+S", 放射時刻は et-ltime_LT (20120723) => 勘違い? (20120927)
    CSPICE_SPKPOS, target_body_name , et-ltime_LT, target_frame_name, 'LT+S',obs_sat_name, $
      Terra_pos_lunar, ltime
    LT_distance=sqrt(total(Terra_pos_lunar^2.))

    ;; 月面上衛星直下点
    CSPICE_SUBPNT, "Intercept: ellipsoid",target_body_name,et,target_frame_name, "LT+S",obs_sat_name,$
      subpoint_terra, subepc, subsrfvec_terra
    ;; Sub Terra lon lat on Moon
    CSPICE_RECLAT, subpoint_terra, t_radii, sub_terra_lon, sub_terra_lat

  endif else begin
    ;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; fixed location test ;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;

    target = 'MOON'
    obsctr = 'EARTH'
    obsref = 'ITRF93'
    outref = 'MOON_ME'
    abcorr = 'CN+S'
    refloc = 'OBSERVER'

    ;; Observer location ;;
    body = 'EARTH'
    cspice_bodvrd, body, 'RADII', 3, radii
    re   = radii[0]
    rp   = radii[2]
    flat = ( re - rp ) / re

    cspice_pgrrec, body, obs_lon_earth, obs_lat_earth, obs_alt_earth, re, flat, obs_pos

    cspice_spkcpo, target, et, outref, refloc, $
      abcorr, obs_pos, obsctr,         $
      obsref, state,  ltime

    ;; StateにはObserverからみた月方向が入っているので月から見たObserver方向にする ;;
    sub_obs_moon_vec = -state[0:2]
    LT_distance = sqrt(total(sub_obs_moon_vec^2.))

    ;; Sub Terra lon lat on Moon
    CSPICE_RECLAT, sub_obs_moon_vec, t_radii, sub_terra_lon, sub_terra_lat

    ;print,"Observer: "
    ;print,sub_terra_lon*180./!dpi
    ;print,sub_terra_lat*180./!dpi

    ;; 月中心から直下点方向は月から見た観測者方向で置き換え ;;
    subpoint_terra = sub_obs_moon_vec
    ;stop
  endelse

  ;; 月面上の地球直下点
  ;; Sub-Earth lat lon ;; pck\pck0010.tls
  ;; observer = EARTHに光が届くケースで月の位置が欲しい
  ;; 時間の矛盾あり? (20120718) => 無いようす (20120723) => 勘違いしていた? (20120818)
  ;; 時刻etに地球が見る月=>et-ltime_LEに光を放つ月 => XLSの役割を勘違いしていた。 (20120927)
  ;; ** 地球の位置は計算に影響しない, 参考まで **
  CSPICE_SUBPNT, "Intercept: ellipsoid",target_body_name,et,target_frame_name, "LT+S","Earth",$
    subpoint_earth, subepc, subsrfvec_earth
  CSPICE_RECLAT, subpoint_earth, radii, sub_earth_lon, sub_earth_lat

  ;; 月-地球間距離
  CSPICE_SPKPOS, target_body_name , et-ltime_LT, target_frame_name, 'LT+S',"Earth", $
    Earth_pos_lunar, ltime
  LE_distance=sqrt(total(Earth_pos_lunar^2.))

  ;; 単位ベクトル化 ;;
  u_subpoint_solar = subpoint_solar/sqrt(total(subpoint_solar^2.))
  u_subpoint_earth = subpoint_earth/sqrt(total(subpoint_earth^2.))
  u_subpoint_terra = subpoint_terra/sqrt(total(subpoint_terra^2.))

  ;; 位相角 (太陽-月-地球) ;;
  tmp1 = acos(total(u_subpoint_earth*u_subpoint_solar))*180./!dpi

  ;; 位相角 (太陽-月-観測者) ;;
  tmp2 = acos(total(u_subpoint_terra*u_subpoint_solar))*180./!dpi

  ;; Scan観測を再現するために後々必要
  phi_deg = 00./3600. ;; 視線中心からのズレ量
  lam_deg = 00. ;; 画像上で視野中心-天体中心を結ぶ線とx軸のなす角

  tags = ['ST_distance','TO_distance','ssl_longitude','ssl_latitude' $
    ,'ssc_longitude', 'ssc_latitude', 'N_azimuth', 'phi', 'lam','phase']
  obs_geo = create_struct(tags $
    ,SL_distance $
    ,LT_distance $
    ,sub_solar_lon*180./!dpi $
    ,sub_solar_lat*180./!dpi $
    ,sub_terra_lon*180./!dpi $
    ,sub_terra_lat*180./!dpi $
    ,NA_deg $
    ,phi_deg $
    ,lam_deg $
    ,tmp2 $
    )

  return, obs_geo
end


;;
;;  係数の読み取り ;;
;;
function get_params_rolo,sel_wav,ifname_ROLO_csv,param_a,param_b,param_d,coeff_c,coeff_p

  param_a = dblarr(4)
  param_b = dblarr(3)
  param_d = dblarr(3)

  ;; 論文掲載パラメータから ;;
  result = read_csv(ifname_ROLO_csv)
  wav = result.field01
  a0 = result.field02
  a1 = result.field03
  a2 = result.field04
  a3 = result.field05

  b1 = result.field06
  b2 = result.field07
  b3 = result.field08

  d1 = result.field09
  d2 = result.field10
  d3 = result.field11

  min_res = min(abs(wav-sel_wav),min_pos)

  param_a[0] = a0[min_pos]
  param_a[1] = a1[min_pos]
  param_a[2] = a2[min_pos]
  param_a[3] = a3[min_pos]

  param_b[0] = b1[min_pos]
  param_b[1] = b2[min_pos]
  param_b[2] = b3[min_pos]

  param_d[0] = d1[min_pos]
  param_d[1] = d2[min_pos]
  param_d[2] = d3[min_pos]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Coefficients for all bands ;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; 元論文でのc1 は c2に、c2はc1に、c3はc4に、c4はc3になるのがよい, based on Yokota and Stone's discussion
  coeff_c = dblarr(4)
  coeff_c[0] = -0.00134252
  coeff_c[1] = 0.000341152
  coeff_c[3] = 0.000959062
  coeff_c[2] = 0.000662286

  ;; 下記は勘違い? ;;
  ;coeff_c[0] = 0.000341152
  ;coeff_c[1] = -0.00134252
  ;coeff_c[2] = 0.000959062
  ;coeff_c[3] =0.000662286

  coeff_p = dblarr(4)
  coeff_p[0] = 4.06054
  coeff_p[1] = 12.8802
  coeff_p[2] = -30.5858
  coeff_p[3] = 16.7498

  error = 0l
  return,error
end

;;
;;
;;
;;
;; Kieffer and Stone, 2005: Equation(10) ;;
;;
function calculate_rolo_reflectance $
  ,param_a,param_b,param_d,coeff_c,coeff_p $
  ,phase_angle_deg,sub_solar_lon $
  ,sub_obs_lat_deg, sub_obs_lon_deg

  term_A = dblarr(9)
  term_A[0] = 0d
  for i=0, 3, 1 do begin
    term_A[0] += param_a[i]*(phase_angle_deg*!dpi/180d)^double(i)
  endfor

  term_A[1] = 0d
  for j=1, 3, 1 do begin
    term_A[1] += param_b[j-1]*sub_solar_lon^double(2*j-1l)
  endfor

  term_A[2] = coeff_c[0]*sub_obs_lat_deg
  term_A[3] = coeff_c[1]*sub_obs_lon_deg
  term_A[4] = coeff_c[2]*sub_solar_lon*sub_obs_lat_deg
  term_A[5] = coeff_c[3]*sub_solar_lon*sub_obs_lon_deg

  term_A[6] = param_d[0]*exp(-phase_angle_deg / coeff_p[0])
  term_A[7] = param_d[1]*exp(-phase_angle_deg / coeff_p[1])
  term_A[8] = param_d[2]*cos((phase_angle_deg - coeff_p[2])/coeff_p[3])
  ;stop
  ;print,exp(total(term_A))
  result = exp(total(term_A))

  return, result
end


;;
;; apolloデータのフィッティング ;;
;;
function generate_rolo_spec,ifname_ROLO_csv, ifname_soil, ifname_rock $
                           ,phase_angle_deg, sub_solar_lon, sub_obs_lat_deg, sub_obs_lon_deg,out_wav

  ;;
  ;; 源泉データ読み込み (ハードコーディング)
  ;;
;  ifldname_spec = 'C:\work\DATA\LunarData\Apollo_soil\'
;  ifname_soil = ifldname_spec + '62231_data.txt'
;  ifname_rock = ifldname_spec + '67445_data.csv'
;
;  ROLO_dir = 'C:\work\ENVI_IDL\LunaCal\Data\'
;  ifname_ROLO_csv = ROLO_dir + 'ROLO_params.csv'

  ;; Soil
  result = read_ascii(ifname_soil,data_Start=1)
  wav_soil = result.field1[0,*]
  spec_soil = result.field1[1,*]

  ;; Rock
  result_rock = read_csv(ifname_rock,n_table_header=1)
  wav_rock = result_rock.field1[*]
  spec_rock = result_rock.field2[*]

  ;; 波長位置を合わせる ;;
  spec_rock_rev = interpol(spec_rock,wav_rock*1000,wav_soil)

  ;; 両者をマージ, 割合はKieffer & Stone 2005による ;;
  spec_merge = 0.95 * spec_soil + 0.05 * spec_rock_rev

  ;; 有効値域に限定 ;;
  spec_pos = where(spec_rock_rev gt 0.2)

  spec_merge = spec_merge[spec_pos]
  spec_rock_rev = spec_rock_rev[spec_pos]
  spec_soil = spec_soil[spec_pos]
  wav_soil = wav_soil[spec_pos]

  ;; ROLO 係数の読み取り ;;
  result = read_csv(ifname_ROLO_csv, n_header_table=1)
  sel_wav = result.field01

  Disk_reflectance_rolo = dblarr(32)
  Apollo_reflectance = dblarr(32)
  for i=0, 31, 1 do begin
    error = get_params_rolo(sel_wav[i],ifname_ROLO_csv,param_a,param_b,param_d,coeff_c,coeff_p)

    ;; Rolo brightness ;;
    Disk_reflectance_rolo[i] = calculate_rolo_reflectance(param_a,param_b,param_d,coeff_c,coeff_p $
      ,phase_angle_deg, sub_solar_lon, sub_obs_lat_deg, sub_obs_lon_deg)

    ;print, sel_wav[i], Disk_reflectance_rolo[i]

    min_res = min(abs(wav_soil-sel_wav[i]),min_pos)
    Apollo_reflectance[i] = spec_merge[min_pos]

  endfor

  ;;
  ;; 線形フィッティング ;;
  ;; L = Σ (rolo_32 - (a + b * wav)*apollo_spec)^2
  ;; dL/da = 0 => ΣAS*(rolo_32 - (a + b * wav)*AS) = Σ-a*AS^2.-b*wav*AS^2.+AS*rolo_32) = 0
  ;; dL/db = 0 => ΣAS*wav*(rolo_32 - (a + b * wav)*AS) = Σ-a*wav*AS^2.-b*wav^2.*AS^2.+AS*wav*rolo_32) = 0
  ;;
  ;; まとめると[[ΣAS^2, Σwav*AS^2],[Σwav*AS^2., Σwav^2.*AS^2.]] [a,b] = [ΣAS*rolo_32, ΣAS*wav*rolo_32]
  ;; => [a,b] = [[ΣAS^2, Σwav*AS^2],[Σwav*AS^2., Σwav^2.*AS^2.]]^(-1) [ΣAS*rolo_32, ΣAS*wav*rolo_32]
  ;;
  sel_pos = where(sel_wav gt 301)

  array = dblarr(2,2)
  array[0,0] = total(Apollo_reflectance[sel_pos]^2)
  array[1,0] = total(Apollo_reflectance[sel_pos]^2*sel_wav[sel_pos])
  array[0,1] = array[1,0]
  array[1,1] = total(Apollo_reflectance[sel_pos]^2*sel_wav[sel_pos]^2.)

  vec = dblarr(2)
  vec[0] = total(Apollo_reflectance[sel_pos]*Disk_reflectance_rolo[sel_pos])
  vec[1] = total(Apollo_reflectance[sel_pos]*Disk_reflectance_rolo[sel_pos]*sel_wav[sel_pos])

  ;;;;;;;;;;;;;;;;;;;;;;;;
  ;; 線形フィットの結果 ;;
  ;;;;;;;;;;;;;;;;;;;;;;;;
  result = invert(array)##vec
  tmp_fitted_reflectance_ROLO = (result[0]+result[1]*wav_soil)*spec_merge
  print,"Fitted a, b: ",result[0],result[1]

  ;;
  ;; logの方で合わせる? ;;
  ;;
;  tmp_as = alog(Apollo_reflectance)
;  tmp_rl = alog(Disk_reflectance_rolo)
;  array = dblarr(2,2)
;  array[0,0] = total(tmp_as[sel_pos]^2)
;  array[1,0] = total(tmp_as[sel_pos]^2*sel_wav[sel_pos])
;  array[0,1] = array[1,0]
;  array[1,1] = total(tmp_as[sel_pos]^2*sel_wav[sel_pos]^2.)
;
;  vec = dblarr(2)
;  vec[0] = total(tmp_as[sel_pos]*tmp_rl[sel_pos])
;  vec[1] = total(tmp_as[sel_pos]*tmp_rl[sel_pos]*sel_wav[sel_pos])
;
;  ;;;;;;;;;;;;;;;;;;;;;;;;
;  ;; 線形フィットの結果 ;;
;  ;;;;;;;;;;;;;;;;;;;;;;;;
;  result = invert(array)##vec
;  tmp_fitted_reflectance_ROLO = exp((result[0]+result[1]*wav_soil)*alog(spec_merge))
;  print,"Fitted a, b: ",result[0],result[1]

  if n_elements(out_wav) eq 0 then begin
    out_wav = wav_soil
    fitted_reflectance_ROLO = tmp_fitted_reflectance_ROLO
  endif else begin
    ;; 線形内挿により指定の波長にあわせる ;;
    fitted_reflectance_ROLO = interpol(tmp_fitted_reflectance_ROLO,wav_soil,out_wav)

  endelse

  return,fitted_reflectance_ROLO
end