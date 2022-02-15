;;
;; Lism標準のデータ
;;

function channel_wavelength_rolo,arr_size=arr_size
  wavelength =  [ 350,  355.1,  405,  412.3,  414.4,  441.6,  465.8,  475 $
    ,486.9, 544,  549.1,  553.8,  665.1,  693.1,  703.6,  745.3 $
    ,763.7,  774.8,  865.3,  872.6,  882,  928.4,  939.3,  942.1 $
    ,  1059.5,  1243.2,  1538.7,  1633.6,  1981.5,  2126.3,  2250.9,  2383.6 $
    ]
  c_size=size(wavelength)
  arr_size=c_size[1]
  return, wavelength
end

;;
;; main
;;
pro get_radiance_ROLO, sol_rad, solfname=solfname, test=test,wav=wav
  ;; Get SP wavelength
  if n_elements(wav) eq 0 then begin
    wav = channel_wavelength_rolo()    
  endif
  
  ;; Read Sol data
  if n_elements(solfname) eq 0 then begin
    ;ifldname = 'C:\work\ENVI_IDL\LunaCal\Data\'
    ;cd, ifldname, current=old_dir
    ifldname = '.\Data\'
    ifname = ifldname+'Whehrli_1985.txt'
  endif else begin
    ifname = solfname
  endelse
  
  count=0l
  tmp_l=''
  tmp_wav = 0.0
  tmp_rad = 0.0
  n_data = 10000l
  obs_rad = fltarr(n_data)
  obs_wav = fltarr(n_data)
  openr,unit,ifname,/get_lun
  readf,unit,tmp_l
  while(not EOF(unit)) do begin
    readf,unit,tmp_wav,tmp_rad
    obs_wav[count] = tmp_wav
    obs_rad[count] = tmp_rad
    count++
  endwhile
  close,unit
  free_lun,unit
  obs_wav = obs_wav[0:count-1]
  obs_rad = obs_rad[0:count-1]
  
  ;; 近い波長の明るさで代表させる ;;
  sol_rad=interpol(obs_rad,obs_wav,wav) ; W/m2/nm
  ;; 台形積分して太陽光スペクトルを得るのはROLOではよくない ;;
  ;sol_rad=f_get_radiance_integ(obs_wav,obs_rad,wav) ; W/m2/nm
  
  if keyword_set(test) ne 1 then begin
    ;cd,old_dir
    return
  endif

  ;; SP modelに対応したチャンネルで太陽光入射量出力 ;;
  ;print,sol_rad
 

  ;cd,old_dir
  return
  
end
