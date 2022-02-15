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
  
  ;; test 用 or 出力用;;
  ;; Band毎の相対波長透過特性 ;;
  band_response,ch_wav_all=ch_wav_all,contrib_all=contrib_all,/no_plot
  s_ch_wav_all = shift(ch_wav_all,-1)
  pos = where((ch_wav_all-s_ch_wav_all) gt 0)
  ch_wav1 = ch_wav_all[0:pos[0]-1]
  ch_wav2 = ch_wav_all[pos[0]:pos[1]-1]
  ch_wav3 = ch_wav_all[pos[1]:pos[2]-1]
  
  contrib1 = contrib_all[0:pos[0]-1]
  contrib2 = contrib_all[pos[0]:pos[1]-1]
  contrib3 = contrib_all[pos[1]:pos[2]-1]
  
  pos1 = where(contrib1 gt 0.05)
  pos2 = where(contrib2 gt 0.05)
  pos3 = where(contrib3 gt 0.05)
  
  ;; Aster での観測幅での太陽光フラックス
  wav_aster = dindgen(1000)/1000.*(min(wav)-min(ch_wav1[pos1])*1000.) + min(ch_wav1[pos1])*1000.
  ;sol_rad_aster=f_get_radiance_lism(ifname, wav_aster,obs_wav=obs_wav,obs_rad=obs_rad) ; W/m2/nm
  sol_rad_aster=f_get_radiance_integ(obs_wav,obs_rad,wav_aster) ; W/m2/nm
  
  ;; Plot ;;
  loadct,0,/silent
  xs = 450 & xe = 950
  ys = 0.8 & ye = 2.2
  plot,obs_wav,obs_rad,xrange=[xs,xe],yrange=[ys,ye],xstyle=1,ystyle=1,thick=1 $
    ,color=196
    
  loadct,39,/silent
  oplot,wav_aster,2+0*sol_rad_aster,linestyle=0,color=254,thick=6
  oplot,ch_wav1[pos1]*1000, 1.86+0*contrib1[pos1],color=82,linestyle=0,thick=6
  oplot,ch_wav2[pos2]*1000, 1.5 +0*contrib2[pos2],color=128,linestyle=0,thick=6
  oplot,ch_wav3[pos3]*1000, 1.1 +0*contrib3[pos3],color=253,linestyle=0,thick=6
  
  ;; SP modelがカバーする範囲
  oplot,wav,sol_rad,linestyle=0,color=0,thick=3
  
  ;; Band特性
  tmp_ys=0 & tmp_ye = 2.
  plot,ch_wav1*1000, contrib1,color=82,linestyle=1,thick=2,/noerase $
    ,xrange=[xs,xe],yrange=[tmp_ys,tmp_ye],xstyle=5,ystyle=5
  plot,ch_wav2*1000, contrib2,color=128,linestyle=1,thick=2,/noerase $
    ,xrange=[xs,xe],yrange=[tmp_ys,tmp_ye],xstyle=5,ystyle=5
  plot,ch_wav3*1000, contrib3,color=253,linestyle=1,thick=2,/noerase $
    ,xrange=[xs,xe],yrange=[tmp_ys,tmp_ye],xstyle=5,ystyle=5
    
  ;; Plot枠
  plot,obs_wav,obs_rad,xrange=[xs,xe],yrange=[ys,ye],xstyle=1,ystyle=1 $
    ,color=0,/nodata,/noerase,title="Solar irradiance spectral"
    
  ;cd,old_dir
  return
  
end
