;;
;; TLE (Two Line Elements)から軌道カーネルSPKを作る
;; TLEが更新されたら回す
;;

forward_function generate_spk_ROLO
forward_function read_cl_file_for_spk_ROLO

pro read_tle_for_AIST_ROLO, ifname_cl_file, log_output=log_output

  if n_elements(ifname_cl_file) eq 0 then begin
    ifname_cl_file = DIALOG_PICKFILE(/READ, FILTER = '*.txt',TITLE='SELECT Cal file')
    if ifname_cl_file eq '' then return
  endif

  ;;;;;;;;;;;;;;;;;;;;;
  ;; TLE to SPK file ;;
  ;;;;;;;;;;;;;;;;;;;;;

  ifldname = file_dirname(ifname_cl_file,/MARK_DIRECTORY)
  cd, ifldname,current=old

  read_error = read_cl_file_for_spk_ROLO(ifname_cl_file,fldname_kernel,metakr $
                                        ,ifname_tle,ofname_SPK)

  error = generate_spk_ROLO(fldname_kernel,metakr,ifname_tle,ofname_spk $
                      ,log_output=log_output)

  ;; ファイル名の出力 ;;
  full_ofname_spk = file_search("./",file_basename(ofname_spk), /FULLY_QUALIFY_PATH)
  if n_elements(log_output) then AIST_ROLO_log_print,"Genarate SPK... OK: " + full_ofname_spk
  print,"Genarate SPK... OK.", full_ofname_spk

  cd, old
  return  
end


function generate_spk_ROLO,fldname_kernel,metakr,ifname_tle,ofname_spk,log_output=log_output
  
  STRLEN = 50
  ;cd, ifldname_tle, 
  cd, fldname_kernel,current = old_dir
  CSPICE_KCLEAR ;; 初期化
  CSPICE_FURNSH, METAKR
  cd,old_dir
  
  ;; TLE読み込み ;;
  tle_ele = strarr(200000l)
  line = ''
  count=0l
  close,1
  openr,1,ifname_tle
  while(not EOF(1) and count lt 200000l) do begin
    readf,1,line
    if strmid(line,0,1) eq '1' or strmid(line,0,1) eq '2' then begin
      tle_ele[count] = line
      count++
    endif
    ;print,tle_ele[count]
  endwhile
  close,1
  ;; print,count
  
  tle_ele = tle_ele[0:count-1]
  
  epoch_x = dblarr(count/2)
  elems_x = dblarr(count/2*10l)
  o_epoch = -1e10
  j=0l
  for i=0, count-2,2 do begin
    lines = [ tle_ele[i], tle_ele[i+1] ]
    ;cspice_getelm,1999L, lines, epoch, elems
    cspice_getelm,1999L, lines, epoch, elems
    
    ;; SPK must have "epoch" with increasing order
    ;if (i ge 5954 and i le 8170) and (epoch gt 0 and epoch gt o_epoch) then begin
    if epoch gt 0 and epoch gt o_epoch then begin
      ;cspice_et2utc, epoch, 'C', 6, utcstr
      ;print,i,epoch,",  ", utcstr
      
      epoch_x[j]           = epoch
      elems_x[0+j*10:9+j*10] = elems
      j                      = j + 1
    endif
    o_epoch = epoch
    
  endfor
  epoch_x = epoch_x[0:j-1]
  elems_x = elems_x[0:10*j-1]
  
  ;;
  ;; Define the beginning and ending time range
  ;;
  
  ;; Code Sample:
  ;; +/- 12 hours from the first and last epochs
  ;; respectively.
  ;first = epoch_x[0] - 0.5D0*cspice_spd()
  ;last  = epoch_x[j-1] + 0.5D0*cspice_spd()
  
  ;;;;;;;;;;;;;;;;;
  ;; 再現期間を設定 ;;
  ;;;;;;;;;;;;;;;;;
  first = epoch_x[0] - 7D0*cspice_spd() ;; TLEの最初のエポックから何日分さかのぼるか
  last  = epoch_x[j-1] + 365D0*cspice_spd() ;; TLEの最終エポックから何日分作るか
  
  cspice_et2utc, epoch_x[0], 'ISOC', 0, tle_first_utcstr
  cspice_et2utc, epoch_x[j-1], 'ISOC', 0, tle_last_utcstr
  print,"TLE coverage epoch: "
  print, " ",tle_first_utcstr
  print, " ",tle_last_utcstr
  
  cspice_et2utc, first, 'ISOC', 0, pass_first_utcstr
  cspice_et2utc, last, 'ISOC', 0, pass_last_utcstr
  print,"Pass coverage epoch: "
  print, " ", pass_first_utcstr
  print, " ", pass_last_utcstr

  if n_elements(log_output) ne 0 then begin
    AIST_ROLO_log_print,"TLE coverage epoch: "
    AIST_ROLO_log_print, " "+tle_first_utcstr
    AIST_ROLO_log_print, " "+tle_last_utcstr

    AIST_ROLO_log_print,"Pass coverage epoch: "
    AIST_ROLO_log_print, " "+pass_first_utcstr
    AIST_ROLO_log_print, " "+pass_last_utcstr
  endif
  
  ;;
  ;; Create a new SPK file.
  ;;
  
  ;SPK = ofldname_spk + tle_object + '.bsp' ;'UNIFORM.bsp'
  if cspice_exists(ofname_spk) then file_delete, ofname_spk
  cspice_spkopn, ofname_spk, 'TEST_FILE', 1000L, handle
  ;; print,"Handle: ",handle
  
  ;; The constants as listed in geophysical.ker.
  ;;
  CONSTS = [  1.082616D-3,  $
    -2.53881D-6,   $
    -1.65597D-6,   $
    7.43669161D-2, $
    120.D,     $
    78.D,      $
    6378.135D,     $
    1.D ]
    
  ;;
  ;; Add the data for the type 10 segment to the new SPK.
  ;;

  ;ID = -142l ;; -142 = Teraa で UNIFORMを偽装  
  ;cspice_spkw10, handle, ID, 399L, 'j2000', first, last, $
  ;  'Terra', consts, j, elems_x, epoch_x

  ;ID = -5l ;; -142 = Teraa で UNIFORMを偽装
  ID = -999l ;; ダミーコード付与
  cspice_spkw10, handle, ID, 399L, 'j2000', first, last, $
    'TEST', consts, j, elems_x, epoch_x
    
  ;;
  ;; Safely close the SPK.
  ;;
  cspice_spkcls, handle
  CSPICE_KCLEAR

  cd, old_dir

  return,0
end

;;
;; 設定ファイルを読み込む
;;
function read_cl_file_for_spk_ROLO,ifname_cl_file,fldname_kernel,metakr $
                                  ,ifname_tle,ofname_SPK
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
          'KERNEL_file_directory': fldname_kernel = tmp_value
          'METAKR_file_name': METAKR = tmp_value
          'TLE_file_name': ifname_tle = tmp_value
          'SPK_path': ofname_spk = tmp_value
          else: print,"!!",tmp_value
        endcase
      endif
    endelse
    
  endwhile
  close,1

  ;; ファイル名からオブジェクト名作成 ;;
  tle_object = file_basename(ifname_tle, '.tle')

  ;; Outputするファイル名 ;;
  ;ofldname_spk =file_dirname(ifname_tle)+'\'
  ;ofname_SPK = ofldname_spk + tle_object + '.bsp'

  if n_elements(fldname_kernel) eq 0 or n_elements(METAKR) eq 0 then err = -1
  
  return,err
end