;;
;; AIST ROLO widget
;; exeへのbuild方法
;; MAKE_RT, 'AIST_ROLO', 'C:\work\ENVI_IDL\LunaCal\tmp', SAVEFILE='AIST_ROLO.sav',/win64,/vm, /OVERWRITE
;; manufest_rt.txtにicy.dlmを加えた
;;
pro AIST_ROLO_log_print,log_string
  common common_log_output, log_text_box
  top_line = WIDGET_INFO(log_text_box,/TEXT_NUMBER)

  WIDGET_CONTROL, log_text_box, SET_VALUE=log_string,/append,SET_TEXT_TOP_LINE=top_line
  return
end

;;
;; main
;;
pro AIST_ROLO
  window_set,2

  ;; 動作に使うファイル名を共有 ;;
  common common_file_name, ifname_cl_file

  common common_log_output, log_text_box

  common common_current_fld, current_fld

  ;;;;;;;;;;;;;;;;;;;;;;;;
  ;; 計算設定ファイルを指定する ;; ここでは設定をよみこまない
  ;;;;;;;;;;;;;;;;;;;;;;;;
  ifname_cl_file = DIALOG_PICKFILE(/READ, FILTER = '*.txt',TITLE='SELECT Cal file')
  if ifname_cl_file eq '' then return

  ifldname = file_dirname(ifname_cl_file,/MARK_DIRECTORY)
  cd, ifldname,current=current_fld

  ;;;;;;;;;;;;;;
  ;; 外見の整形 ;;
  ;;;;;;;;;;;;;;

  ;;;;;;;;;;;;;;;;;;;;
  ;; Base Frameを作る ;;
  ;;;;;;;;;;;;;;;;;;;;
  base_window = WIDGET_BASE(xsize=512*2,ysize=512+128,xoffset=50,yoffset=50)

  ;;;;;;;;;;;
  ;; 1 列目 ;;
  ;;;;;;;;;;;
  ;; プログラムボタンの配置 ;;
  bottun_a = WIDGET_BUTTON(base_window, VALUE='1. TLE to SPK for SPICE',xoffset=50,yoffset=50,xsize=200,ysize=50)
  bottun_b = WIDGET_BUTTON(base_window, VALUE='2. AIST_ROLO',xoffset=50,yoffset=125,xsize=200,ysize=50)

  ;; Exit ボタン ;;
  exit_button = WIDGET_BUTTON(base_window, VALUE='Exit'   ,xoffset=50,yoffset=425+150,xsize=400,ysize=50)

  ;; 空のText box の配置 ;;
  base_text_box = WIDGET_TEXT(base_window,value='AIST_ROLO programs <TK soft>',uname='base_text_box',xsize=28)
  TLE_to_SPK_text_box  = WIDGET_TEXT(base_window,value='',uname='TLE_to_SPK_text_box',xoffset=275,yoffset=65)
  Obs_date_text_box    = WIDGET_TEXT(base_window,value='',uname='obs_date_text_box',xoffset=275,yoffset=140)

  ;; 空のText box の配置 ;;
  log_text_box = WIDGET_TEXT(base_window,uname='log_text_box',/editable $
                             ,xoffset=512,yoffset=25,xsize=65,ysize=40,/scroll,font="*14")

  ;;;;;;;;;;;;;;;;
  ;; Realize?する ;;
  ;;;;;;;;;;;;;;;;
  WIDGET_CONTROL, /REALIZE, base_window
  
  ;; 各ボタンにイベントをアサイン ;;
  XMANAGER, 'TLE_to_SPK',bottun_a, event_handler='TLE_to_SPK_event',/NO_BLOCK
  XMANAGER, 'AIST_ROLO',bottun_b, event_handler='AIST_ROL_event',/NO_BLOCK
  XMANAGER, 'exit', exit_button , /NO_BLOCK

  return
end

;;
;;
;;
;;
;; TLE to SPK
;;
PRO TLE_to_SPK_event,ev
  common common_file_name, ifname_cl_file
  AIST_ROLO_log_print,"--- Genarate SPK LOG ---"

  text_box = WIDGET_INFO(ev.top,/all_children,find_by_uname='base_text_box')
  WIDGET_CONTROL, text_box, SET_VALUE='Processing...'

  read_tle_for_AIST_ROLO, ifname_cl_file,/log_output
  WIDGET_CONTROL, text_box, SET_VALUE='Done'

  AIST_ROLO_log_print,""

  return
end


;;
;; ROLO_main
;;
PRO AIST_ROL_event, ev
  common common_file_name, ifname_cl_file
  AIST_ROLO_log_print,"--- AIST ROLO LOG ---"


  text_box = WIDGET_INFO(ev.top,/all_children,find_by_uname='base_text_box')
  WIDGET_CONTROL, text_box, SET_VALUE='Processing...'
 
  estimate_ROLO_brignthess_band,ifname_cl_file,/log_output

  box_text = "Done"
  WIDGET_CONTROL, text_box, SET_VALUE=box_text

  AIST_ROLO_log_print,""

  return
END

;;
;; Programを終わらせる
;;
PRO exit_event, ev
  print,"Exit"

  common common_current_fld, current_fld
  ;cd, current_fld

  WIDGET_CONTROL, ev.top,/destroy


  ;; window close
  device, window_state = window_state
  w_pos = where(window_state ne 0)
  w_number = 0
  if w_pos[0] ne -1 then begin
    for i=w_number, n_elements(w_pos)-1, 1 do begin
      wdelete,w_pos[i]
    endfor
  endif

  return
END
