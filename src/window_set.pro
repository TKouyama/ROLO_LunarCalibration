;;
;; 初期設定＋windowを任意の数出す
;;

pro window_set,w_number,xs=xs,ys=ys,xpos=xpos,ypos=ypos,ct=ct
  if n_elements(w_number) eq 0 then w_number = 1
  
  ;; Delete window
  device, window_state = window_state
  w_pos = where(window_state ne 0)
  if w_pos[0] ne -1 then begin
    for i=w_number, n_elements(w_pos)-1, 1 do begin
      wdelete,w_pos[i]
    endfor
  endif

  if n_elements(xs) eq 0 then xs = 640 ;512
  if n_elements(ys) eq 0 then ys = 640 ;512

  display_size = get_screen_size(resolution = resolution)
  if n_elements(xpos) eq 0 then xpos = display_size[0]-xs-50 ;512
  if n_elements(ypos) eq 0 then ypos = 0 ;512

  ;xs = 640
  ;ys = 640

  ;;
  ;; Windowを表示
  ;;
  for i=w_number-1, 0,-1 do begin
    ;window,xs=xs,ys=ys,xpos=1200,ypos=50*i,i
    window,xs=xs,ys=ys,xpos=xpos,ypos=50*i,i
  endfor

  ;window,xs=640,ys=640,xpos=850,ypos=0,0
  device, decomposed=0

  ;!p.position = [0.247215,0.247215,0.890134,0.890134]
  ;!p.position = [0.15,0.15,0.95,0.95]
  !p.position = [0.15,0.15,0.925,0.925]

  ;; Mitsuyama's
  device,retain=2,decomposed=0
  if n_elements(ct) eq 0 then begin
    loadct,39
  endif else begin
    loadct,ct
  endelse

  !p.background = 255
  !p.color = 0
  !p.font = -1
  !p.charsize =1.25
  !p.charthick =1.
  !p.thick =1;2.5
  !p.linestyle = 0

end