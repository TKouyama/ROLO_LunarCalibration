;;
;; Kieffer and Stone, 2005: Equation(10) ;;
;;
function calculate_rolo_brightness $
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
