function calcula_tau, r, opa
	n = n_elements(r)
	tau = fltarr(n)
	tau(n-1) = 0.d0
	for i = n-2, 0, -1 do begin
		tau(i) = tau(i+1) + 0.5d0*(opa(i+1)+opa(i)) * abs(r(i+1)-r(i))
	endfor
	return, tau
end