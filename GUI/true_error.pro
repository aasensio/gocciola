pro true, n_lev, n_atm, n_itr
	 a=dblarr(n_lev,n_atm,n_itr)
	 openr,2,'../fort.12'
	 readf,2,a
	 close,2
	 
	 dpop = dblarr(n_itr-1)
	 for i = 0, n_itr-2 do begin
	 	  dpop(i) = max(abs((a(*,*,i)-a(*,*,n_itr-1))) / abs(a(*,*,n_itr-1)))
	 endfor
	 stop
end
