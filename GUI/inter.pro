pro inter
	rmin=1.0006484e+16
	rmax=4.6283938e+17
	r=findgen(100)/99.*(rmax-rmin)+rmin
	 
   x=findgen(100)/10.+1.
	y=alog(x)/max(alog(x))
	r = y*(rmax-rmin)+rmin
;stop
	a=''
	openr,2,'../Atmos/hco+.atmos'
	readf,2,a
	readf,2,a
	readf,2,a
	readf,2,a
	f= fltarr(8,50)
	readf,2,f
	close,2
	
	T=interpol(f(1,*),f(0,*),r)
	nh2= interpol(f(2,*),f(0,*),r)
	n=interpol(f(3,*),f(0,*),r)
	nel=interpol(f(4,*),f(0,*),r)
	vmic=interpol(f(5,*),f(0,*),r)
	Tdust=interpol(f(6,*),f(0,*),r)
	vmac=interpol(f(7,*),f(0,*),r)
	vmac(where(abs(vmac) lt 1.d-5)) = 0.d0	
	print, transpose([[r],[T],[nh2],[n],[nel],[vmic],[Tdust],[vmac]])
	
	openw,2,'../Atmos/hco+interp.atmos',width=132
	printf,2,a
	printf,2,a
	printf,2,a
	printf,2, 100

	printf, 2, transpose([[r],[T],[nh2],[n],[nel],[vmic],[Tdust],[vmac]])
	close,2
	

end
