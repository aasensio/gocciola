pro gfunct, x, par, f, pder
	 I0 = par(0)
	 n = par(1)
	 re = par(2)
	 bn = 2.d0*n - 0.324d0
	 f = I0 * exp(-bn*(x/re)^(1.d0/n))
end

pro ajusta
	 analyze, out=a
	 
	 x = reform(a.geometry.r)
	 y = reform(a.atom.pop(3,*)/a.atom.b(3,*))
	 stop
	 	 
	 w = fltarr(n_elements(y))
	 w = replicate(1.d0,n_elements(y)) ;/ y
	 par = [6.d-4,1.4d0,2.0d17]
	 	 
	 yfit = curvefit(x,y,w,par,sigma,function_name='gfunct',/noderivative)
	 
	 plot, x, y, /ylog
	 oplot, x, yfit, psym=4
	 print,par
	 stop
	 
end
