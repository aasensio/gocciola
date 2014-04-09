pro draw, fich,titulo
	b = ''
	openr,2,fich
	for i = 0, 21 do begin
		readf, 2, b
	endfor
	readf,2,n
	energy = fltarr(4,n)
	for i = 0, 4 do begin
		readf, 2, b
	endfor
	readf,2,n_lines
	for i = 0, 4 do begin
		readf, 2, b
	endfor
	readf,2,energy
	for i = 0, 1 do begin
		readf, 2, b
	endfor
	lines = fltarr(3,n_lines)
	readf,2,lines
	close,2
	!p.multi = 0

	emax = 1.1 * max(energy(0,*))
	plot,[0,0],[0,0],/nodata,xrange=[0,1],yrange=[-1.,emax],ystyle=1,xstyle=1,$
		XTICKFORMAT='(A1)',ytitle='E [cm!e-1!n]',title=titulo
	for i = 0, n-1 do begin
		fromx = 0.05 + i/(n-1.d0)*0.5
		tox = 0.5 + i/(n-1.d0)*0.5
		oplot,[fromx,tox],[energy(0,i),energy(0,i)]
	endfor

	aleat1 = randomu(3123L,n_lines)
	aleat2 = randomu(3124L,n_lines)
	for i = 0, n_lines-1 do begin
		up = fix(lines(0,i))
		low = fix(lines(1,i))
		fromx1 = 0.05 + up/(n-1.d0)*0.5
		fromx2 = 0.5 + up/(n-1.d0)*0.5
		fromx = 0.5*(fromx1+fromx2) + 0.1 * (1.-2.*aleat1(i))
		tox1 = 0.05 + low/(n-1.d0)*0.5
		tox2 = 0.5 + low/(n-1.d0)*0.5
		tox = 0.5*(tox1+tox2) + 0.1 * (1.-2.*aleat2(i))
		arrow, fromx, energy(0,up-1), tox, energy(0,low-1),/data

	endfor

	stop
end
