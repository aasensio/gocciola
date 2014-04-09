tam = 39
;abre_ps, 'dickel.ps',/todo
!p.multi=[0,2,3]
plot,geometry.r,atom.pop(0,0:tam)/total(atom.pop,1),/ylog,/xlog,$
	 xstyle=1,yrange=[0.01,1], xtitle='R [cm]', ytitle='frac. population'
	 
plot,geometry.r,atom.pop(1,0:tam)/total(atom.pop,1),/ylog,/xlog,$
	 xstyle=1,yrange=[0.1,1], xtitle='R [cm]', ytitle='frac. population'
	 
plot,geometry.r,atom.pop(2,0:tam)/total(atom.pop,1),/ylog,/xlog,$
	 xstyle=1,yrange=[0.1,1], xtitle='R [cm]', ytitle='frac. population'
	 
plot,geometry.r,atom.pop(3,0:tam)/total(atom.pop,1),/ylog,/xlog,$
	 xstyle=1,yrange=[0.1,1], xtitle='R [cm]', ytitle='frac. population'
	 
plot,geometry.r,atom.pop(4,0:tam)/total(atom.pop,1),/ylog,/xlog,$
	 xstyle=1,yrange=[0.01,1], xtitle='R [cm]', ytitle='frac. population'
	 
plot,geometry.r,atom.pop(5,0:tam)/total(atom.pop,1),/ylog,/xlog,$
	 xstyle=1,yrange=[1.d-3,1], xtitle='R [cm]', ytitle='frac. population'
	 
plot,geometry.r,atom.pop(6,0:tam)/total(atom.pop,1),/ylog,/xlog,$
	 xstyle=1,yrange=[1.d-8,1], xtitle='R [cm]', ytitle='frac. population'
	 
plot,geometry.r,atom.pop(7,0:tam)/total(atom.pop,1),/ylog,/xlog,$
	 xstyle=1,yrange=[1.d-10,1], xtitle='R [cm]', ytitle='frac. population'
	 
plot,geometry.r,atom.pop(8,0:tam)/total(atom.pop,1),/ylog,/xlog,$
	 xstyle=1,yrange=[1.d-12,1], xtitle='R [cm]', ytitle='frac. population'
	 
plot,geometry.r,atom.pop(9,0:tam)/total(atom.pop,1),/ylog,/xlog,$
	 xstyle=1,yrange=[1.d-15,1], xtitle='R [cm]', ytitle='frac. population'	 	 
	 
!p.multi=0
;cierra_ps

;abre_ps, 'dickel_depart.ps',/todo
!p.multi=[0,2,3]
plot,geometry.r,atom.b(0,0:tam),/xlog,$
	 xstyle=1, xtitle='R [cm]', ytitle='Departure coeff.'
	 
plot,geometry.r,atom.b(1,0:tam),/xlog,$
	 xstyle=1, xtitle='R [cm]', ytitle='Departure coeff.'
	 
plot,geometry.r,atom.b(2,0:tam),/xlog,$
	 xstyle=1, xtitle='R [cm]', ytitle='Departure coeff.'
	 
plot,geometry.r,atom.b(3,0:tam),/xlog,$
	 xstyle=1, xtitle='R [cm]', ytitle='Departure coeff.'
	 
plot,geometry.r,atom.b(4,0:tam),/xlog,$
	 xstyle=1, xtitle='R [cm]', ytitle='Departure coeff.'
	 
plot,geometry.r,atom.b(5,0:tam),/xlog,$
	 xstyle=1, xtitle='R [cm]', ytitle='Departure coeff.'
	 
plot,geometry.r,atom.b(6,0:tam),/xlog,$
	 xstyle=1, xtitle='R [cm]', ytitle='Departure coeff.'
	 
plot,geometry.r,atom.b(7,0:tam),/xlog,$
	 xstyle=1, xtitle='R [cm]', ytitle='Departure coeff.'
	 
plot,geometry.r,atom.b(8,0:tam),/xlog,$
	 xstyle=1, xtitle='R [cm]', ytitle='Departure coeff.'
	 
plot,geometry.r,atom.b(9,0:tam),/xlog,$
	 xstyle=1, xtitle='R [cm]', ytitle='Departure coeff.'	 	 
	 
!p.multi=0
;cierra_ps

;abre_ps, 'dickel_texc.ps',/todo
!p.multi=[0,2,3]
plot, geometry.r, transition.texc(0,0:tam),/xlog,$
	 xstyle=1, yrange=[0,2000], xtitle='R [cm]', ytitle='T!Iexc!N', ystyle=1
xyouts, 1.e19, 201, 'HCO!E+!N J=1-0', /data

plot, geometry.r, transition.texc(1,0:tam),/xlog,$
	 xstyle=1, yrange=[0,50], xtitle='R [cm]', ytitle='T!Iexc!N', ystyle=1
xyouts, 1.e19, 12.5, 'HCO!E+!N J=2-1', /data

plot, geometry.r, transition.texc(2,0:tam),/xlog,$
	 xstyle=1, yrange=[0,50], xtitle='R [cm]', ytitle='T!Iexc!N', ystyle=1
xyouts, 1.e19, 8.5, 'HCO!E+!N J=3-2', /data

plot, geometry.r, transition.texc(3,0:tam),/xlog,$
	 xstyle=1, yrange=[0,50], xtitle='R [cm]', ytitle='T!Iexc!N', ystyle=1
xyouts, 1.e19, 8.5, 'HCO!E+!N J=4-3', /data

plot, geometry.r, transition.texc(4,0:tam),/xlog,$
	 xstyle=1, yrange=[0,50], xtitle='R [cm]', ytitle='T!Iexc!N', ystyle=1
xyouts, 1.e19, 25, 'HCO!E+!N J=5-4', /data

plot, geometry.r, transition.texc(5,0:tam),/xlog,$
	 xstyle=1, yrange=[0,50], xtitle='R [cm]', ytitle='T!Iexc!N', ystyle=1
xyouts, 1.e19, 13, 'HCO!E+!N J=6-5', /data

plot, geometry.r, transition.texc(6,0:tam),/xlog,$
	 xstyle=1, yrange=[0,50], xtitle='R [cm]', ytitle='T!Iexc!N', ystyle=1
xyouts, 1.e19, 11, 'HCO!E+!N J=7-6', /data

plot, geometry.r, transition.texc(7,0:tam),/xlog,$
	 xstyle=1, yrange=[0,50], xtitle='R [cm]', ytitle='T!Iexc!N', ystyle=1
xyouts, 1.e19, 13, 'HCO!E+!N J=8-7', /data

plot, geometry.r, transition.texc(8,0:tam),/xlog,$
	 xstyle=1, yrange=[0,50], xtitle='R [cm]', ytitle='T!Iexc!N', ystyle=1
xyouts, 1.e19, 17, 'HCO!E+!N J=9-8', /data

!p.multi = 0
;cierra_ps

;abre_ps, 'dickel_lvg.ps',/todo
!p.multi=[0,2,3]
plot,geometry.r,atom.pop(0,0:tam)/atom.popi(0,0:tam),/xlog,$
	 xstyle=1,yrange=[0,1.2], xtitle='R [cm]', ytitle='n!INLTE!N/n!ILVG!N'
	 
plot,geometry.r,atom.pop(1,0:tam)/atom.popi(1,0:tam),/xlog,$
	 xstyle=1,yrange=[0,1.2], xtitle='R [cm]', ytitle='n!INLTE!N/n!ILVG!N'
	 
plot,geometry.r,atom.pop(2,0:tam)/atom.popi(2,0:tam),/xlog,$
	 xstyle=1,yrange=[0,1.2], xtitle='R [cm]', ytitle='n!INLTE!N/n!ILVG!N'
	 
plot,geometry.r,atom.pop(3,0:tam)/atom.popi(3,0:tam),/xlog,$
	 xstyle=1,yrange=[0,3], xtitle='R [cm]', ytitle='n!INLTE!N/n!ILVG!N'
	 
plot,geometry.r,atom.pop(4,0:tam)/atom.popi(4,0:tam),/xlog,$
	 xstyle=1,yrange=[0,4], xtitle='R [cm]', ytitle='n!INLTE!N/n!ILVG!N'
	 
plot,geometry.r,atom.pop(5,0:tam)/atom.popi(5,0:tam),/xlog,$
	 xstyle=1,yrange=[0,3], xtitle='R [cm]', ytitle='n!INLTE!N/n!ILVG!N'
	 
plot,geometry.r,atom.pop(6,0:tam)/atom.popi(6,0:tam),/xlog,$
	 xstyle=1,yrange=[0,2], xtitle='R [cm]', ytitle='n!INLTE!N/n!ILVG!N'
	 
plot,geometry.r,atom.pop(7,0:tam)/atom.popi(7,0:tam),/xlog,$
	 xstyle=1,yrange=[0,2], xtitle='R [cm]', ytitle='n!INLTE!N/n!ILVG!N'
	 
plot,geometry.r,atom.pop(8,0:tam)/atom.popi(8,0:tam),/xlog,$
	 xstyle=1,yrange=[0,2], xtitle='R [cm]', ytitle='n!INLTE!N/n!ILVG!N'
	 
plot,geometry.r,atom.pop(9,0:tam)/atom.popi(9,0:tam),/xlog,$
	 xstyle=1,yrange=[0,2], xtitle='R [cm]', ytitle='n!INLTE!N/n!ILVG!N'	 	 
	 
	 
!p.multi=0
;cierra_ps
