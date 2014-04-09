;abre_ps, 'hco+2LVG.ps',/todo
!p.multi=[0,2,3]
plot,geometry.r,atom.pop(0,*)/atom.popi(0,*),/xlog,xrange=[7.e15,8.e17],$
         xstyle=1,yrange=[0,5], xtitle='R [cm]', ytitle='n!INLTE!N/n!ILVG!N',ystyle=1
         
plot,geometry.r,atom.pop(1,*)/atom.popi(1,*),/xlog,xrange=[7.e15,8.e17],$
         xstyle=1,yrange=[0,2], xtitle='R [cm]', ytitle='n!INLTE!N/n!ILVG!N',ystyle=1
         
plot,geometry.r,atom.pop(2,*)/atom.popi(2,*),/xlog,xrange=[7.e15,8.e17],$
         xstyle=1,yrange=[0,1.5], xtitle='R [cm]', ytitle='n!INLTE!N/n!ILVG!N',ystyle=1
         
plot,geometry.r,atom.pop(3,*)/atom.popi(3,*),/xlog,xrange=[7.e15,8.e17],$
         xstyle=1,yrange=[0,1.5], xtitle='R [cm]', ytitle='n!INLTE!N/n!ILVG!N',ystyle=1
         
plot,geometry.r,atom.pop(4,*)/atom.popi(4,*),/xlog,xrange=[7.e15,8.e17],$
         xstyle=1,yrange=[0,1.5], xtitle='R [cm]', ytitle='n!INLTE!N/n!ILVG!N',ystyle=1
         
plot,geometry.r,atom.pop(5,*)/atom.popi(5,*),/xlog,xrange=[7.e15,8.e17],$
         xstyle=1,yrange=[0,1], xtitle='R [cm]', ytitle='n!INLTE!N/n!ILVG!N',ystyle=1
         
;plot,geometry.r,atom.pop(6,*)/atom.popi(6,*),/xlog,xrange=[7.e15,8.e17],$
;         xstyle=1,yrange=[0,1.2], xtitle='R [cm]', ytitle='n!INLTE!N/n!ILVG!N',ystyle=1
         
;plot,geometry.r,atom.pop(7,*)/atom.popi(7,*),/xlog,xrange=[7.e15,8.e17],$
;         xstyle=1,yrange=[0,1.4], xtitle='R [cm]', ytitle='n!INLTE!N/n!ILVG!N',ystyle=1
         
;plot,geometry.r,atom.pop(8,*)/atom.popi(8,*),/xlog,xrange=[7.e15,8.e17],$
;         xstyle=1,yrange=[0,1.4], xtitle='R [cm]', ytitle='n!INLTE!N/n!ILVG!N',ystyle=1
         
;plot,geometry.r,atom.pop(9,*)/atom.popi(9,*),/xlog,xrange=[7.e15,8.e17],$
;         xstyle=1,yrange=[0,1.4], xtitle='R [cm]', ytitle='n!INLTE!N/n!ILVG!N',ystyle=1 
         
!p.multi=0
;cierra_ps
