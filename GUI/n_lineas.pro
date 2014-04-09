function n_lineas, fichero
 a=' '
 tam=0L
 get_lun, unit
 openr,unit, fichero
 while not eof(unit) do begin
   readf,unit,a
   tam = tam+1
 endwhile
 free_lun, unit
 return, tam

end
