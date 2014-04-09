; This function reads n non-important lines 
pro lb, unit, n
	a = ''
	for i = 1, n do begin
		readf, unit, a
	endfor
end 


