shot=70000
sequence=583
userid='jnv7243'


;Reading from JSP
;------------------------------------------------------------------------
dtype=['ne','te','ti']

for i=0,n_elements(dtype)-1 do begin
  data=dtype(i)+'_jsp'
  cmd="ppfread,shot="+string(shot,format='(I5)')+",dda='jsp',dtype='"+dtype(i)+"',x=x_jsp,t=t_jsp,data="+data+",seq="+string(sequence,format='(I3)')+",ppfuid='"+userid+"'"
  res = execute(cmd)
  print,res
endfor


;Reading from SSP1 (Ar)
;------------------------------------------------------------------------
dtype=['n0','n1','n2','n3','n4','n5','n6','n7','n8','n9']

for i=0,n_elements(dtype)-1 do begin
  data=dtype(i)+'_ssp1'
  cmd="ppfread,shot="+string(shot,format='(I5)')+",dda='ssp1',dtype='"+dtype(i)+"',x=x_ssp1,t=t_ssp1,data="+data+",seq="+string(sequence,format='(I3)')+",ppfuid='"+userid+"'"
  res = execute(cmd)
  print,res
endfor

;Reading from SSP2 (W)
;------------------------------------------------------------------------
dtype=['n0','n1','n2','n3','n4','n5','n6']

for i=0,n_elements(dtype)-1 do begin
  data=dtype(i)+'_ssp2'
  cmd="ppfread,shot="+string(shot,format='(I5)')+",dda='ssp2',dtype='"+dtype(i)+"',x=x_ssp2,t=t_ssp2,data="+data+",seq="+string(sequence,format='(I3)')+",ppfuid='"+userid+"'"
  res = execute(cmd)
  print,res
endfor

save,filename='jetto_output.sav'

end
  
