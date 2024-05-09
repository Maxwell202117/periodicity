pro GLSP,fbeg=fbeg,fend=fend,pbeg=pbeg,pend=pend,np=np,outarr=outarr,outarrdown=outarrdown,file=file
;Thiis program is used to calculate the significance range and add to the GLS


w=strmid(file,0,9)
set_plot,'ps'
device,filename='GLSN_Sigma_123'+w+'.ps',/color
loadct,39
;;This part is used to simu


f0=file_search(file,count=n0)
f1=file_search('A_*.txt',count=n1)

if ~keyword_set(np) then np=10000
if ~keyword_set(fbeg) then  fbeg=0.0005
if ~keyword_set(fend) then fend=0.5

;print,n1,n0
; simuarr=dblarr(n1,np)
 significance=[0.6826,0.9545,0.9973]
 a_cor=dblarr(1,np)



nline0=file_lines(f0[0])


arr0=dblarr(3,nline0)
  openr,lun0,f0[0],/get_lun
   readf,lun0,arr0
   free_lun,lun0
b=dblarr(nline0)

;GLS,arr0[0,*],arr0[1,*],arr0[2,*],f,p
;h=n_elements(p)
;simuarr=dblarr(n1,h)
;print,h
nline3=file_lines(f1[0])
arr3=dblarr(3,nline3)
   openr,lun3,f1[0],/get_lun
   readf,lun3,arr3
   free_lun,lun3
time=2

for k=0l,nline0-1l do begin
   c=abs(arr0[0,k]-arr3[0,*])
   d=where(c le time,count)

   e=mean(d)
   g=long(e)
   b(k)=g

   endfor
;print,b

GLSN,arr3[0,b],arr3[1,b],arr3[2,b],f,p
h=n_elements(p)
simuarr=dblarr(n1,h)

for j=0l,n1-1l do begin
   nline1=file_lines(f1[j])
   
arr2=dblarr(3,nline1)
   openr,lun2,f1[j],/get_lun
   readf,lun2,arr2
   free_lun,lun2 
   GLSN,arr2[0,b],arr2[1,b],arr2[2,b],f,p

;o=n_elements(f)
;simuarr=dblarr(n1,h)
;print,h,n1,f
simuarr[j,*]=p
 endfor
tsimuarr=transpose([simuarr])
;print,tsimuarr


outarr=dblarr(n_elements(significance),h)
outarrdown=dblarr(n_elements(significance),h)
 for k=0,h-1l do begin
 ; a_cor[0,k]=mean(simuarr[*,k])
  arrtemp=simuarr[*,k]

temp_finite=where(finite(arrtemp),c_fin)
if c_fin eq 0 then begin
outarr[*,k]=!values.d_nan
outarrdown[*,k]=!values.d_nan
endif else begin

arrtemp=arrtemp(temp_finite)
seq=sort(arrtemp)

 arrtemp=arrtemp(seq)
;print,arrtemp
indup=n_elements(arrtemp)*(1+significance)/2-1
outarr[*,k]=arrtemp[indup]
;print,arrtemp[indup]
index_sigdown=n_elements(arrtemp)*(1-significance)/2
outarrdown[*,k]=arrtemp[index_sigdown]
endelse
endfor


;This part is used to print file
 GLSN,arr0[0,*],arr0[1,*],arr0[2,*],f,p
;q=max(p[100:h])
;print,q
;print,n_elements(f)
;print,h
openw,lunout,'GLSN_obs_'+w+'.txt',/get_lun
; for j=0l,n_elements(h)-1l do begin
for j=0l,h-1l do begin

printf,lunout,f[j],p[j],1/f[j],outarr[0:2,j],format='(6f15.7)'
 ;printf,lunout,outarr[0:2,j],format='(3f15.7)'
 endfor
free_lun,lunout

;This part is used to plot
plot,alog10(f),p,xrange=[alog10(fbeg),alog10(0.05)],yrange=[0,1],xtitle='Log(f)',ytitle='Power'
;oplot,f,p
;oplot,!dpi*2./om,outarrdown[0,*],psym=1
;plot,alog10(f),outarr[0,*],psym=2,color=230
oplot,alog10(f),outarr[0,*],psym=2,color=230
;oplot,!dpi*2./om,outarrdown[1,*],psym=2,color=230
oplot,alog10(f),outarr[1,*],psym=5
;oplot,f,outarr[1,*],psym=5
oplot,alog10(f),outarr[2,*],psym=4,color=100
;oplot,f,outarr[2,*],psym=4,color=100
;oplot,!dpi*2./om,outarrdown[2,*],psym=5
device,/close
end

