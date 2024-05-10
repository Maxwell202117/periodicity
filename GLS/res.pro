pro res4,fbeg=fbeg,fend=fend,pbeg=pbeg,pend=pend,np=np,outarr=outarr,outarrdown=outarrdown,$
file=file




w=strmid(file,5,9)
set_plot,'ps'
device,filename='res'+w+'.ps',/color
loadct,39

f0=file_search(file,count=n0)
f1=file_search('A_*.txt',count=n1)



nline0=file_lines(f0[0])
arr0=dblarr(3,nline0)
  openr,lun0,file,/get_lun
   readf,lun0,arr0
   free_lun,lun0
b=dblarr(nline0)

;openr,lun1,file3,/get_lun
; n1=file_lines(file3)
; arr1=dblarr(3,n1)
; readf,lun1,arr1
;free_lun,lun1


nline3=file_lines(f1[0])
arr3=dblarr(3,nline3)
   openr,lun3,f1[0],/get_lun
   readf,lun3,arr3
   free_lun,lun3
time=1

for k=0l,nline0-1l do begin
   c=abs(arr0[0,k]-arr3[0,*])
   d=where(c le time,count)

   e=mean(d)
   g=long(e)
   b(k)=g
endfor


GLSN,arr3[0,b],arr3[1,b],arr3[2,b],f3,p3

h=n_elements(p3)
xx=dblarr(1,h)
yy=dblarr(1,h)
zz=dblarr(1,h)
for j=0l,n1-1l do begin
   nline1=file_lines(f1[j])
   
arr2=dblarr(3,nline1)
   openr,lun2,f1[j],/get_lun
   readf,lun2,arr2
   free_lun,lun2 
   GLSN,arr2[0,b],arr2[1,b],arr2[2,b],f2,p2
   GLSN,arr2[0,*],arr2[1,*],arr2[2,*],f4,p4


for m=0l,h-1l do begin
xx[*,m]=xx[*,m]+p2[m]-p4[m]
yy[*,m]=yy[*,m]+p2[m]
zz[*,m]=zz[*,m]+p4[m]
endfor

endfor
;print,yy
;print,p3
print,n1
openw,lunout,'1'+w+'.txt',/get_lun,/append
 for j=0,n_elements(f2)-1 do begin
 printf,lunout,f3[j],xx[j]/n1,format='(f15.5,f15.7)'
endfor

plot,alog10(f3),xx/n1,xrange=[-4,-2],yrange=[-0.5,0.5],linestyle=0,color=2
oplot,alog10(f3),zz/n1,linestyle=1,color=230
oplot,alog10(f3),yy/n1,linestyle=2,color=100
device,/close

free_lun,lunout
end
