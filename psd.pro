function myfunct,x,a

       model= A(0) * x^A(1) + A(2)

       ;dF(x)/dA(0)= x^A(1)
       ;dF(x)/dA(1)= A(0) * Alog(x)*x^A(1)
       ;dF(x)/dA(2)= 1.0

        bx=A(0)*x^A(1)
        return,[[bx+A(2)],[x^A(1)],[bx*Alog(x)],[1.0]]
 end

Pro PSD,inname=inname
f=inname
n=file_lines(f)
a=dblarr(3,n)
openr,lun,f,/get_lun
readf,lun,a
 free_lun,lun

;readcol,'D:\Kepler Data\PDCSAP_Kepler_w2R\data\kplr006690887-2013131215648_llc.txt',t,c,format='f',count=m
;t=float(t)
;c=float(c)

han = HANNING(n, /DOUBLE)
scargle,a[0,*],a[1,*],om,px
;x=om ; angular frequency of PSD
;x0=om/(2*!DPI);normal frequency  (nu=omega/(2*!DPI))
y0=px   ;the psd-values corresponding to omega
X0=om
x=Alog10(x0)
y=Alog10(y0)
;Gp1=PLOT(x0,y0,symbol='circle',sym_size=1,color='black',linestyle=6,Xtitle='om',Ytitle='psd ')
x1=where(x le -2 )
x2=where(x ge -2 )
;print,y(x1)
a=[0.0017,1.36,0.0]   ;Initial guess
y0fit=linfit(x(x1),y(x1),/double,prob=p2,sigma=sigma)
y1fit=linfit(x(x2),y(x2),/double,prob=p2,sigma=sigma)
y2fit=linfit(x,y,/double,prob=p2,sigma=sigma)
print,y0fit[0],sigma[0],y0fit[1],sigma[1]
print,y1fit[0],sigma[0],y1fit[1],sigma[1]
print,y2fit[0],sigma[0],y2fit[1],sigma[1]

yfit=lmfit(x,y,a,FUNCTION_NAME='myfunct')
yy=A(0) * x^A(1) + A(2)
;Gp1=PLOT(x,yy,symbol='circle',sym_size=1,color='black',linestyle=6)
for i=0,2 do print,i,a(i),format='("A (",i1,")= ",F6.2)'


set_plot,'ps'
device,filename='PSD_'+f+'.ps',/color
loadct,39
plot,x,y
oplot,x,yy,linestyle=1
oplot,x,y0fit[1]*x+y0fit[0],linestyle=2
oplot,x,y1fit[1]*x+y1fit[0],linestyle=2
oplot,x,y2fit[1]*x+y2fit[0],linestyle=2
device,/close
;A (0)=  35.09
;A (1)=  -0.78
;A (2)=  -1.70


openw,lunout,'scargle_'+f+'.txt',/get_lun
 for j=0,n_elements(om)-1 do begin
 printf,lunout,om[j],px[j],!dpi*2./om[j],format='(3f15.7)'
 endfor
free_lun,lunout

End


