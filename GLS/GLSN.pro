PRO GLSN,time,ydata,dy,f,p,fbeg=fbeg, fend=fend,  pbeg=pbeg, pend=pend, np=np, ls=ls, nostat=nostat, noplot=noplot, ofac=ofac
;;;2021_11_22 corrected by Jiang, Zhaohao Chen
;;; THis version, we normalize the flux and flux error by dividing it with the maximal value of flux



twopi = 2.0*!DPI
IF NOT KEYWORD_SET(ofac) THEN ofac = 10.   ; set default value

ymax=max(ydata)
ydata=ydata/ymax
dy=dy/ymax
; prepare input
t = time - MIN(time)  &  nt = N_ELEMENTS(t) ; shift the time serie
tbase = MAX(t)                       ; time base
w = DOUBLE(1.+0.*ydata)              ; initialise unit weights
w = 1./dy^2  ; use weights if given
w = w / TOTAL(w)                     ; normalize weights, now TOTAL(w)=1
Y  = TOTAL(w*ydata)    ; Eq.(7), Y is also the (weighted) average of ydata
wy = ydata-Y           ; shifting data to zero mean
;Y  = TOTAL(w*wy)       ; Y should be now zero
YY = TOTAL(w*wy^2)     ; Eq.(10), variance for the zero mean data
wy = w*wy              ; data with attached weights
delta_f = 1./tbase/ofac     ; frequency sampling depends on the time span
fmin    = delta_f           ; default for start frequency
delta_t = tbase/(nt-1)      ; mean time sampling
fmax    = 0.5/delta_t       ; mean Nyquist frequency as default for end frequency

   
; check if user has specified other range
IF KEYWORD_SET(pend) THEN fmin=1./pend
IF KEYWORD_SET(fbeg) THEN fmin=fbeg
fbeg = fmin
IF KEYWORD_SET(pbeg) THEN fmax=1./pbeg
IF KEYWORD_SET(fend) THEN fmax=fend
fend = fmax

; check limits
IF (fend LT fbeg) THEN BEGIN
   PRINT, 'WARNING: fend < fbeg; swapping both values!'
   fmin=fend & fend=fbeg & fbeg=fmin
ENDIF
nout = LONG((fend-fbeg)/delta_f)   ; size of frequency grid
; general infomation
;IF NOT KEYWORD_SET(nostat) THEN BEGIN
  ; PRINT,'  number of input points:     ', nt
   ;PRINT,'  weighted mean of dataset:   ', Y
  ; PRINT,'  weighted rms of dataset:    ', SQRT(YY)
   ;PRINT,'  time base:                  ', tbase
  ; PRINT,'  number of frequency points: ', nout
;ENDIF
f = FINDGEN(nout)*delta_f+fbeg            ; set up frequency grid
P = DOUBLE(0*f) & A=P & B=P & off=P       ; initialize power grid

; calculate the periodogram
FOR k=0L,nout-1 DO BEGIN
   omega = twopi*f[k]    ; angular frequency
   cosx = COS(omega*t)   ; to call the cos-function only one time
   sinx = SIN(omega*t)
   C = TOTAL(w*cosx)     ; Eq.(8)
   S = TOTAL(w*sinx)     ; Eq.(9)
   IF KEYWORD_SET(ls) THEN BEGIN ; set C=S=0 for the classical Lomb-Scargle
      C=0 & S=0
   ENDIF
   YC = TOTAL(wy*cosx)   ;   & YC=YC-Y*C     ; Eq.(11) simplifies, since Y should be 0 (mean was subtracted)
   YS = TOTAL(wy*sinx)   ;   & YS=YS-Y*S     ; Eq.(12)   -"-
   CC = TOTAL(w *cosx^2)                     ; Eq.(13)
   SS = 1.-CC                & SS=SS-S*S     ; Eq.(14)
   CC = CC-C*C                               ; Eq.(13)
   CS = TOTAL(w *cosx*sinx)  & CS=CS-C*S     ; Eq.(15)
   D  = CC*SS-CS^2                           ; Eq.(6)

   A[k] = (YC*SS-YS*CS) / D
   B[k] = (YS*CC-YC*CS) / D
   off[k] = -A[k]*C-B[k]*S
   p[k]=(SS*YC^2/D+CC*YS^2/D-2.*CS*YC*YS/D)/YY  ; Eq.(5)
ENDFOR

   pmax = MAX(p,i) ; returns also the index i of the maximum power
   rms = SQRT(YY*(1-P[i]))
   ; get the curvature in the power peak by fitting a parabola
   ;yh = p[i-1:i+1]-p[i] & xh=(f[i-1:i+1]-f[i])^2  ; shift the parabola origin to power peak
   ;aa = TOTAL(yh*xh)/TOTAL(xh*xh)                 ; calculate the curvature
   Amp = SQRT(A[i]^2+B[i]^2)
   ph  = ATAN(A[i],B[i]) / twopi
   T0  = MIN(time)-ph/f[i]+1./f[i]
   off[i] = off[i]+Y  ; re-add the mean

; statistics:
IF NOT KEYWORD_SET(nostat) THEN BEGIN
;   PRINT,''
  PRINT,'  maximum power p:       ', p[i]
;   PRINT,'  rms of residuals:   ', SQRT(YY*(1-P[i]))
   IF KEYWORD_SET(dy) THEN $
;      PRINT,'mean weighted internal error:     ', SQRT(nt/TOTAL(1./dy^2))
   ;PRINT,'  best sine frequency f: ', f[i],' +/-',SQRT(-2./nt* p[i]/aa*(1-p[i])/p[i])
  ; PRINT,'  best sine period Psin: ', 1./f[i],' +/-',SQRT(-2./nt* p[i]/aa*(1-p[i])/p[i])/f[i]^2
 ;  PRINT,'   amplitude: ', Amp,   ' +/-',SQRT(2./nt)*rm
  ; PRINT,'   phase ph:  ', ph,    ' +/-',SQRT(2./nt)*rms/Amp/twopi
   PRINT,'   phase T0:  ', T0,    ' +/-',SQRT(2./nt)*rms/Amp/twopi/f[i]
   ;PRINT,'   offset:    ', off[i],' +/-',SQRT(1./nt)*rms
;   PRINT, '   Psin_err:', SQRT(6./nt)*2/twopi/tbase*rms/Amp/f[i]^2
;   PRINT, '   f_err:   ', SQRT(6./nt)*2/twopi/tbase*rms/Amp
;   PRINT, '   f_err:   ', SQRT(-2./nt* p[i]/aa*(1-p[i])/p[i]) ; error from peak curvature
ENDIF
END
