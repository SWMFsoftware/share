;  Copyright (C) 2002 Regents of the University of Michigan, 
;  portions used with permission 
;  For more information, see http://csem.engin.umich.edu/tools/swmf
;===========================================================================
function funcdef, xx, w, func

; Originally developed for the Versatile Advection Code by G. Toth (1996-99).
; Rewritten for the BATSRUS code by G. Toth (2000).
;
; This is the function called by ".r animate" and ".r plotfunc".
;
; "xx"        array contains the "ndim" components of the coordinates.
; "w"         array contains the "nw" variables.
; "func"      string describes the function to be returned.
;
; The "func" string is interpreted by the following rules:
;
;   If the first character is '-' it is stripped off and the function
;       is multiplied by -1 before returned.
;
;   A single number between 0 and nw-1 returns the variable indexed by
;       that number, e.g. '2' in 3D returns w(*,*,*,2). 
;
;   Variable names in the "variables" array mean the appropriate variable.
;
;   Function names defined in the functiondef array or listed in the
;   case statement below are calculated and returned.
;
;   Expressions formed from the 
;   1. standard variables:  rho, ux, uy, uz, uu, u, p, bx, by, bz, bb, b
;   2. standard coordinates: x, y, z, r
;   3. scalar parameters: gamma, gammae, rbody, c0, mi, me, ...
;   4. standard IDL functions and operators
;   5. {names of coordinates}, {names of variables}, {names of equation parameters}
;      Examples: {phi}, {bx1}, {eta} will be replaced with
;                xx(*,*,1), w(*,*,10), eqpar(2)
;   6. {names of functions defined in the functiondef array}
;      Examples: {Ti}, {Mfast}
;
; Examples for valid strings: 
;   '3', 'rho', '{Mfast}+ux', '-T', 'rho*ux', '{bx1}^2+{by1}^2'...
;
; One can use "funcdef" interactively too, e.g.
; 
; ekin =funcdef(x,w,'(ux^2+uy^2)*rho')
; cfast=funcdef(x1,w1,'cfast')
;
;===========================================================================

  common debug_param & on_error, onerror

  ;; Define various functions of the basic MHD variables
  ;; The functions names are evaluated in lower case
  functiondef = $
     strlowcase(transpose([ $
     ['mxB'    , 'rho*ux+(bb*ux-(ux*bx+uy*by+uz*bz)*bx)/c0^2'], $ ; Boris momenta
     ['myB'    , 'rho*uy+(bb*uy-(ux*bx+uy*by+uz*bz)*by)/c0^2'], $
     ['mzB'    , 'rho*uz+(bb*uz-(ux*bx+uy*by+uz*bz)*bz)/c0^2'], $
     ['mx'       , 'rho*ux'                                  ], $ ; momenta
     ['my'       , 'rho*uy'                                  ], $
     ['mz'       , 'rho*uz'                                  ], $
     ['uH'       , 'uH0*sqrt({jx}^2+{jy}^2+{jz}^2)/rho'      ], $ ; Hall velocity
     ['uHx'      , 'uH0*{jx}/rho'                            ], $
     ['uHy'      , 'uH0*{jy}/rho'                            ], $ 
     ['uHz'      , 'uH0*{jz}/rho'                            ], $ 
     ['uex'      , 'ux-uH0*{jx}/rho'                         ], $ ; electron velocity
     ['uey'      , 'uy-uH0*{jy}/rho'                         ], $ 
     ['uez'      , 'uz-uH0*{jz}/rho'                         ], $ 
     ['ue'       , 'sqrt({uex}^2+{uey}^2+{uez}^2)'           ], $
     ['ur'       , '(x*ux+y*uy+z*uz)/r'                      ], $ ; radial u
     ['uxrot'    , 'ux+y*xSI*omegaSunSI/uSI'                 ], $ ; rotational velocity in HGR
     ['uyrot'    , 'uy-x*xSI*omegaSunSI/uSI'                 ], $
     ['uphi'     , '(uy*x-ux*y)/r'                           ], $                       ; uphi
     ['ulon'     , '-sin(Lon)*ux+cos(Lon)*uy'                ], $ ; ulon
     ['ulat'     , '-sin(Lat)*(cos(Lon)*ux+sin(Lon)*uy)+cos(Lat)*uz'], $ ; ulat
     ['ulonrot'  , '-sin(Lon)*{uxrot}+cos(Lon)*{uyrot}'      ], $ ; ulon witout rotation
     ['ulatrot'  , '-sin(Lat)*(cos(Lon)*{uxrot}+sin(Lon)*{uyrot})+cos(Lat)*{uzrot}'], $ ; ulatrot
     ['Br'       , '(x*bx+y*by+z*bz)/r'                      ], $ ; Br
     ['Bt'       , 'sin(Lat)*(cos(Lon)*bx+sin(Lon)*by)-cos(Lat)*bz'], $ ; Btheta
     ['Bp'       , '-sin(Lon)*bx+cos(Lon)*by'                ], $ ; Bphi
     ['Blon'     , '{Bp}'                                    ], $ ; Blon
     ['Blat'     , '-{Bt}'                                   ], $ ; Blat
     ['B1r'      , '(x*{b1x}+y*{b1y}+z*{b1z})/r'             ], $ ; B1r
     ['B1lon'    , '-sin(Lon)*{b1x}+cos(Lon)*{b1y}'          ], $ ; B1lon
     ['B1lat'    , 'cos(Lat)*{b1z}-sin(Lat)*(cos(Lon)*{b1x}+sin(Lon)*{b1y})'], $ ; B1lat
     ['Bphi'     , '(by*x-bx*y)/r'                           ], $ ; Bphi
     ['SignB'    , '{br}/(abs({br})>1e-30)'                  ], $ ; sign(Br)
     ['B0x'      , 'bx-{b1x}'                                ], $ ; B0x
     ['B0y'      , 'by-{b1y}'                                ], $ ; B0y
     ['B0z'      , 'bz-{b1z}'                                ], $ ; B0z
     ['B0'       , 'sqrt({B0x}^2+{B0y}^2+{B0z}^2)'           ], $ ; B0
     ['B0r'      , '(x*{B0x}+y*{B0y}+z*{B0z})/r'             ], $ ; radial B0
     ['B1'       , 'sqrt({b1x}^2+{b1y}^2+{b1z}^2)'           ], $ ; B1
     ['jr'       , '(x*{jx}+y*{jy}+z*{jz})/r'                ], $ ; radial current
     ['jlat'     , 'cos(Lat)*{jz}-sin(Lat)*(cos(Lon)*{jx}+sin(Lon)*{jy})'], $ ; Jlat
     ['jlon'     , '-sin(Lon)*{jx}+cos(Lon)*{jy}'            ], $ ; Jlon
     ['j'        , 'sqrt({jx}^2+{jy}^2+{jz}^2)'              ], $ ; current density
     ['j1x'      , '{jx}-{j0x}'                                ], $ ; j1x
     ['j1y'      , '{jy}-{j0y}'                                ], $ ; j1y
     ['j1z'      , '{jz}-{j0z}'                                ], $ ; j1z
     ['j1'       , 'sqrt({j1x}^2+{j1y}^2+{j1z}^2)'           ], $ ; j1
     ['jxBx'     , '{jy}*bz-{jz}*by'                         ], $ ; Lorentz force x
     ['jxBy'     , '{jz}*bx-{jx}*bz'                         ], $ ; Lorentz force y
     ['jxBz'     , '{jx}*by-{jy}*bx'                         ], $ ; Lorentz force z
     ['jxB'      , 'sqrt({jxBx}^2+{jxBy}^2+{jxBz}^2)'        ], $ ; JxB magnitude
     ['j0xB0x'     , '{j0y}*{b0z}-{j0z}*{b0y}'                   ], $ ; j0xb0 x
     ['j0xB0y'     , '{j0z}*{b0x}-{j0x}*{b0z}'                   ], $ ; j0xb0 y
     ['j0xB0z'     , '{j0x}*{b0y}-{j0y}*{b0x}'                   ], $ ; j0xb0 z
     ['j0xB0'      , 'sqrt({j0xB0x}^2+{j0xB0y}^2+{j0xB0z}^2)'], $ ; J0xB0 magnitude
     ['j0xB1x'     , '{j0y}*{b1z}-{j0z}*{b1y}'                   ], $ ; j0xb1 x
     ['j0xB1y'     , '{j0z}*{b1x}-{j0x}*{b1z}'                   ], $ ; j0xb1 y
     ['j0xB1z'     , '{j0x}*{b1y}-{j0y}*{b1x}'                   ], $ ; j0xb1 z
     ['j0xB1'      , 'sqrt({j0xB1x}^2+{j0xB1y}^2+{j0xB1z}^2)'], $ ; J0xB1 magnitude
     ['j1xB0x'     , '{j1y}*{b0z}-{j1z}*{b0y}'                   ], $ ; j1xb0 x
     ['j1xB0y'     , '{j1z}*{b0x}-{j1x}*{b0z}'                   ], $ ; j1xb0 y
     ['j1xB0z'     , '{j1x}*{b0y}-{j1y}*{b0x}'                   ], $ ; j1xb0 z
     ['j1xB0'      , 'sqrt({j1xB0x}^2+{j1xB0y}^2+{j1xB0z}^2)'], $ ; J1xB0 magnitude
     ['j1xB1x'     , '{j1y}*{b1z}-{j1z}*{b1y}'                   ], $ ; j1xb1 x
     ['j1xB1y'     , '{j1z}*{b1x}-{j1x}*{b1z}'                   ], $ ; j1xb1 y
     ['j1xB1z'     , '{j1x}*{b1y}-{j1y}*{b1x}'                   ], $ ; j1xb1 z
     ['j1xB1'      , 'sqrt({j1xB1x}^2+{j1xB1y}^2+{j1xB1z}^2)'], $ ; J1xB1 magnitude
     ['jxBr'     , '({jxBx}*x+{jxBy}*y+{jxBz}*z)/r'          ], $ ; JxB in r direction
     ['divbxy'   , 'div(bx,by,x,y)'                          ], $ ; div(B) in 2D
     ['divb1xy'  , 'div({b1x},{b1y},x,y)'                    ], $ ; div(B1) in 2D
     ['Ex'       , 'by*uz-uy*bz'                             ], $ ; electric field
     ['Ey'       , 'bz*ux-uz*bx'                             ], $
     ['Ez'       , 'bx*uy-ux*by'                             ], $
     ['e'        , 'p/(gamma-1)+0.5*(rho*uu + bb)'           ], $ ; energy density
     ['pbeta'    , '2*mu0*p/bb'                              ], $ ; plasma beta
     ['s'        , 'p/rho^gamma'                             ], $ ; entropy
     ['se'       , '{pe}/rho^gammae'                         ], $ ; electron entropy
     ['spar'     , '{ppar}*bb/rho^3'                         ], $ ; parallel entropy
     ['sperp'    , '{pperp}/b/rho'                           ], $ ; perp. entropy
     ['Ti'       , 'ti0*p/rho'                               ], $ ; ion temperature [K]
     ['Te'       , 'ti0*{pe}/rho'                            ], $ ; electron temp. [K]
     ['calfvenx' , 'bx/sqrt(rho*mu0A)'                       ], $ ; Alfven velocity
     ['calfveny' , 'by/sqrt(rho*mu0A)'                       ], $
     ['calfvenz' , 'bz/sqrt(rho*mu0A)'                       ], $
     ['calfven'  , 'b /sqrt(rho*mu0A)'                       ], $
     ['Malfvenx' , 'ux/bx*sqrt(rho*mu0A)'                    ], $ ; Alfven Mach number
     ['Malfveny' , 'uy/by*sqrt(rho*mu0A)'                    ], $
     ['Malfvenz' , 'uz/bz*sqrt(rho*mu0A)'                    ], $
     ['Malfven'  , 'u /b *sqrt(rho*mu0A)'                    ], $
     ['csound'   , 'sqrt(gs*p/rho)'                          ], $ ; ion sound speed
     ['csounde'  , 'sqrt(gs*pe/rho*mi/me)'                   ], $ ; electron sound speed
     ['mach'     , 'u /sqrt(gs*p/rho)'                       ], $ ; Mach number
     ['machx'    , 'ux/sqrt(gs*p/rho)'                       ], $
     ['machy'    , 'uy/sqrt(gs*p/rho)'                       ], $
     ['machz'    , 'uz/sqrt(gs*p/rho)'                       ], $
     ['cfast'    , 'sqrt(cc/rho)'                            ], $ ; fast magnetosonic speed
     ['cfastx'   , 'sqrt((cc+sqrt(cc^2-c4*p*bx^2))/2/rho)'   ], $
     ['cfasty'   , 'sqrt((cc+sqrt(cc^2-c4*p*by^2))/2/rho)'   ], $
     ['cfastz'   , 'sqrt((cc+sqrt(cc^2-c4*p*bz^2))/2/rho)'   ], $
     ['cslowx'   , 'sqrt((cc-sqrt(cc^2-c4*p*bx^2))/2/rho)'   ], $ ; slow speed
     ['cslowy'   , 'sqrt((cc-sqrt(cc^2-c4*p*by^2))/2/rho)'   ], $
     ['cslowz'   , 'sqrt((cc-sqrt(cc^2-c4*p*bz^2))/2/rho)'   ], $
     ['Mfast'    , 'sqrt(rho*uu/cc)'                         ], $ ; fast Mach number
     ['Mfastx'   , 'ux/sqrt((cc+sqrt(cc^2-c4*p*bx^2))/2/rho)'], $
     ['Mfasty'   , 'uy/sqrt((cc+sqrt(cc^2-c4*p*by^2))/2/rho)'], $
     ['Mfastz'   , 'uz/sqrt((cc+sqrt(cc^2-c4*p*bz^2))/2/rho)'], $
     ['Mslowx'   , 'ux/sqrt((cc-sqrt(cc^2-c4*p*bx^2))/2/rho)'], $ ; slow Mach number
     ['Mslowy'   , 'uy/sqrt((cc-sqrt(cc^2-c4*p*by^2))/2/rho)'], $
     ['Mslowz'   , 'uz/sqrt((cc-sqrt(cc^2-c4*p*bz^2))/2/rho)'], $
     ['uth'      , 'sqrt(cs0*p/rho)'                         ], $ ; ion thermal speed
     ['uthe'     , 'sqrt(cs0*{pe}/rho*mi/me)'                ], $ ; electron thermal speed
     ['omegapi'  , 'op0*sqrt(rho)'                           ], $ ; ion plasma frequency
     ['omegape'  , 'op0*sqrt(rho*mi/me)'                     ], $ ; electron plasma freq.
     ['omegaci'  , 'oc0*b'                                   ], $ ; ion gyro frequency
     ['omegace'  , 'oc0*b*mi/me'                             ], $ ; electron gyro freq.
     ['rgyro'    , 'rg0*sqrt(p/rho)/(b>1e-30)'               ], $ ; gyro radius  
     ['rgSI'     , 'rg0*sqrt(p/rho)/(b>1e-30)*xSI'           ], $ ; gyro radius in SI
     ['rgyroe'   , 'rg0*sqrt(p/rho*me/mi)/(b>1e-30)'         ], $ ; electron gyro radius  
     ['rgeSI'    , 'rg0*sqrt(p/rho*me/mi)/(b>1e-30)*xSI'     ], $ ; electron gyro radius in SI
     ['dinertial', 'di0/sqrt(rho)'                           ], $ ; inertial length
     ['diSI'     ,' di0/sqrt(rho)*xSI'                       ], $ ; ion inertial length in SI
     ['skindepth',' di0/sqrt(rho*mi/me)'                     ], $ ; electron skin depth
     ['deSI'     ,' di0/sqrt(rho*mi/me)*xSI'                 ], $ ; electron skin depth in SI
     ['ldebye'   , 'ld0/c0*sqrt(p)/rho'                      ], $ ; Debye length
     ['ldSI'     , 'ld0/c0*sqrt(p)/rho*xSI'                  ], $ ; Debye length in SI
     ['n0'       , '{rhos0}/mS0'                             ], $ ; number density of species 0
     ['n1'       , '{rhos1}/mS1'                             ], $ ; number density of species 1
     ['qtot'     , 'qS1*{n1}+qS0*{n0}'                       ], $ ; total charge
     ['dqtot'    , '{qtot}-div({ex},{ey},x,y)*eps0'          ], $ ; error in net charge
     ['ni'       , 'n1'                                      ], $ ; ion number density
     ['ne'       , 'n0'                                      ], $ ; electron number density
     ['npcgs'    , 'Rho/1.6726e-24'                           ], $ ; ion number density in cgs
     ['dqtot1d'  , '{qtot}-diff1({ex},x)*eps0'               ], $ ; error in net charge in 1D
     ['jpx'      , 'qS1*{n1}*{uxs1}+qS0*{n0}*{uxs0}'           ], $ ; jx from particles
     ['jpy'      , 'qS1*{n1}*{uys1}+qS0*{n0}*{uys0}'           ], $ ; jy from particles
     ['jpz'      , 'qS1*{n1}*{uzs1}+qS0*{n0}*{uzs0}'           ], $ ; jz from particles
     ['jex'      , 'qS0*{n0}*{uxs0}'                          ], $ ; jx from electrons
     ['jey'      , 'qS0*{n0}*{uys0}'                          ], $ ; jy from electrons
     ['jez'      , 'qS0*{n0}*{uzs0}'                          ], $ ; jz from electrons
     ['jp'       , 'sqrt({jpx}^2+{jpy}^2+{jpz}^2)'           ], $ ; j from particles
     ['je'       , 'sqrt({jex}^2+{jey}^2+{jez}^2)'           ], $ ; j from electrons 
     ['jppar'    , '({jpx}*{bx}+{jpy}*{by}+{jpz}*{bz})/b'    ], $ ; j parallel to field line
     ['jpperp'   , 'sqrt({jp}^2-{jppar}^2)'                  ], $ ; j perpendicular to field line
     ['jpxbx'    , '{jpy}*{bz}-{jpz}*{by}'                   ], $ ; (j x b)_x
     ['jpxby'    , '{jpz}*{bx}-{jpx}*{bz}'                   ], $ ; (j x b)_y
     ['jpxbz'    , '{jpx}*{by}-{jpy}*{bx}'                   ], $ ; (j x b)_z
     ['p11S0'    , '{pXXS0}*x1^2+{pyyS0}*y1^2+{pzzS0}*z1^2+2*({pxyS0}*x1*y1+{pxzS0}*x1*z1+{pyzS0}*y1*z1)'], $           ;
     ['p22S0'    , '{pXXS0}*x2^2+{pyyS0}*y2^2+{pzzS0}*z2^2+2*({pxyS0}*x2*y2+{pxzS0}*x2*z2+{pyzS0}*y2*z2)'], $           ;
     ['p33S0'    , '{pXXS0}*x3^2+{pyyS0}*y3^2+{pzzS0}*z3^2+2*({pxyS0}*x3*y3+{pxzS0}*x3*z3+{pyzS0}*y3*z3)'], $           ;
     ['p12S0'    , '{pXXS0}*x1*x2+{pyyS0}*y1*y2+{pzzS0}*z1*z2+{pxyS0}*(x1*y2+y1*x2)+{pxzS0}*(x1*z2+z1*x2)+{pyzS0}*(y1*z2+z1*y2)'], $ ;
     ['p13S0'    , '{pXXS0}*x1*x3+{pyyS0}*y1*y3+{pzzS0}*z1*z3+{pxyS0}*(x1*y3+y1*x3)+{pxzS0}*(x1*z3+z1*x3)+{pyzS0}*(y1*z3+z1*y3)'], $ ;
     ['p23S0'    , '{pXXS0}*x2*x3+{pyyS0}*y2*y3+{pzzS0}*z2*z3+{pxyS0}*(x2*y3+y2*x3)+{pxzS0}*(x2*z3+z2*x3)+{pyzS0}*(y2*z3+z2*y3)'], $ ;
     ['p11S1'    , '{pXXS1}*x1^2+{pyyS1}*y1^2+{pzzS1}*z1^2+2*({pxyS1}*x1*y1+{pxzS1}*x1*z1+{pyzS1}*y1*z1)'], $           ;
     ['p22S1'    , '{pXXS1}*x2^2+{pyyS1}*y2^2+{pzzS1}*z2^2+2*({pxyS1}*x2*y2+{pxzS1}*x2*z2+{pyzS1}*y2*z2)'], $           ;
     ['p33S1'    , '{pXXS1}*x3^2+{pyyS1}*y3^2+{pzzS1}*z3^2+2*({pxyS1}*x3*y3+{pxzS1}*x3*z3+{pyzS1}*y3*z3)'], $           ;
     ['p12S1'    , '{pXXS1}*x1*x2+{pyyS1}*y1*y2+{pzzS1}*z1*z2+{pxyS1}*(x1*y2+y1*x2)+{pxzS1}*(x1*z2+z1*x2)+{pyzS1}*(y1*z2+z1*y2)'], $ ;
     ['p13S1'    , '{pXXS1}*x1*x3+{pyyS1}*y1*y3+{pzzS1}*z1*z3+{pxyS1}*(x1*y3+y1*x3)+{pxzS1}*(x1*z3+z1*x3)+{pyzS1}*(y1*z3+z1*y3)'], $ ;
     ['p23S1'    , '{pXXS1}*x2*x3+{pyyS1}*y2*y3+{pzzS1}*z2*z3+{pxyS1}*(x2*y3+y2*x3)+{pxzS1}*(x2*z3+z2*x3)+{pyzS1}*(y2*z3+z2*y3)'],  $ ;
     ['dBn7'     , '{dBnmhd}+{dBnfac}+smooth({dBnhal}+{dBnped},7)'], $ ; dBn smooth 
     ['dBe7'     , '{dBemhd}+{dBefac}+smooth({dBehal}+{dBeped},7)'], $ ; dBn smooth
     ['dBh7'     , 'sqrt({dBn7}^2+{dBe7}^2)'], $  ; dB horizontal with smooth 5
     ['dBh'      , 'sqrt({dBn}^2+{dBe}^2)'],    $  ; dB horizontal
     ['ne14ux'   , '({neuux}*{neurho}+{ne4ux}*{ne4rho})/({neurho}+{ne4rho})'], $
     ['ne14uz'   , '({neuuz}*{neurho}+{ne4uz}*{ne4rho})/({neurho}+{ne4rho})'], $
     ['uxne14'   , '({uxneu}*{rhoneu}+{uxne4}*{rhone4})/({rhoneu}+{rhone4})'], $
     ['uzne14'   , '({uzneu}*{rhoneu}+{uzne4}*{rhone4})/({rhoneu}+{rhone4})']  $
                          ]))

  common file_head
  common phys_units
  common phys_convert
  common phys_const
  common plot_param             ; rcut, vec0

  if n_elements(xx) eq 0 or n_elements(w) eq 0 then begin
     print,'ERROR in funcdef: xx or w are not defined'
     help,xx,w
     retall
  endif

  if n_elements(rcut) eq 0 then rcut = -1

  ;; In 1D xx(n1), in 2D xx(n1,n2,2), in 3D xx(n1,n2,n3,3)
  siz=size(xx)
  ndim=siz(0)-1
  if ndim eq 0 then ndim=1
  n1=siz(1)
  if ndim gt 1 then n2=siz(2)
  if ndim gt 2 then n3=siz(3)

  ;; For 1 variable: w(n1), w(n1*n2), w(n1,n2),  w(n1,n2,n3)
  ;; for more      : w(n1,nw),w(n1,n2,nw),w(n1,n2,n3,nw)
  siz=size(w)
  if siz(0) le ndim then nw=1 else nw=siz(ndim+1)

  ;; Define gamma for the sound speed = sqrt(gs*p/rho) with units
  gs = gamma*cs0

  ;; Define electric permittivity of vacuum
  eps0 = 1/(mu0*c0^2)

  ;; Variable names
  if n_elements(variables) eq 0 then variables=strarr(ndim+nw+neqpar)
  wnames = strlowcase(variables(ndim:ndim+nw-1))

  ;; Check for a negative sign in func
  if strmid(func,0,1) eq '-' then begin 
     f=strmid(func,1,strlen(func)-1)
     sign=-1
  endif else begin
     f=func
     sign=1
  endelse

  ;; Check if f is among the variable names listed in wnames or if it is a number
  for iw=0,nw-1 do $
     if f eq strtrim(string(iw),2) or strlowcase(f) eq wnames(iw) then $
        case ndim of
     1:result = w(*,iw)
     2:result = w(*,*,iw)
     3:result = w(*,*,*,iw)
  endcase

  if n_elements(result) gt 0 and rcut le 0 then return, sign*result

  ;; set radial distance (assuming cartesian coordinates)
  case ndim of
     1: r = abs(xx)
     2: r = sqrt(xx(*,*,0)^2 + xx(*,*,1)^2)
     3: r = sqrt(xx(*,*,*,0)^2 + xx(*,*,*,1)^2 + xx(*,*,*,2)^2)
  endcase

  if n_elements(result) gt 0 then begin
                                ; Return result after cutting out at rcut 
                                ; set value to the minimum of the remaining values
     loc = where(r le rcut, count)
     loc1= where(r gt rcut)
     if count gt 0 then result(loc) = min(result(loc1))
     return, sign*result
  endif

  ;; set the coordinate arrays x, y, z if they occur in the variable names
  x = 0 & y = 0 & z = 0
  for idim = 0, ndim-1 do case ndim of
     1: case variables(idim) of
        'x': x = xx
        'y': y = xx
        'z': z = xx
        else: x = xx
     endcase
     2: case variables(idim) of
        'x': x = xx(*,*,idim)
        'y': y = xx(*,*,idim)
        'z': z = xx(*,*,idim)
        'r'        : Radius = xx(*,*,idim)
        'R'        : Radius = xx(*,*,idim)
        'Radius'   : Radius = xx(*,*,idim)
        'Longitude': Lon = xx(*,*,idim)
        'Lon'      : Lon = xx(*,*,idim)
        'lon'      : Lon = xx(*,*,idim)
        'Latitude' : Lat = xx(*,*,idim)
        'Lat'      : Lat = xx(*,*,idim)
        'lat'      : Lat = xx(*,*,idim)
        else: case idim of
           0: x = xx(*,*,idim)
           1: y = xx(*,*,idim)
        endcase
     end
     3: case variables(idim) of
        'x': x = xx(*,*,*,idim)
        'y': y = xx(*,*,*,idim)
        'z': z = xx(*,*,*,idim)
        'logRadius': Radius = 10^xx(*,*,*,idim)
        'logr'     : Radius = 10^xx(*,*,*,idim)
        'r'        : Radius = xx(*,*,*,idim)
        'R'        : Radius = xx(*,*,*,idim)
        'Radius'   : Radius = xx(*,*,*,idim)
        'Longitude': Lon = xx(*,*,*,idim)   
        'Lon'      : Lon = xx(*,*,*,idim)
        'lon'      : Lon = xx(*,*,*,idim)
        'Latitude' : Lat = xx(*,*,*,idim)   
        'Lat'      : Lat = xx(*,*,*,idim)
        'lat'      : Lat = xx(*,*,*,idim)
        else: case idim of
           0: x = xx(*,*,*,idim)
           1: y = xx(*,*,*,idim)
           2: z = xx(*,*,*,idim)
        endcase
     end
  endcase

  if n_elements(Lon) gt 0 then $
     if max(Lon) gt 2*!pi then Lon = !dtor*Lon

  if n_elements(Lat) gt 0 then $
     if max(abs(Lat)) gt 0.5*!pi then Lat = !dtor*Lat

  if n_elements(Lon) gt 0 and n_elements(Lat) gt 0 then begin
     if n_elements(Radius) eq 0 then Radius = 1.0
     r = Radius
     x = r*cos(Lat)*cos(Lon)
     y = r*cos(Lat)*sin(Lon)
     z = r*sin(Lat)
  end

  ;; Extract primitive variables for calculating MHD type functions

  ;; initialize all the variables as scalars 
  rho=0 & ux=0 & uy=0 & uz=0 & bx=0 & by=0 & bz=0 & p=0 & e=0

  ;; set the variables from the w
  for iw = 0, nw-1 do case ndim of
     1: case wnames(iw) of
        'rho': rho=w(*,iw)
        'ux' : ux=w(*,iw)
        'uy' : uy=w(*,iw)
        'uz' : uz=w(*,iw)
        'bx' : bx=w(*,iw)
        'by' : by=w(*,iw)
        'bz' : bz=w(*,iw)
        'p'  : p=w(*,iw)
        'pth': p=w(*,iw)
        'mx' : ux=w(*,iw)/rho
        'my' : uy=w(*,iw)/rho
        'mz' : uz=w(*,iw)/rho
        'e'  : e=w(*,iw)
        else :
     endcase
     2: case wnames(iw) of
        'rho': rho=w(*,*,iw)
        'ux' : ux=w(*,*,iw)
        'uy' : uy=w(*,*,iw)
        'uz' : uz=w(*,*,iw)
        'bx' : bx=w(*,*,iw)
        'by' : by=w(*,*,iw)
        'bz' : bz=w(*,*,iw)
        'p'  : p=w(*,*,iw)
        'pth': p=w(*,*,iw)
        'mx' : ux=w(*,*,iw)/rho
        'my' : uy=w(*,*,iw)/rho
        'mz' : uz=w(*,*,iw)/rho
        'e'  : e=w(*,*,iw)
        else :
     endcase
     3: case wnames(iw) of
        'rho': rho=w(*,*,*,iw)
        'ux' : ux=w(*,*,*,iw)
        'uy' : uy=w(*,*,*,iw)
        'uz' : uz=w(*,*,*,iw)
        'bx' : bx=w(*,*,*,iw)
        'by' : by=w(*,*,*,iw)
        'bz' : bz=w(*,*,*,iw)
        'p'  : p=w(*,*,*,iw)
        'pth': p=w(*,*,*,iw)
        'mx' : ux=w(*,*,*,iw)/rho
        'my' : uy=w(*,*,*,iw)/rho
        'mz' : uz=w(*,*,*,iw)/rho
        'e'  : e=w(*,*,*,iw)
        else :
     endcase
  endcase

  ;; Extra variables
  uu = ux^2 + uy^2 + uz^2       ; velocity squared, 
  bb = bx^2 + by^2 + bz^2       ; magnetic field squared
  u  = sqrt(uu)                 ; speed
  b  = sqrt(bb)                 ; magnetic field strength

  ;; Change energy into pressure
  if n_elements(p) le 1 and n_elements(e) gt 1 then $
     p = (gamma-1)*(e - 0.5*(rho*uu + bb))

  ;; Calculate gamma*p+bb if needed
  if strpos(f,'fast') ge 0 or strpos(f, 'slow') ge 0 then begin
     c4 = 4*gs/mu0A
     cc = gs*p + bb/mu0A
  end

  ;; Calculate x1, x2, ..., z3 field aligned basis vectors if needed
  if stregex(f,'p[123][123]') ge 0 then begin
     x0 = vec0(0)      & y0 = vec0(1)       & z0 = vec0(2)     ; vec0 from common phys_convert
     x1 = bx/(b>1e-30) & y1 = by/(b>1e-30)  & z1 = bz/(b>1e-30) ; vec1 = (bx, by, bz)/b
     x2 = y0*z1-z0*y1  & y2 = z0*x1-x0*z1   & z2 = x0*y1-y0*x1  ; vec2 = vec0 x vec1 / norm
     z3 = sqrt(x2^2 + y2^2 + z2^2) > 1e-30                      ; use z3 for normalization
     x2 /= z3          & y2 /= z3           & z2 /= z3          ; normalize vec2
     x3 = y1*z2-z1*y2  & y3 = z1*x2-x1*z2   & z3 = x1*y2-y1*x2  ; vec3 = vec1 x vec2
  endif

  ;; Add functions to the basic variable list
  functions = [strlowcase(variables), functiondef(*,0)]

  ;;==== Functions that cannot be expressed from the basic variables
  case f of
     'Btheta'   : result=atan(by,sqrt(bx^2+bz^2))
                                ; sound speed, Mach number
     'jz': case ndim of
        1: result=deriv(x,by)
        2: result=curl(bx,by,x,y)
        3: begin
           print,'Error in funcdef: jz is not implemented for 3D yet'
           retall
        end
     endcase

     ;; Magnetic vector potential A in slab symmetry
     ;; Density of contourlines is proportional to B
     'Ax': begin
        result=dblarr(n1,n2)
                                ; Integrate along the first row
        for i1=1,n1-1 do result(i1,0)=result(i1-1,0) $
           +(bz(i1,0)+bz(i1-1,0))*(y(i1,0)-y(i1-1,0))*0.5 $
           -(by(i1,0)+by(i1-1,0))*(z(i1,0)-z(i1-1,0))*0.5
                                ; Integrate all columns vertically
        for i2=1,n2-1 do result(*,i2)=result(*,i2-1) $
           +(bz(*,i2)+bz(*,i2-1))*(y(*,i2)-y(*,i2-1))*0.5 $
           -(by(*,i2)+by(*,i2-1))*(z(*,i2)-z(*,i2-1))*0.5
     end
     'Ay': begin
        result=dblarr(n1,n2)
                                ; Integrate along the first row
        for i1=1,n1-1 do result(i1,0)=result(i1-1,0) $
           -(bz(i1,0)+bz(i1-1,0))*(x(i1,0)-x(i1-1,0))*0.5 $
           +(bx(i1,0)+bx(i1-1,0))*(z(i1,0)-z(i1-1,0))*0.5
                                ; Integrate all columns vertically
        for i2=1,n2-1 do result(*,i2)=result(*,i2-1) $
           -(bz(*,i2)+bz(*,i2-1))*(x(*,i2)-x(*,i2-1))*0.5 $
           +(bx(*,i2)+bx(*,i2-1))*(z(*,i2)-z(*,i2-1))*0.5
     end
     'Az': begin
        result=dblarr(n1,n2)
                                ; Integrate along the first row
        for i1=1,n1-1 do result(i1,0)=result(i1-1,0) $
           +(by(i1,0)+by(i1-1,0))*(x(i1,0)-x(i1-1,0))*0.5 $
           -(bx(i1,0)+bx(i1-1,0))*(y(i1,0)-y(i1-1,0))*0.5
                                ; Integrate all columns vertically
        for i2=1,n2-1 do result(*,i2)=result(*,i2-1) $
           +(by(*,i2)+by(*,i2-1))*(x(*,i2)-x(*,i2-1))*0.5 $
           -(bx(*,i2)+bx(*,i2-1))*(y(*,i2)-y(*,i2-1))*0.5
     end

                                ; If "f" has not matched any function,
                                ; try evaluating it as an expression
     else:begin
                                ; check if f is the name of a function
        iFunc = where(functiondef(*,0) eq STRLOWCASE(f))
        iFunc = iFunc(0)

        if iFunc ge 0 then f = functiondef(iFunc,1)

        if sign eq -1 then begin
           f='-'+f
           sign=1
        endif

        ;; Create * or *,* or *,*,* for {var} --> w(*,*,iVar) replacements
        case nDim of
           1: Stars = '*,'
           2: Stars = '*,*,'
           3: Stars = '*,*,*,'
        endcase

        ;; replace {name} with appropriate variable
        while strpos(f,'{') ge 0 do begin
           iStart = strpos(f,'{')
           iEnd   = strpos(f,'}')
           Name   = strmid(f,iStart+1,iEnd-iStart-1)
           iVar   = where(functions eq STRLOWCASE(Name))
           iVar   = iVar(0)
           if iVar lt 0 then begin
              print,'Error in funcdef: cannot find variable "{',Name, $
                    '}" among variables'
              print,variables
              retall
           endif
           fStart = strmid(f,0,iStart)
           fEnd   = strmid(f,iEnd+1,strlen(f)-iEnd-1)
           if iVar lt ndim then $
              f = fStart + 'xx(' + Stars + strtrim(iVar,2) + ')' + fEnd $
           else if iVar lt ndim+nw then $
              f = fStart + 'w(' + Stars + strtrim(iVar-nDim,2) + ')' + fEnd $
           else if iVar lt ndim+nw+neqpar then $
              f = fStart + 'eqpar(' + strtrim(iVar-nDim-nW,2) + ')' + fEnd $
           else $
              f = fStart + '(' + functiondef(iVar-nDim-nW-nEqpar,1) + ')' + fEnd

        endwhile

        if not execute('result='+f) then begin
           print,'Error in funcdef: cannot evaluate function=',func
           retall
        endif

     end


  endcase

  if n_elements(result) gt 0 then begin
     
     if rcut gt 0 then begin
                                ; exclude r < rcut
        loc = where(r le rcut, count) 
        loc1= where(r gt rcut)
        if count gt 0 then result(loc) = min(result(loc1))
     endif
     return,sign*result
  endif else begin
     print,'Error in funcdef: function=',func,' was not calculated ?!'
     retall
  endelse

end
