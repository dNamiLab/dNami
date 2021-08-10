#
# With Sympy :
#
import re

varname      = {'rho' : 1,
                  'u' : 2,  
                  'v' : 3, 
                  'et': 4}

varloc       = {'p': 'gamma_m1*rho*(e)',
                'e': '(et-0.5_wp*(u*u+v*v))',
                'T': 'CvInv*(e)'}


import sympy as sym
import numpy as np

x      = sym.Symbol('x')
y      = sym.Symbol('y')
z      = sym.Symbol('z')
t      = sym.Symbol('t')

def apply_dNamivar(expr,varname):

  expr = sym.sympify(expr,evaluate=False)
  for v in varname:
    expr = expr.subs(v,sym.Function(v)(x,y,z,t),simultaneous=True,evaluate=False)
  
  # for b in expr.atoms(sym.Subs):
  #   expr = expr.xreplace({b: b.doit()})  

  return expr   

def dNami_to_sympy(expr,varname):
  # convert derivatives:
  dirvars = ['x','y','z']
  for v in dirvars: 
    expr = expr.replace('{','Derivative(')
    expr = expr.replace('}_1'+v,','+v+')')
    expr = expr.replace('[','Derivative(')
    expr = expr.replace(']_1'+v,','+v+')')
  
  expr = expr.replace('_wp','')

  # convert symbols into functions when needed:
  expr = apply_dNamivar(expr,varname) 
  
  return expr

def sympy2dNami(expr):

  fltcheck = re.compile(r'_wp\d')

  dervSympy = []
  dervdNami = []

  floatSympy = []
  floatdNami = []

  x      = sym.Symbol('x')
  y      = sym.Symbol('y')
  z      = sym.Symbol('z')
  t      = sym.Symbol('t')
  u = sym.Function('u')(x,y,z,t)

  tst = [(x,1),(y,1),(z,1)]

  dummyderv = sym.Derivative(u,z)
  
  crossder = False
  crosscount= 0

  for i in sym.preorder_traversal(expr):
      #convert derivatives  
      if isinstance(i, type(dummyderv)):
        if not crossder:
          crossder   = True 
          dervSympy.append(i)
          if (len(i.args) == 3):
            dervdNami.append('[ {'+str(i.args[0])+' }_1'+str(i.args[1][0])+']_1'+str(i.args[2][0])+' ')
          elif (len(i.args) == 2):
            dervdNami.append('[ '+str(i.args[0])+' ]_1'+str(i.args[1][0])+' ') 
        else:
          dervSympy.append(i)
          if (len(i.args) == 3) or crosscount > 2:
            print("Derivation Order > 2 NRY")
            import sys
            sys.exit()
            dervdNami.append('[ {'+str(i.args[0])+' }_1'+str(i.args[1][0])+']_1'+str(i.args[2][0])+' ')
          elif (len(i.args) == 2):
            dervdNami.append('{ '+str(i.args[0])+' }_1'+str(i.args[1][0])+' ')
        crosscount = crosscount + len(i.args) - 1
        if crosscount > 2:
          print("Derivation Order > 2 NRY",i)
          import sys
          sys.exit()
      if i in tst:
          crossder   = False
          crosscount = crosscount - 1

      # convert floats:    
      if isinstance(i, sym.Float):        
        floatSympy.append(i)
        floatdNami.append(str(i)+'_wp')                       
      # print(floatdNami)  
  expStr = str(expr)

  for dsy,ddn in zip(dervSympy,dervdNami):
    expStr = expStr.replace(str(dsy),ddn)

  fltsyDone = []
  for fltsy,fltdn in zip(floatSympy,floatdNami):
    # Very dirty brute force !!!!
    replaceTrue = True
    for d in range(1,16):                    
      if replaceTrue and fltsy not in fltsyDone:
        expStrOld = expStr   
        expStrtest = expStr.replace(str(fltsy.n(d)),str(fltsy.n(d))+'_wp')
        if fltcheck.search(expStrtest) == None:
          expStr = expStrtest
          replaceTrue = expStrOld == expStr

    fltsyDone.append(fltsy)    

  expStr = expStr.replace('(x, y, z, t)','')
  expStr = expStr.replace('(x, y, t)','')

  return expStr

rho   = sym.Function('rho')(x,y,z,t)
u     = sym.Function('u')(x,y,z,t)
v     = sym.Function('v')(x,y,z,t)
w     = sym.Function('w')(x,y,z,t)
et    = sym.Function('et')(x,y,z,t)
P     = sym.Function('P')(x,y,z,t)


def LODI(dir,dim):

  vel_dir    = {'x':'u','y':'v','z':'w'}

  delta_xdir = {'x':'1.0','y':'0.0','z':'0.0'}
  delta_ydir = {'x':'0.0','y':'1.0','z':'0.0'}
  delta_zdir = {'x':'0.0','y':'0.0','z':'1.0'}

  A   = [[vel_dir[dir]  ,'rho'+'*'+delta_xdir[dir]     ,'rho'+'*'+delta_ydir[dir]     ,'rho'+'*'+delta_zdir[dir]     ,'0.0'                        ],
        ['0.0'          ,vel_dir[dir]                  ,'0.0'                         ,'0.0'                         ,'1.0/rho'+'*'+delta_xdir[dir]],
        ['0.0'          ,'0.0'                         ,vel_dir[dir]                  ,'0.0'                         ,'1.0/rho'+'*'+delta_ydir[dir]],
        ['0.0'          ,'0.0'                         ,'0.0'                         ,vel_dir[dir]                  ,'1.0/rho'+'*'+delta_zdir[dir]],
        ['0.0'          ,'rho*c**2'+'*'+delta_xdir[dir],'rho*c**2'+'*'+delta_ydir[dir],'rho*c**2'+'*'+delta_zdir[dir],vel_dir[dir]                 ]] 

  Phi = [[ '[rho]_1'+dir+' '],
         [ '[ u ]_1'+dir+' '],
         [ '[ v ]_1'+dir+' '],
         [ '[ w ]_1'+dir+' '],
         [ '[ p ]_1'+dir+' ']]    

  if dim <= 2:
     A   = np.delete(A  ,3,1)
     A   = np.delete(A  ,3,0)
     Phi = np.delete(Phi,3,0)
  if dim == 1:
     A   = np.delete(A  ,2,1)
     A   = np.delete(A  ,2,0)
     Phi = np.delete(Phi,2,0)        
  
  A0     = [[dNami_to_sympy(a,{*varname,*varloc}) for a in row] for row in A]
  Asym   = sym.Matrix(A0)   
  PhiSym = sym.Matrix([[dNami_to_sympy(a,{*varname,*varloc}) for a in row] for row in Phi])

  (Sm1,D) = Asym.diagonalize()
  S = Sm1.inv()

  LeftEV = [ v[2] for v in Asym.left_eigenvects()]
  
  EV     = []
  for v in Asym.left_eigenvects():
    for i in range(0,v[1]):
      EV.append(v[0])
  
  
  Pleft = []
  for row in LeftEV:
    for a in row: 
      Pleft.append(a)

  Pleft = sym.Matrix(Pleft)

  LiSym = Pleft*PhiSym

  LiFinal = sym.Matrix([e*l for e,l in zip(EV,LiSym)])

  Li        = []
  Lambda_i  = []

  for ev,l in zip(EV,LiFinal):
    Li.append(sympy2dNami(l))
    Lambda_i.append(sympy2dNami(ev))

  return Lambda_i,Li,Pleft