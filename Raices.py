import sympy as sy

x = sy.Symbol('x')

def biseccion(fun,ext,TOL,N):

	'''
    Devuelve la aproximacion de la raiz de la funcion en el intervalo [a,b] y el numero de iteraciones
    necesarias en caso de ser posible, o la aproximacion y 0 en caso de exceder las iteraciones maximas:

            Parametros:
                    fun (x): funcion 
                    ext (list): intervalo [a,b]
                    TOL (float): tolerancia de error
                    N (int): numero maximo de iteraciones en caso de fallo

            Returns:
                    pn (float): aproximacion de la raiz del polinomio
                    i (int): numero de iteraciones necesarias en caso exitoso, 0 sino 
    '''

	a = ext[0]
	b = ext[-1]
	i=1
	FA = f.subs(x,a)
	while i<=N:
		p = a+(b-a)/2
		FP = f.subs(x,p)
		if FP==0 or abs(b-a)/2<TOL:
			return (p,i)
		i+=1
		if FA*FP>0: 
			a=p
			FA=FP
		else:  b=p
	return (p,0)

def PuntoFijo(fun,p0,TOL,N):

	'''
    Devuelve la aproximacion de la raiz de la funcion en el intervalo [a,b] y el numero de iteraciones
    necesarias en caso de ser posible, o la aproximacion y 0 en caso de exceder las iteraciones maximas:

            Parametros:
                    fun (x): funcion 
                    p0 (float): aproximacion inicial
                    TOL (float): tolerancia de error
                    N (int): numero maximo de iteraciones en caso de fallo

            Returns:
                    pn (float): aproximacion de la raiz del polinomio
                    i (int): numero de iteraciones necesarias en caso exitoso, 0 sino 
    '''

	i=1
	while i<=N:
		p = fun.subs(x,p0)
		if abs(p-p0)<TOL:
			return (p,i)
		p0=p	
		i+=1
	return(p,0)

def NewtonRaphson(fun,p0,TOL,N):

	'''
    Devuelve la aproximacion de la raiz de la funcion en el intervalo [a,b] y el numero de iteraciones
    necesarias en caso de ser posible, o la aproximacion y 0 en caso de exceder las iteraciones maximas:

            Parametros:
                    fun (x): funcion 
                    p0 (float): aproximacion inicial
                    TOL (float): tolerancia de error
                    N (int): numero maximo de iteraciones en caso de fallo

            Returns:
                    pn (float): aproximacion de la raiz del polinomio
                    i (int): numero de iteraciones necesarias en caso exitoso, 0 sino 
    '''

	i=1
	while i<=N:
		p = p0 - fun.subs(x,p0)/fun.diff(x).subs(x,p0)
		if abs(p-p0)<TOL:
			return (float(p),i)
		p0=p	
		i+=1
	return(p,0)

def NewtonRaphsonMod(fun,p0,TOL,N):

	'''
    Devuelve la aproximacion de la raiz de la funcion en el intervalo [a,b] y el numero de iteraciones
    necesarias en caso de ser posible, o la aproximacion y 0 en caso de exceder las iteraciones maximas:

            Parametros:
                    fun (x): funcion 
                    p0 (float): aproximacion inicial
                    TOL (float): tolerancia de error
                    N (int): numero maximo de iteraciones en caso de fallo

            Returns:
                    pn (float): aproximacion de la raiz del polinomio
                    i (int): numero de iteraciones necesarias en caso exitoso, 0 sino 
    '''

	i=1
	while i<=N:
		p = p0 - fun.subs(x,p0)*fun.diff(x).subs(x,p0)/(fun.diff(x).subs(x,p0)**2-fun.subs(x,p0)*fun.diff(x,2).subs(x,p0))
		if abs(p-p0)<TOL:
			return (float(p),i)
		p0=p	
		i+=1
	return(p,0)

def Secante(fun,p0,p1,TOL,N):

	'''
    Devuelve la aproximacion de la raiz de la funcion en el intervalo [a,b] y el numero de iteraciones
    necesarias en caso de ser posible, o la aproximacion y 0 en caso de exceder las iteraciones maximas:

            Parametros:
                    fun (x): funcion 
                    p0,p1 (float): aproximaciones iniciales
                    TOL (float): tolerancia de error
                    N (int): numero maximo de iteraciones en caso de fallo

            Returns:
                    pn (float): aproximacion de la raiz del polinomio
                    i (int): numero de iteraciones necesarias en caso exitoso, 0 sino 
    '''

	i=2
	q0 = fun.subs(x,p0)
	q1 = fun.subs(x,p1)
	while i<=N:
		p = p1-q1*(p1-p0)/(q1-q0)
		if abs(p-p1)<TOL:
			return (float(p),i)
		p0 = p1
		q0 = q1
		p1 = p
		q1 = fun.subs(x,p)	
		i+=1
	return(p,0)
