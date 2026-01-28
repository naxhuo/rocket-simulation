import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sin
import sympy as sp
import math


'''
-----------------------------------------------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------------------------------------------
DATOS
-----------------------------------------------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------------------------------------------
'''


'''
-----------------------------------------------------------------------------------------------------------------------------------------------
DATOS COHETES SOLO 1 ETAPA
-----------------------------------------------------------------------------------------------------------------------------------------------
'''

M = 2290000  # kg # Masa total inicial
mc = 130000  # kg #Masa cohete sin combustible
mu_0 = M - mc  # kg # Masa de combustible inicial
T = 168  # s # tiempo que tarda en consumir toda el combustible
t = np.linspace(0, T, 10**3)
u = 2000  # m/s # Velocidad de salidad de los gases


'''
-----------------------------------------------------------------------------------------------------------------------------------------------
VARIABLES QUE NOS SERVIRÁN EN EL PROGRAMA
-----------------------------------------------------------------------------------------------------------------------------------------------
'''

#n_list = np.array([0.25, 0.5, 0.75, 1, 2, 4]) #ctes de consumo que utilizaremos
n_list = np.array([0.25, 0.5, 0.75, 1, 2, 4])  # ctes de consumo que utilizaremos
models = np.array(["singrav", "gravcte", "gravar", "rozcte"]
                  )  # modelos a calcular
methods = ['RK45', 'RK23', 'DOP853', 'Radau',
           'BDF', 'LSODA']  # métodos usados por solve_ivp
k = 10**2  # cte >>1
v_0 = 0  # velocidad inicial
y_0 = 0  # posición inicial
g = 9.80665  # gravedad
cin = np.array([y_0, v_0])  # condiciones iniciales

t_in = np.array([0, T*2])  # intervalo de tiempo
# puntos donde se valorarán los datos
tspan = np.arange(t_in[0], t_in[-1], 0.01)


'''
-----------------------------------------------------------------------------------------------------------------------------------------------
MODELO DE ROZAMIENTO
-----------------------------------------------------------------------------------------------------------------------------------------------
'''

p0 = 1.2266  # Densidad del aire en la superficie
gamma = 1/8243  # Cte para el modelo exponencial del cálculo de la densidad del aire
A = np.pi*(10/2)**2  # Sección del cohete
cd = 0.65  # Cte que depende de la forma del cohete # en este caso con cono con apertura 40º


"""
-----------------------------------------------------------------------------------------------------------------------------------------------
DATOS COHETE VARIAS ETAPAS
-----------------------------------------------------------------------------------------------------------------------------------------------
"""

"""
------------------------------------
Masa de las etapas y del combustible
------------------------------------
"""

mseco_et1 = 137000  # mseco_et1: Masa en seco de la etapa 1
mcom_et1 = 2077000  # mcom_et1: Masa de combustible de la etapa 1

mseco_et2 = 36000  # mseco_et2: Masa en seco de la etapa 2
mcom_et2 = 444000  # mcom_et2: Masa de combustible de la etapa 2

mseco_et3 = 10000  # mseco_et3: Masa en seco de la etapa 3
mcom_et3 = 109000  # mcom_et3: Masa de combustible de la etapa 3

msecolist = [mseco_et1, mseco_et2, mseco_et3]
mulist = [mcom_et1, mcom_et2, mcom_et3]


"""
--------------------------------------------
Tiempo de quema de combustible de las etapas
--------------------------------------------
"""
T1 = 161  # Tiempo de quema de combustible de la primera etapa
T2 = 360  # Tiempo de quema de combustible de la segunda etapa
T3 = 150  # Tiempo de quema de combustible de la tercera etapa

Tlist_etap = [T1, T2, T3]

tspan_etap = np.arange(0, T2+T1+T3+500)
t_in_etap = [0, T1+T2+T3+500]


"""
-------------------------------------------------
Velocidad de expulsión de los gases de las etapas
-------------------------------------------------
"""
u1 = 2580  # Velocidad de expulsión de gases de la primera etapa
u2 = 4130  # Velocidad de expulsión de gases de la segunda etapa
u3 = 4130  # Velocidad de expulsión de gases de la tercera etapa
ulist_etap = [u1, u2, u3]


'''
-----------------------------------------------------------------------------------------------------------------------------------------------
Lista de distintas latitudes en famosas bases espaciales
-----------------------------------------------------------------------------------------------------------------------------------------------
'''

latcape_canaveral = math.radians(28.396837)
latjiuquan = math.radians(40.95778)
latguiana = math.radians(5.05)
lathawkebay = math.radians(-39.26071)

'''
-----------------------------------------------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------------------------------------------
FUNCIONES
-----------------------------------------------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------------------------------------------
'''
'''
-----------------------------------------------------------------------------------------------------------------------------------------------
FUNCIONES DE MASA (Funciones que nos servirán para calcular el consumo y la masa en cada momento)
-----------------------------------------------------------------------------------------------------------------------------------------------
'''


def mu_t(t, n, mu_0, T, k=10**5):

    # Ecuaciones 1.1 del guion (hoja 4)
    '''

    Parameters:
     - t (float,array): tiempos
     - n (float): velocidad de combustión
     - mu_0 (float): masa de combustible inicial
     - T (float): tiempo en el que se consume todo el combustible
     - k (float): cte >> 1
    
    Objective:
     Calcular el combustible restante en el cohete en función del tiempo. 

    Returns:
     Devuelve la transformación de t con el mu correspondiente, teniendo en cuenta que el resto de valores son constantes. 
    '''
    if n <= 0:
        print('Valor de n no válido, n > 0')
        return

    elif n >= 1:
        tau = t/T
        return mu_0*(1-tau**n)

    else:
        tau = (k*t+T)/(k*T+T)
        return (1+1/((k+1)**n-1))*mu_0*(1-tau**n)

def m_t(t, n, mu_0, T, mc, k=10**5):
    '''
    Parameters:
     - t (float,array): tiempos
     - n (float): velocidad de combustión
     - mu_0 (float): masa de combustible inicial
     - T (float): tiempo en el que se consume todo el combustible
     - mc (float): masa en seco del cohete
     - k (float): cte >> 1

    Objective:
     Calcular masa del cohete en el tiempo dado.

    Returns:
     Devuelve la transformación de t con el peso correspondiente, teniendo en cuenta que el resto de valores son constantes.
    '''
    return mu_t(t, n, mu_0, T, k) + mc  # Masa comb. + masa en seco

def dm_t(t, n, mu_0, T, k=10**5):
    '''
    Parameters:
     - t (float,array): tiempos
     - n (float): velocidad de combustión
     - mu_0 (float): masa de combustible inicial 
     - T (float): tiempo en el que se consume todo el combustible 
     - k (float): cte >> 1

    Objective:
     Calcular la variación de masa del cohete en función del tiempo.

    Returns:
     Devuelve la transformación de t con la variación de masa correspondiente, teniendo en cuenta que el resto de valores son constantes.
    '''

    if n <= 0:
        print('Valor de n no válido, n > 0')
        return

    elif n >= 1:
        tau = t/T
        return -mu_0*n/T*tau**(n-1)

    else:
        tau = (k*t+T)/(k*T+T)
        return -(1+(1/(((k+1)**n)-1)))*((mu_0*n)/T)*(k/(k+1))*tau**(n-1)


'''
Funciones para que calcule la masa después del final de la quema
'''

def mu_t_mod(t, n, mu_0, T, k=10**5):

    # Ecuaciones 1.1 del guion (hoja 4)
    '''

    Parameters:
     - t (float,array): tiempos
     - n (float): velocidad de combustión
     - mu_0 (float): masa de combustible inicial
     - T (float): tiempo en el que se consume todo el combustible
     - k (float): cte >> 1
    
    Objective:
     Calcular el combustible restante en el cohete en función del tiempo. 

    Returns:
     Devuelve la transformación de t con el mu correspondiente, teniendo en cuenta que el resto de valores son constantes. 
    '''
    if t <= T:
        if n <= 0:
            print('Valor de n no válido, n > 0')
            return

        elif n >= 1:
            tau = t/T
            return mu_0*(1-tau**n)

        else:
            tau = (k*t+T)/(k*T+T)
            return (1+1/((k+1)**n-1))*mu_0*(1-tau**n)

    else: return 0

def m_t_mod(t, n, mu_0, T, mc, k=10**5):
    '''
    Parameters:
     - t (float,array): tiempos
     - n (float): velocidad de combustión
     - mu_0 (float): masa de combustible inicial
     - T (float): tiempo en el que se consume todo el combustible
     - mc (float): masa en seco del cohete
     - k (float): cte >> 1

    Objective:
     Calcular masa del cohete en el tiempo dado.

    Returns:
     Devuelve la transformación de t con el peso correspondiente, teniendo en cuenta que el resto de valores son constantes.
    '''
    return mu_t_mod(t, n, mu_0, T, k) + mc  # Masa comb. + masa en seco


def dm_t_mod(t, n, mu_0, T, k=10**5):
    '''
    Parameters:
     - t (float,array): tiempos
     - n (float): velocidad de combustión
     - mu_0 (float): masa de combustible inicial 
     - T (float): tiempo en el que se consume todo el combustible 
     - k (float): cte >> 1

    Objective:
     Calcular la variación de masa del cohete en función del tiempo.

    Returns:
     Devuelve la transformación de t con la variación de masa correspondiente
    '''
    if t<= T:
        if n <= 0:
            print('Valor de n no válido, n > 0')
            return

        elif n >= 1:
            tau = t/T
            return -mu_0*n/T*tau**(n-1)

        else:
            tau = (k*t+T)/(k*T+T)
            return -(1+(1/(((k+1)**n)-1)))*((mu_0*n)/T)*(k/(k+1))*tau**(n-1)
    else: return 0

'''
Funciones de masa para el modelo de varias etapas
'''

def m_t_etap(t, n, mu_0s, Ts, mcs, k=10**5):
    '''
    Parameters:
     - t (float,array): tiempos
     - n (float): velocidad de combustión
     - mu_0s (lista): masas de combustible inicial
     - Ts (lista): tiempos en los que se consume todo el combustible de la etapa correspondiente
     - mcs (lista): masas en seco del cohete
     - k (float): cte >> 1

    Objective:
     Calcular masa del cohete en el tiempo dado, para varias etapas

    Returns:
     Devuelve la transformación de t con el peso correspondiente, teniendo en cuenta la etapa en la que se encuentra
    '''

    T1, T2, T3 = Ts
    mc1, mc2, mc3 = mcs
    mu01, mu02, mu03 = mu_0s

    if t <= T1:
        # Masa comb. + masa en seco
        return mu_t(t, n, mu01, T1, k) + np.sum(mcs) + mu02 + mu03
    elif (t > T1) & (t <= T1 + T2):
        return mu_t(t-T1, n, mu02, T2, k) + mc2 + mc3 + mu03
    elif (t > T1+T2) & (t<=T1+T2+T3):
        return mu_t(t-T1-T2, n, mu03, T3, k) + mc3
    else: 
        return mc3


def dm_t_etap(t, n, mu_0s, Ts, k=10**5):
    '''
    Parameters:
     - t (float,array): tiempos
     - n (float): velocidad de combustión
     - mu_0s (list): masas de combustible inicial 
     - Ts (list): tiempos en los que se consume todo el combustible de la etapa correspondiente
     - k (float): cte >> 1

    Objective:
     Calcular la variación de masa del cohete en función del tiempo.

    Returns:
     Devuelve la transformación de t con la variación de masa correspondiente, teniendo en cuenta que el resto de valores son constantes.
    '''

    T1, T2, T3 = Ts
    mu01, mu02, mu03 = mu_0s

    if t <= T1:
        mu_0 = mu01
        T = T1
        tet = t
        if n <= 0:
            print('Valor de n no válido, n > 0')
            return
        elif n >= 1:
            tau = tet/T
            return -mu_0*n/T*tau**(n-1)
        else:
            tau = (k*tet+T)/(k*T+T)
            return -(1+(1/(((k+1)**n)-1)))*((mu_0*n)/T)*(k/(k+1))*tau**(n-1)

    elif (t > T1) & (t <= T1 + T2):
        mu_0 = mu02
        T = T2
        tet = t-T1
        if n <= 0:
            print('Valor de n no válido, n > 0')
            return
        elif n >= 1:
            tau = tet/T
            return -mu_0*n/T*tau**(n-1)

        else:
            tau = (k*tet+T)/(k*T+T)
            return -(1+(1/(((k+1)**n)-1)))*((mu_0*n)/T)*(k/(k+1))*tau**(n-1)
    elif (t > T1+T2) & (t <= T1+T2+T3):
        mu_0 = mu03
        T = T3
        tet = t-T1-T2
        if n <= 0:
            print('Valor de n no válido, n > 0')
            return
        elif n >= 1:
            tau = tet/T
            return -mu_0*n/T*tau**(n-1)
        else:
            tau = (k*tet+T)/(k*T+T)
            return -(1+(1/(((k+1)**n)-1)))*((mu_0*n)/T)*(k/(k+1))*tau**(n-1)

    else:
        return 0


def u_etap(t, Ts = Tlist_etap, us = ulist_etap):

    '''
    Parameters:
     - t(float, array): tiempos
     - Ts(list): tiempos en los que se consume todo el combustible de la etapa correspondiente
     - us (list): velocidad de salida de los gases en cada etapa

    Returns:
     Devuelve la velocidad de salida de los gases dependiendo de la etapa en la que este
    ''' 

    T1,T2,T3 = Ts
    u1,u2,u3 = ulist_etap

    if t <= T1:
        return u1   
    elif (t > T1) & (t <= T1 + T2):
        return u2
    elif (t > T1+T2) & (t <= T1+T2+T3):
        return u3
    else:
        return 0


'''
-----------------------------------------------------------------------------------------------------------------------------------------------
FUNCIÓN PARA LA REPRESENTACIÓN DE LA EVOLUCIÓN DE LA MASA
----------------------------------------------------------------------------------------------------------------------------------------------- 
'''

def RepresMasas(n_list=n_list, mu_0=mu_0, T=T, k=k):

    '''
    Parameters:
     - n_list (lista): cte de consumo que se quieren comparar
     - mu_0 (float): masa de combustible inicial
     - T (float): tiempo en el que se consume todo el combustible
     - k (float): cte >> 1

    Returns:
     Representación de las masas en función del tiempo para todas las ctes de consumo en n_list
    
    '''

    plt.figure(1, figsize=(12, 8), clear=True)
    plt.suptitle('Quema de Combustible')

    # plt.rcParams['text.usetex'] = True

    ''' Evolución de la masa de combustible'''
    ax1 = plt.subplot(1, 2, 1)
    # plt.title(r'$\textbf{\frac{1}{m_0}}$', fontsize=15)
    plt.title(r'mu vs t', fontsize=15)
    plt.xlabel("t/T")
    plt.grid()

    ''' Evolución del gasto de combustible'''
    ax2 = plt.subplot(1, 2, 2)
    # plt.title(r'$\textbf{$\frac{1}{m_0}\frac{d}{dt}m(t)}$', fontsize=15)
    plt.title(r'dm/dt vs t', fontsize=15)
    plt.xlabel("t/T")
    plt.grid()

    for n in n_list:

        mut = mu_t(t, n, mu_0, T, k)
        dmt = dm_t(t, n, mu_0, T, k)

        ax1.plot(t/T, mut/mu_0, label=f'n = {n}')
        ax1.legend()
        ax2.plot(t/T, dmt/mu_0, label=f'n = {n}')
        # ax2.legend()

    plt.show()


'''
-----------------------------------------------------------------------------------------------------------------------------------------------
FUNCIONES PARA EL CÁLCULO DE LA GRAVEDAD
-----------------------------------------------------------------------------------------------------------------------------------------------
'''

def g0lat(lat, g0=9.780327, a=0.0053024, b=0.0000058):
    '''
    Parameters:
     - lat (float): latitud del punto de despegue
     - g0 (float): constante número 1
     - a (float): constante número 2
     - b (float): constante número 3

    Objective:
     Calcular la gravedad en función de la latitud.

    Returns:
     Gravedad (según el modelo dependiente de la latitud del lanzamiento mostrada en los apuntes del trabajo)
    '''

    return g0*(1+a*(np.sin(lat))**2-b*(np.sin(2*lat))**2)


def g0latWGS(lat, g0=9.7803253359, a=0.001931852652, b=0.0066943799901):
    '''
    Parameters:
     - lat (float): latitud del punto de despegue
     - g0 (float): constante número 1
     - a (float): constante número 2
     - b (float): constante número 3

    Objective:
     Calcular la gravedad en función de la latitud, tiene en cuenta la forma elipsoidal de la Tierra (modelo WGS).

    Returns:
     Gravedad (según el modelo dependiente de la latitud del lanzamiento mostrada en los apuntes del trabajo)
    '''

    return g0*((1+a*(np.sin(lat))**2)/np.sqrt(1-(b*(np.sin(lat))**2)))


def g_y(h, Rp=6.3710e6, g0=9.80665):
    '''
    Parameters:
     - h (float): altura a la que se mide la gravedad
     - Rp (int): radio del planeta
     - g0 (float): gravedad en la superficie

    Objective:
     Calcular la gravedad en función de la altura.

    Returns:
     Gravedad (en este caso esta tomada la ecuacion 1.13 (pag 14) con los datos terrestres)
    '''

    return g0/(1+h/Rp)**2


def g_y2(h, lat, Rp=6.3710e6):
    '''
    Parameters:
     - h (float): altura a la que se mide la gravedad
     - lat (float): latitud del punto de despegue
     - Rp (int): radio del planeta

    Objective:
     Calcular la gravedad en función de la altura y de la latitud.

    Returns:
     Gravedad (en este caso esta tomada la ecuacion 1.13 (pag 14) con los datos terrestres y modificada según la latitud)
    '''
    g0 = g0lat(lat)
    return g0/(1+h/Rp)**2


def g_y3(h, lat, Rp=6.3710e6):
    '''
    Parameters:
     - h (float): altura a la que se mide la gravedad
     - lat (float): latitud del punto de despegue
     - Rp (int): radio del planeta

    Objective:
     Calcular la gravedad en función de la altura y de la latitud.

    Returns:
     Gravedad (en este caso esta tomada la ecuacion 1.13 (pag 14) con los datos terrestres y modificada según la latitud)
    '''
    g0 = g0latWGS(lat)
    return g0/(1+h/Rp)**2


'''
-----------------------------------------------------------------------------------------------------------------------------------------------
¡FUNCIÓN PARA LA CREACIÓN DE PLOTS
-----------------------------------------------------------------------------------------------------------------------------------------------   
'''

def CreacionPlot(plot_masa=True,suptitle=''):

    '''
    Parameters:
     - plot_masa (boolean): si quieres o no la representación de la masa frente al tiempo
     - suptitle (string): Titulo de la figura

    Objective:
     Crea la base de un plot para las futuras representaciones

    Returns:
     Los distintos subplots creados (2 si plot_masa False y 3 si True)
    '''

    plt.figure(2, figsize=(16, 12), clear=True)
    plt.suptitle(suptitle)

    if plot_masa:
        ''' Altura respecto al tiempo'''
        ax1 = plt.subplot(2, 2, 1)
        # plt.title(r'$\textbf{\frac{1}{m_0}}$', fontsize=15)
        plt.title(r'Altura', fontsize=15)
        plt.ylabel('y (m)')
        plt.xlabel('t (s)')
        plt.grid()

        ''' Velocidad respecto al tiempo'''
        ax2 = plt.subplot(2, 2, 3)
        # plt.title(r'$\textbf{$\frac{1}{m_0}\frac{d}{dt}m(t)}$', fontsize=15)
        plt.title(r'Velocidad', fontsize=15)
        plt.ylabel('v (m/s)')
        plt.xlabel('t (s)')
        plt.grid()

        ''' Porcentaje de masa respecto al tiempo'''
        ax3 = plt.subplot(2, 2, (2, 4))
        # plt.title(r'$\textbf{$\frac{1}{m_0}\frac{d}{dt}m(t)}$', fontsize=15)
        plt.title(r'Consumo de la masa', fontsize=15)
        plt.ylabel('Fracción de mu restante')
        plt.xlabel('t (s)')
        plt.grid()
        return ax1,ax2,ax3
        
    else:
        ''' Altura respecto al tiempo'''
        ax1 = plt.subplot(2, 1, 1)
        # plt.title(r'$\textbf{\frac{1}{m_0}}$', fontsize=15)
        plt.title(r'Altura', fontsize=15)
        plt.ylabel('y (m)')
        plt.xlabel('t (s)')
        plt.grid()

        ''' Velocidad respecto al tiempo'''
        ax2 = plt.subplot(2, 1, 2)
        # plt.title(r'$\textbf{$\frac{1}{m_0}\frac{d}{dt}m(t)}$', fontsize=15)
        plt.title(r'Velocidad', fontsize=15)
        plt.ylabel('v (m/s)')
        plt.xlabel('t (s)')
        plt.grid()

        return ax1,ax2


'''
-----------------------------------------------------------------------------------------------------------------------------------------------
FUNCIÓN PARA EL CÁLCULO DE TRAYECTORIAS EN UNA ETAPA CON INTEGRALES
-----------------------------------------------------------------------------------------------------------------------------------------------   
'''

def pos_t(t, n, mu_0, T, mc, u, v_0=0, y_0=0, k=10**5, model='SinGravedad'):
    '''
    Parameters:
     - t (float,array): tiempos
     - n (float): velocidad de combustión
     - mu_0 (float): masa de combustible inicial
     - T (float): tiempo en el que se consume todo el combustible
     - mc (float): masa del cohete en seco
     - u (float): velocidad a la que sale el gas del cohete
     - v_0 (float): velocidad inicial del cohete
     - y_0 (float): altura inicial del cohete
     - k (float): cte >> 1
     - model (string): modelo usado. Válidos: "singrav", "gravcte", "gravar"

    Objective:
     Calcular la posición del cohete (altura) y su velocidad (eje y) respecto al tiempo, asi como la masa restante, con el método de integración

    Returns:
     y (array): vector de posiciones en los tiempos donde se pide medir
     v (array): vector de las velocidades
     porcm (array): vector del porcentaje de masa restante 
    '''

    if model == 'singrav':

        mu = mu_t(t, n, mu_0, T, k)
        porcm = mu/mu_0
        v = v_0+u*np.log((mc+mu_0)/(mc+mu))
        y = np.array([])

        for i in range(len(t)):
            ti = t[:i+1]
            vi = v[:i+1]
            integ = sin.simps(vi, ti)
            yi = y_0 + integ
            y = np.append(y, yi)
        return y, v, porcm

    elif model == 'gravcte':
        mu = mu_t(t, n, mu_0, T, k)
        porcm = mu/mu_0
        v = v_0-g*t+u*np.log((mc+mu_0)/(mc+mu))
        y = np.array([])

        for i in range(len(t)):
            ti = t[:i+1]
            vi = v[:i+1]
            integ = sin.simps(vi, ti)
            yi = y_0 + integ
            y = np.append(y, yi)
        return y, v, porcm

    elif model == 'gravar':
        mu = mu_t(t, n, mu_0, T, k)
        porcm = mu/mu_0

        v = np.array([])
        y = np.array([])
        gs = np.array([])
        y = np.append(y, y_0)
        v = np.append(v, v_0)

        for i in range(1, len(t)):
            ti = t[:i+1]
            gi = g_y(y[i-1])
            gs = np.append(gs, gi)
            vins = v_0-gi*t[i-1]+u*np.log((mc+mu_0)/(mc+mu[i]))
            v = np.append(v, vins)
            vi = v[:i+1]
            integ = sin.simps(vi, ti)
            yi = y_0 + integ
            y = np.append(y, yi)

        gi = g_y3(y[-1], latjiuquan)
        gs = np.append(gs, gi)

        return y, v, porcm

    else:
        print('Method not valid')


'''
-----------------------------------------------------------------------------------------------------------------------------------------------
FUNCIONES PARA LA REPRESENTACIÓN DE UNA ETAPA CON INTEGRALES
----------------------------------------------------------------------------------------------------------------------------------------------- 
'''

def IntegRepresN(pos_t, ax1,ax2,ax3, model='singrav', n_list=n_list, mu_0=mu_0, T=T, mc=mc, u=u, v_0=v_0, y_0=y_0, k=k):
    '''
    Parameters:
     - pos_t (function): función especifica para el movimiento en cada modelo
     - ax1,ax2,ax3: subplots creadas anteriormente con CreacionPlot
     - n_list (lista): cte de consumo que se quieren representar
     - + Variables antes descritas
     

    Objective:
        Representacion de posicion, velocidad y masa calculados con integracion para todas las n introducidas 

    Returns:
        Una plot sin representar
     
    '''
   
    for n in n_list:

        y, v, porcm = pos_t(t, n, mu_0, T, mc, u, v_0, y_0, k,model=model)

        ax1.plot(t, y, label=f'n={n}')
        ax2.plot(t, v, label=f'n={n}')
        ax3.plot(t, porcm, label=f'n={n}')
        ax1.legend()

    
def IntegRepresModel(pos_t, ax1,ax2,ax3,models=['singrav'], n=1, mu_0=mu_0, T=T, mc=mc, u=u, v_0=v_0, y_0=y_0, k=k):

    '''
    Parameters:
     - pos_t (function): función especifica para el movimiento en cada modelo
     - ax1,ax2,ax3: subplots creadas anteriormente con CreacionPlot
     - n (int): cte que se quiere aplicar
     - models (list): modelos que se quieren estudiar
     - + Variables antes descritas
     

    Objective:
        Representacion de posicion, velocidad y masa calculados con integracion para cada modelo pasado   

    Returns:
        Una plot sin representar
     
    '''

    for model in models:

        y, v, porcm = pos_t(t, n, mu_0, T, mc, u, v_0, y_0, k,model = model)

        ax1.plot(t, y, label=f'model={model}')
        ax2.plot(t, v, label=f'model={model}')
        ax3.plot(t, porcm, label=f'model={model}')
        ax1.legend()


'''
-----------------------------------------------------------------------------------------------------------------------------------------------
FUNCIONES PARA LA CREACION DE LOS SEDO'S (una y 3 etapas)
-----------------------------------------------------------------------------------------------------------------------------------------------  
'''

def odefunGen(t, Y, modelo='singrav', n=1, mu_0=mu_0, mc=mc,T=T,k=k,lat=latcape_canaveral, A=A, cd =cd, gamma = gamma, u=u):

    '''
    Parameters:
     - t (int): tiempos
     - Y (list): [posición,velocidad]
     - Modelo de fricción: A: sección, cd: cte dependiente de la forma, gammma: cte usada en el modelo exponencial de densidad del aire
     - + Variables antes descritas

    Objective:
        Calcular el SEDO que nos interese, para una etapa solamente

    Returns:
        Devuelve el SEDO de posición y velocidad correspondiente a cada modelo

    '''

    if modelo == "singrav":
        return [Y[1], -u*dm_t_mod(t, n, mu_0, T, k)/m_t_mod(t, n, mu_0, T, mc, k)]
    elif modelo == "gravcte":
        return [Y[1], (-u*dm_t_mod(t, n, mu_0, T, k) - m_t_mod(t, n, mu_0, T, mc, k)*g)/m_t_mod(t, n, mu_0, T, mc, k)]
    elif modelo == "gravar":
        return [Y[1], (-u*dm_t_mod(t, n, mu_0, T, k) - m_t_mod(t, n, mu_0, T, mc, k)*g_y2(Y[0], lat))/m_t_mod(t, n, mu_0, T, mc, k)]
    elif modelo == "rozcte":
        return [Y[1],((-p0*np.exp(-gamma*Y[0]))*0.5*A*cd*Y[1]**2 - u*dm_t_mod(t, n, mu_0, T, k) - m_t_mod(t, n, mu_0, T, mc, k)*g_y3(Y[0], lat=latcape_canaveral))/m_t_mod(t, n, mu_0, T, mc, k)]


def odefunGen_etap(t, Y, modelo='singrav', n=1, mu_0=mulist, mc=msecolist, T=Tlist_etap, k=k, lat=latcape_canaveral, A=A, cd=cd, gamma=gamma, u_etap=u_etap):
    
    '''
    Parameters:
     - t (int): tiempos
     - Y (list): [posición,velocidad]
     - Modelo de fricción: A: sección, cd: cte dependiente de la forma, gammma: cte usada en el modelo exponencial de densidad del aire
     - + Variables antes descritas: nu_0, mc, T en formato lista

    Objective:
        Calcular el SEDO que nos interese, para tres etapas 

    Returns:
        Devuelve el SEDO de posición y velocidad correspondiente a cada modelo

    '''
    if modelo == "singrav":
        return [Y[1], -u_etap(t)*dm_t_etap(t, n, mu_0, T, k)/m_t_etap(t, n, mu_0, T, mc, k)]
    elif modelo == "gravcte":
        return [Y[1], (-u_etap(t)*dm_t_etap(t, n, mu_0, T, k) - m_t_etap(t, n, mu_0, T, mc, k)*g)/m_t_etap(t, n, mu_0, T, mc, k)]
    elif modelo == "gravar":
        return [Y[1], (-u_etap(t)*dm_t_etap(t, n, mu_0, T, k) - m_t_etap(t, n, mu_0, T, mc, k)*g_y2(Y[0], lat))/m_t_etap(t, n, mu_0, T, mc, k)]
    elif modelo == "rozcte":
        return [Y[1], ((-p0*np.exp(-gamma*Y[0]))*0.5*A*cd*Y[1]**2 - u_etap(t)*dm_t_etap(t, n, mu_0, T, k) - m_t_etap(t, n, mu_0, T, mc, k)*g_y3(Y[0], lat=latcape_canaveral))/m_t_etap(t, 1, mu_0, T, mc, k)]


'''
-----------------------------------------------------------------------------------------------------------------------------------------------
FUNCIONES PARA LA REPRESENTACIÓN DE UNA ETAPA CON SEDO'S
-----------------------------------------------------------------------------------------------------------------------------------------------   
'''

def OdeSystemRepres(ax1,ax2,ax3,modelo = 'singrav',t_in=t_in, cin=cin, tspan=tspan, n_list=n_list, method='RK45'):
    '''
    Parameters:
     - ax1,ax2,ax3: subplots creadas anteriormente con CreacionPlot
     - method (string): método que utilizará solve_ivp para hacer los cálculos
     - t_in (lista): intervalo donde se hacen los cálculos
     - cin (lista): condiciones iniciales
     - tspan (array): puntos donde se valoran los distintos sistemas
     - + Variables antes descritas

    Objective:
     Representar la posicion y velocidad del ODE dado para distintos n y un modelo

    Returns:
     Plot sin representar      
    '''

    for n in n_list:

        def odefun(t, Y): return odefunGen(t,Y,modelo=modelo,n=n)

        result = sin.solve_ivp(odefun, t_in, cin, t_eval=tspan, method=method)
        tsol = result.t
        v = result.y[1]
        y = result.y[0]

        porcm = [mu_t_mod(t, n, mu_0, T, k)/mu_0 for t in tspan]        

        ax1.plot(tsol, y, label=f'n={n}')
        ax2.plot(tsol, v, label=f'n={n}')
        if n == n_list[-1]:
            ax2.axvline(x=T, color='k', linestyle='--', label='Fin combustible')
            ax1.axvline(x=T, color='k', linestyle='--', label='Fin combustible')
        ax3.plot(tspan, porcm, label=f'n={n}')
        ax1.legend()


def OdeSystemRepresModels(ax1,ax2,models=['singrav'], n=1, t_in=t_in, cin=cin, tspan=tspan, method='RK45'):
    '''
    Parameters:
     - ax1,ax2,ax3: subplots creadas anteriormente con CreacionPlot 
     - models (list): modelos que se quieren estudiar
     - n (int): cte que se quiere aplicar
     - t_in (lista): intervalo donde se hacen los cálculos
     - cin (lista): condiciones iniciales
     - tspan (array): puntos donde se valoran los distintos sistemas
     - method (string): método que utilizará solve_ivp para hacer los cálculos

    Objective:
     Representar la posicion y velocidad del ODE dado para distintos modelos y un n

    Returns:
     Plot sin representar      
    '''

    for modelo in models:

        def odefun(t, Y): return odefunGen(t, Y, modelo=modelo, n=n)

        result = sin.solve_ivp(odefun, t_in, cin, t_eval=tspan, method=method)
        tsol = result.t
        v = result.y[1]
        y = result.y[0]
        ax1.plot(tsol, y, label=f'model={modelo}')
        ax2.plot(tsol, v, label=f'model={modelo}')
        if modelo == models[-1]:
            ax2.axvline(x=T, color='k', linestyle='--',label='Fin combustible')
            ax1.axvline(x=T, color='k', linestyle='--',label='Fin combustible')
        ax1.legend()


def OdeSystemRepresMethods(ax1,ax2,model = 'singrav',n=1, t_in=t_in, cin=cin, tspan=tspan, methods=['RK45']):
    '''
    Parameters:
     - ax1,ax2,ax3: subplots creadas anteriormente con CreacionPlot 
     - model (string): modelo que se quiere estudiar
     - n (int): cte que se quiere aplicar
     - t_in (lista): intervalo donde se hacen los cálculos
     - cin (lista): condiciones iniciales
     - tspan (array): puntos donde se valoran los distintos sistemas
     - methods (list): método que utilizará solve_ivp para hacer los cálculos

    Objective:
     Representar la posicion y velocidad del ODE dado para distintos métodos, un modelo y un n

    Returns:
     Plot sin representar      
    '''

    for met in methods:

        def odefun(t, Y): return odefunGen(t, Y, modelo=model, n=n)

        result = sin.solve_ivp(odefun, t_in, cin, t_eval=tspan, method=met)
        tsol = result.t
        v = result.y[1]
        y = result.y[0]

        # mu = mu_t(tspan, n, mu_0, T, k)
        # porcm = mu/mu_0

        ax1.plot(tsol, y, label=f'method={met}')
        ax2.plot(tsol, v, label=f'method={met}')
        if met == methods[-1]:
            ax2.axvline(x=T, color='k', linestyle='--', label='Fin combustible')
            ax1.axvline(x=T, color='k', linestyle='--', label='Fin combustible')
        # ax3.plot(tspan, porcm, label=f'method={met}')
        ax1.legend()
       

def ComparacMetodos(mu_0=mu_0, T=T, k=k, cin=cin, suptitle='', methods=['RK45']):

    '''
    Parameters:
     - methods (list): métodos a comparar. Válidos: 'RK45', 'RK23', 'DOP853', 'Radau', 'BDF', 'LSODA'

    Objective:
     Comparar los distintos métodos existentes en solve_ivp con la solución analítica en el caso n = 1 y con el modelo sin grav. 

    Returns:
     Gráfica con escalas logarítmicas con el error de los distintos métodos en la posición y la velocidad

    '''

    t_simb = np.linspace(0, T, 10**2)
    ts = sp.symbols("ts", real=True)
    ys = sp.symbols("ys", cls=sp.Function)
    cis = {ys(0): 0, ys(ts).diff(ts).subs(ts, 0): 0}

    ode = sp.Eq(ys(ts).diff(ts, 2), u*mu_0/(m_t(ts, 1, mu_0, T, mc, k)*T))
    sol = sp.dsolve(ode, ys(ts), ics=cis)
    Yexac = sol.rhs
    DYexac = Yexac.diff(ts)
    # Vector de posiciones solución analítica
    YY = np.array([Yexac.subs(ts, z) for z in t_simb])
    # Vector de velocidades solución analítica
    DY = np.array([DYexac.subs(ts, z) for z in t_simb])

    t_in = np.array([0, T])

    plt.figure(1, figsize=(12, 12), clear=True)
    plt.suptitle(suptitle)

    # plt.rcParams['text.usetex'] = True

    '''Altura respecto al tiempo'''
    ax1 = plt.subplot(2, 1, 1)
    # plt.title(r'$\textbf{\frac{1}{m_0}}$', fontsize=15)
    plt.title(r'Errores en la altura', fontsize=15)
    plt.xlabel("t(s)")
    plt.grid()

    ''' Velocidad respecto al tiempo'''
    ax2 = plt.subplot(2, 1, 2)
    # plt.title(r'$\textbf{$\frac{1}{m_0}\frac{d}{dt}m(t)}$', fontsize=15)
    plt.title(r'Errores en la velocidad', fontsize=15)
    plt.xlabel("t(s)")
    plt.grid()

    def odefun(t, Y): return odefunGen(t, Y, modelo='singrav')

    for met in methods:

        result = sin.solve_ivp(odefun, t_in, cin, t_eval=t_simb, method=met)
        tsol = result.t
        v = result.y[1]
        y = result.y[0]
        ax1.semilogy(tsol, abs(YY-y), label=str(met))
        ax2.semilogy(tsol, abs(DY-v))

    ax1.legend()
    plt.show()


'''
-----------------------------------------------------------------------------------------------------------------------------------------------
FUNCIONES PARA LA REPRESENTACIÓN DE VARIAS ETAPAS
-----------------------------------------------------------------------------------------------------------------------------------------------  
'''

def RepresEtapas(ax1,ax2,modelo='singrav',n_list=n_list):
    
    '''
    Parameters:
     - ax1,ax2: subplots creadas anteriormente con CreacionPlot
     - modelo (string): modelo que se quiere estudiar
     - n_list (lista): n que se quieren representar

    Objective:
     Representar la posicion y velocidad del cohete de tres etapas para cada n

    Returns:
     Plot sin representar
    '''

    for n in n_list:

        def odefun(t, Y): return odefunGen_etap(t, Y, modelo=modelo, n=n)

        result = sin.solve_ivp(odefun, t_in_etap, cin, t_eval=tspan_etap)
        # result = sin.solve_ivp(odefun, t, [0], t_eval=tspan)
        tsol = result.t
        y = result.y[0]
        v = result.y[1]

        ax1.plot(tsol, y, label=f'n={n}')
        ax2.plot(tsol, v, label=f'n={n}')

        if n==n_list[-1]:
            ax2.axvline(x=T1,color='k',linestyle='-.',label='Fin Etapa 1')
            ax2.axvline(x=T2+T1, color='k', linestyle='--',label='Fin Etapa 2')
            ax2.axvline(x=T3+T2+T1, color='k', linestyle='-', label='Fin Etapa 3')
            ax1.axvline(x=T1, color='k',linestyle='-.', label='Fin Etapa 1')
            ax1.axvline(x=T1+T2, color='k', linestyle='--',label='Fin Etapa 2')
            ax1.axvline(x=T1+T2+T3, color='k',linestyle='-', label='Fin Etapa 3')

        ax1.legend()


def ComparacionEtap(ax1,ax2,modelo='singrav',n=1):

    '''
    Parameters:
     - ax1,ax2: subplots creadas anteriormente con CreacionPlot
     - modelo (string): modelo que se quiere estudiar
     - n (int): cte de consumo usada

    Objective:
     Comparar la trayectoria del cohete de 1 etapa con el de 3

    Returns:
     Plot sin representar
    '''

    def odefun(t, Y): return odefunGen_etap(t, Y, modelo=modelo, n=n)

    result = sin.solve_ivp(odefun, t_in_etap, cin, t_eval=tspan_etap)
    # result = sin.solve_ivp(odefun, t, [0], t_eval=tspan)
    tsol = result.t
    y = result.y[0]
    v = result.y[1]

    ax1.plot(tsol, y, label=f'Varias Etapas')
    ax2.plot(tsol, v)

    def odefun(t, Y): return odefunGen(t, Y, modelo=modelo, n=n, T = T1+T2+T3, mu_0=sum(mulist),mc = sum(msecolist))

    result = sin.solve_ivp(odefun, t_in_etap, cin, t_eval=tspan_etap)
    tsol = result.t
    v = result.y[1]
    y = result.y[0]

    ax1.plot(tsol, y, label=f'Una etapa')
    ax2.plot(tsol, v)

    ax2.axvline(x=T1, color='k', linestyle='-.', label='Fin Etapa 1')
    ax2.axvline(x=T2+T1, color='k', linestyle='--', label='Fin Etapa 2')
    ax2.axvline(x=T3+T2+T1, color='k',linestyle='-', label='Fin Etapa 3')
    ax1.axvline(x=T1, color='k', linestyle='-.', label='Fin Etapa 1')
    ax1.axvline(x=T1+T2, color='k', linestyle='--',label='Fin Etapa 2')
    ax1.axvline(x=T1+T2+T3, color='k',linestyle='-', label='Fin Etapa 3')

    ax1.legend()


