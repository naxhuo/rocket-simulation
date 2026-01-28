

import matplotlib.pyplot as plt
from GlobalFunctionsFin import *


'''
-----------------------------------------------------------------------------------------------------------------------------------------------
ESTUDIO DEL COMPORTAMIENTO DE LA MASA CON DISTINTAS CTES DE CONSUMO
-----------------------------------------------------------------------------------------------------------------------------------------------
'''

RepresMasas()


'''
-----------------------------------------------------------------------------------------------------------------------------------------------
DISTINTOS MODELOS POR INTEGRACIÓN
-----------------------------------------------------------------------------------------------------------------------------------------------
'''

# Representación en cada modelo por separado

for model in models[:-1]:
	ax1,ax2,ax3 = CreacionPlot(suptitle='Resuelto por Integración ' + model)
	IntegRepresN(pos_t,ax1,ax2,ax3, model = model)
	plt.show()

# ax1, ax2, ax3 = CreacionPlot(suptitle='Integración y gravcte')
# IntegRepresN(pos_t, ax1, ax2, ax3, model=models[1])
# plt.show()

# Representación de todos los modelos juntos

ax1, ax2, ax3 = CreacionPlot(suptitle='Comparación de los modelos usando integración')
IntegRepresModel(pos_t, ax1,ax2,ax3,n=1,models = models[:-1])
plt.show()


'''
-----------------------------------------------------------------------------------------------------------------------------------------------
REPRESENTACIÓN ODE'S
-----------------------------------------------------------------------------------------------------------------------------------------------
'''

#Representación de cada modelo por separado
for model in models:
	ax1, ax2, ax3 = CreacionPlot(suptitle="ODE's "+model)
	OdeSystemRepres(ax1,ax2,ax3,modelo=model,n_list=[0.25,0.5,0.75,1])
	plt.show()

# ax1, ax2, ax3 = CreacionPlot(suptitle='rozcte')
# OdeSystemRepres(ax1, ax2, ax3, modelo=models[1])
# plt.show()

ax1, ax2, ax3 = CreacionPlot(suptitle="ODE's rozcte")
OdeSystemRepres(ax1, ax2, ax3, modelo='rozcte',n_list=n_list[:-2])
plt.show()

# Representación de todos los modelos juntos
ax1,ax2 = CreacionPlot(False,suptitle='Comparación de los modelos usando EDOs')
OdeSystemRepresModels(ax1,ax2,models=models)
plt.show()

# Representación de todos los métodos juntos
ax1, ax2 = CreacionPlot(False,suptitle='Comparación distintos métodos con rozcte')
OdeSystemRepresMethods(ax1,ax2,model = 'rozcte', methods=methods)
plt.show()


'''
------------------------------------------------------------------------------------------------------------------------------------------------
COMPARACIÓN MÉTODOS DE RESOLUCIÓN DE EDO NUMÉRICOS CON SOLUCIÓN ANALÍTICA PARA EL MODELO SIN GRAVEDAD
------------------------------------------------------------------------------------------------------------------------------------------------
'''

ComparacMetodos(methods = methods)


'''
-----------------------------------------------------------------------------------------------------------------------------------------------
COHETE POR ETAPAS (3 ETAPAS)
-----------------------------------------------------------------------------------------------------------------------------------------------
'''

# Representación del modelo de 3 etapas para distintas n con el modelo elegido

ax1,ax2 = CreacionPlot(False,suptitle='Cohete tres etapas (Modelo: Rozamiento variable)')
#RepresEtapas(ax1, ax2, modelo='singrav',n_list=n_list)
RepresEtapas(ax1, ax2, modelo='rozcte',n_list=n_list[:-2])#,n_list=n_list)
plt.show()


# Comparación de 1 etapa con 3 para n=1
ax1, ax2 = CreacionPlot(False, suptitle='Comparación una etapa vs tres etapas (Modelo: Sin Gravedad)')
ComparacionEtap(ax1,ax2,modelo='singrav')
plt.show()



'''
-----------------------------------------------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------------------------------------------
CREACIÓN DE UN COHETE CAPAZ DE LLEVAR UN SATÉLITE DE 1000 KG A ÓRBITA BAJA
-----------------------------------------------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------------------------------------------
'''

print('''
-------------------------------------------------------------------------------
Ahora diseña tu un cohete que llegue a órbita baja con un satélite de 1000 kg!!
Nota: Todo cáracter no correspondiente a un número entero será tomado como un dato ya preestablecido, de tal manera que todos los preestablecidos funcionan
-------------------------------------------------------------------------------
''')



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


def check_int(intro, other):
    if intro.isdigit():
        return int(intro)
    else:
        return other


# mseco_et1: Masa en seco de la etapa 1
mseco_et1 = check_int(input('Masa en seco (kg) de la etapa 1: '), 13700)
# mcom_et1: Masa de combustible de la etapa 1
mcom_et1 = check_int(input('Masa combustible (kg) de la etapa 1: '), 200000)

# mseco_et2: Masa en seco de la etapa 2
mseco_et2 = check_int(input('Masa en seco (kg) de la etapa 2: '), 3600)
# mcom_et2: Masa de combustible de la etapa 2
mcom_et2 = check_int(input('Masa combustible (kg) de la etapa 2: '), 6000)

mseco_et3base = check_int(input('Masa en seco (kg) de la etapa 3: '), 3000)
msatelite = 1000
mseco_et3 = mseco_et3base + msatelite  # mseco_et3: Masa en seco de la etapa 3
# mcom_et3: Masa de combustible de la etapa 3
mcom_et3 = check_int(input('Masa combustible (kg) de la etapa 3: '), 1090)

msecolist = [mseco_et1, mseco_et2, mseco_et3]
mulist = [mcom_et1, mcom_et2, mcom_et3]


"""
--------------------------------------------
Tiempo de quema de combustible de las etapas
--------------------------------------------
"""
T1 = check_int(input('Tiempo (s) consumo de la etapa 1: '),
               150)  # Tiempo de quema de combustible de la primera etapa
# Tiempo de quema de combustible de la segunda etapa
T2 = check_int(input('Tiempo (s) consumo de la etapa 2: '), 300)
# Tiempo de quema de combustible de la tercera etapa
T3 = check_int(input('Tiempo (s) consumo de la etapa 3: '), 400)

Tlist_etap = [T1, T2, T3]

tspan_etap = np.arange(0, T2+T1+T3+500)
t_in_etap = [0, T1+T2+T3+500]


"""
-------------------------------------------------
Velocidad de expulsión de los gases de las etapas
-------------------------------------------------
"""
u1 = check_int(input('Velocidad salida gases (m/s) etapa 1: '),
               5000)  # Velocidad de expulsión de gases de la primera etapa
# Velocidad de expulsión de gases de la segunda etapa
u2 = check_int(input('Velocidad salida gases (m/s) etapa 2: '), 3000)
# Velocidad de expulsión de gases de la tercera etapa
u3 = check_int(input('Velocidad salida gases (m/s) etapa 3: '), 3000)


ulist_etap = [u1, u2, u3]

print(f'''
----------------------------------
Datos introducidos:

Masas en seco: {msecolist}
Masas combustible: {mulist}
Períodos: {Tlist_etap}
Salida Gases: {ulist_etap}
----------------------------------
''')



def odefunGen_etap_mivicho(t, Y, modelo='singrav', n=0.25, mu_0=mulist, mc=msecolist, T=Tlist_etap, k=k, lat=latcape_canaveral, A=A, cd=cd, gamma=gamma, u_etap=u_etap):
    '''
    Parameters:
     - t (int): tiempos
     - Y (list): [posición,velocidad]
     - Modelo de fricción: A: sección, cd: cte dependiente de la forma, gammma: cte usada en el modelo exponencial de densidad del aire
     - + Variables antes descritas: nu_0, mc, T en formato lista

    Objective:
        Calcular el SEDO que nos interese, para tres etapas con los datos personalizados del cohete

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
        return [Y[1], ((-p0*np.exp(-gamma*Y[0]))*0.5*A*cd*Y[1]**2 - u_etap(t, Ts=T, us=ulist_etap)*dm_t_etap(t, n, mu_0, T, k) - m_t_etap(t, n, mu_0, T, mc, k)*g_y2(Y[0], lat=latcape_canaveral))/m_t_etap(t, 1, mu_0, T, mc, k)]


def RepresEtapas_mivicho(ax1, ax2, modelo='singrav', n=0.25):
    '''
    Parameters:
     - ax1,ax2: subplots creadas anteriormente con CreacionPlot
     - modelo (string): modelo que se quiere estudiar


    Objective:
     Representar la posicion y velocidad del cohete de tres etapas con los datos personalizados del cohete

    Returns:
     Plot sin representar
    '''

    def odefun(t, Y): return odefunGen_etap_mivicho(t, Y, modelo=modelo, n=n)

    result = sin.solve_ivp(odefun, t_in_etap, cin, t_eval=tspan_etap)
    # result = sin.solve_ivp(odefun, t, [0], t_eval=tspan)
    tsol = result.t
    y = result.y[0]
    v = result.y[1]

    if np.any(y < 0):
        print('Llega a alturas negativas!!')
    else:
        if (y[T1+T2+T3+200] > 150000) & (y[T1+T2+T3+200] < 2000000):
            print('Has conseguido que llegue a órbita baja!!')
    ax1.plot(tsol, y, label=f'n={n}')
    ax2.plot(tsol, v, label=f'n={n}')

    ax2.axvline(x=T1, color='k', linestyle='-.', label='Fin Etapa 1')
    ax2.axvline(x=T2+T1, color='k', linestyle='--', label='Fin Etapa 2')
    ax2.axvline(x=T3+T2+T1, color='k', linestyle='-', label='Fin Etapa 3')
    ax1.axvline(x=T1, color='k', linestyle='-.', label='Fin Etapa 1')
    ax1.axvline(x=T1+T2, color='k', linestyle='--',  label='Fin Etapa 2')
    ax1.axvline(x=T1+T2+T3, color='k',  linestyle='-', label='Fin Etapa 3')

    ax1.axhline(y=150000, color='g',  linestyle='-') 
    ax1.axhline(y=2000000, color='g',  linestyle='-')
    ax2.axhline(y=0, color='k')
    ax1.legend()


ax1, ax2 = CreacionPlot(False, suptitle='Varias etapas')
RepresEtapas_mivicho(ax1, ax2, modelo='rozcte')
plt.show()
