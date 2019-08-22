# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from numpy import *
from scipy.integrate import odeint
from scipy import signal
from copy import deepcopy
#from SALib.sample import saltelli
#from SALib.analyze import sobol
import time
import pickle
import xlrd
from xlwt import Workbook
from scipy.stats import linregress
from mpl_toolkits.mplot3d import Axes3D


def eqs(Y, t, TabConst, bug):
    
    """Résout les equa-diff : Y contient les variables à calculer (concentrations,
    recepteurs occupés...), t représente le temps, TabConst contient les valeurs des constantes
    présentes dans les équations, bug est juste là pour éviter un bug"""
    
    #On récupère les valeurs
    c1, c21, c22, c3, c41, c42, c5, c6, ros, rom = Y
    clp, V, alpha, clm , fpe , ka1, ka21 , ka22 , lag1 , lag2 , bio, MMI , frac , add_p , prop_p , add_m , prop_m, b0 , a1 , phi1 , a2 , phi2 ,kons , koffs , Emax , sigmaadd ,aon , boff = TabConst
    ka3 = clp / V
    ka4 = clm / V
    konm = aon*kons
    koffm = boff*koffs
    
    #On dérive 
    dc1 = -ka1*c1
    dc21 = -ka21*c21
    dc22 = -ka22*c22
    dc3 = -ka1*c3
    dc41 = -ka21*c41
    dc42 = -ka22*c42
    dc5 = ka1*c1 + ka21*c21 + ka22*c22 - ka3*c5
    dc6 = ka1*c3 + ka21*c41 + ka22*c42 - ka4*c6 + alpha*ka3*c5
    dros = kons*c5*(1 - ros - rom) - koffs*ros
    drom = konm*c6*(1 - ros - rom) - koffm*rom
    
    #On renvoie la valeur des dérivées
    return [dc1, dc21, dc22, dc3, dc41, dc42, dc5, dc6, dros, drom]


def base1(t, TabConst):
    
    """Renvoie le rhythme cardiaque de base du patient, selon le groupe SBT07 """
    b0 , a1 , phi1 , a2 , phi2 =TabConst[17:22]
    Emax=TabConst[24]
    return b0*(1 + a1*cos(2*pi*t*(1/24) + phi1) + a2*cos(4*pi*t*(1/6) + phi2))


def base2(t, TabConst):
    
    """La base a priori telle que donnée dans le rapport"""
    b0 , a1 , phi1 , a2 , phi2 =TabConst[17:22]
    Emax=TabConst[24]
    return b0*(1 + a1*cos(2*pi*t/24 + phi1) + a2*cos(4*pi*t/24 + phi2))
    
#Qui choisir entre base1 et base2 ??? A priori ça n'a pas d'influence pour l'étude de la 
#variabilité/sensibilité, ça change surtout l'allure...

def effet(Y, T, TabConst):
    
    """Renvoie le tableau contenant les battements du patient avec traitement"""
    
    Emax = TabConst[24]
    Ros = extraire_var(Y, "ros")
    Rom = extraire_var(Y, "rom")
    B = [ base1(T[i], TabConst)*(1 - Emax*(Ros[i] + Rom[i])) for i in range(len(Y)) ]
    return B 


def calcul(Yinit, t1, t2, TabConst):
    
    """Calcule l'evolution des variables entre t1 et t2 à partir de la config Yinit"""
    
    T = linspace(t1, t2, floor((t2-t1)*ptsperh))
    return T, odeint(eqs, Yinit, T, args = (TabConst, 0))
    

def lag(Inj, TabConst):
    
    """Renvoie un tableau des injections prenant en compte le lag et la proportion de
    chaque dose allant dans chaque compartiment"""
    lag1 , lag2 , bio, MMI , frac = TabConst[8:13]
    fpe=TabConst[4]
    Inj2 = []
    
    for (tinj, cinj) in Inj:
        #On répartit la dose dans les 4 compartiments à des temps différents
        comp1 = [tinj + lag1, bio*frac*cinj, 0]
        comp21 = [tinj + lag2, bio*(1-frac)*cinj, 1]
        comp3 = [tinj + lag1, bio*fpe*frac*cinj, 3]
        comp41 = [tinj + lag2, bio*fpe*(1-frac)*cinj, 4]
        #On prévoit le changement de compartiment
        comp22 = [tinj + lag2 + MMI, -1, 2]
        comp42 = [tinj + lag2 + MMI, -1, 5]
        Inj2.extend([comp1, comp21, comp3, comp41, comp22, comp42])
        
    Inj2.sort()
    return Inj2


def simul(yinit, Inj, TabConst):
    
    """Renvoie l'evolution de toutes les variables (dont l'effet) durant la simulation"""
    
    Yinit = deepcopy(yinit)
    T = [0]
    Res = [Yinit]
    t = 0
    Inj2 = lag(Inj, TabConst)#Prise en compte du lag/proportions de la dose dans les compartiments
    
    for i in range(len(Inj2)):
    
        tinj, cinj, cible = Inj2[i]
        
        if t < tinj : #On fait la simulation jusqu'à la prochaine injection
            T2, Res2 = calcul(Yinit, t, tinj, TabConst)
            Res = concatenate((Res, Res2))
            T = concatenate((T, T2))
        
        if cinj != -1:
            Res[-1][cible] += cinj#Le compartiment-cible reçoit la dose
        else:
            Res[-1][cible] += Res[-1][cible - 1]#on fait passer c21 dans c22
            Res[-1][cible - 1] = 0
        
        t = tinj#On rétablit les bonnes valeurs
        Yinit = deepcopy(Res[-1])
        
    if t < tmax:#on finit la simulation
        T2, Res2 = calcul(Yinit, t, tmax, TabConst)
        Res = concatenate((Res, Res2))
        T = concatenate((T, T2))
        
    return T, Res, effet(Res, T, TabConst)


def extraire_var(Res, nom_var):
    
    """Mettre en argument "c1", "c2", "c3", "c4", "c5", "c6", "ros", "rom" pour obtenir
    le tableau de valeurs correspondant à cette variable"""
    
    n = len(Res)
    var = ["c1", "c2", "c3", "c4", "c5", "c6", "ros", "rom"]
    indices = [[0], [1, 2], [3], [4, 5], [6], [7], [8], [9]]
    ind = indices[var.index(nom_var)]
    if len(ind) == 1:
        return [Res[i][ind[0]] for i in range(n)]
    else:
        return [Res[i][ind[1]] + Res[i][ind[1]] for i in range(n)]


#Toutes les constantes clp, v, alpha, clm , fpe , ka1, ka21 , ka22 , lag1 , lag2 , bio, MMI , frac , add_p , prop_p , add_m , prop_m,b0 , a1, phi1 , a2 , phi2 ,kons , koffs , Emax , sigmaadd ,aon , boff
TabConst=[160,201,0.031,55.6,0.22,0.158,0.0073,0.00271,0.754,4.97,1,5.59,0.426,0,24.7,0.3,16.2,64.1,-0.0332,0.987,0.0204,0.588,0.0175,0.085,0.364,3.97,6,10]

# Tableau des variances IIV et IOV 
IIV=[0,.471,0,.246,0,0,.638,0,0,0,.438,0,0,0,0,0,0,.0718,0,0,0,0,.885,.581,0,0,0,0]
IOV=[0,.322,0,.183,0.333,0,.66,1.41,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]


#Les valeurs de toutes les variables à l'état initial (si injection à t=0, mettre dans Inj)
Yinit = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

#Temps de simulation (en h)
global tmax
tmax = 24*20

#Nombre de pts calculés par h
global ptsperh
ptsperh = 10

#Calendrier des injections : (instant de l'injection en h, dose injectée en mg (?))
Inj = [(24*i, 50) for i in range(0, tmax//24 + 1)]
Inj = [(0, 50)]


#Script typique pour faire et tracer une simulation
#T le tableau des temps
#Res contient l'evolution de [c1, c21, c22, c3, c41, c42, c5, c6, ros, rom]
#Pour choper l'une des variables, utiliser extraire_var
#E contient l'évolution des battements du patient avec traitement
'''plt.close()
T, Res, E = simul(Yinit, Inj, tmax, TabConst)
plt.plot(T, E, linewidth = 0.5)
plt.show()'''

def varConst():
    
    TabConst2 = deepcopy(TabConst)
    for j in range(len(TabConst)):
        if IIV[j] != 0 :
            TabConst2[j]=TabConst2[j]*exp(random.normal(0,IIV[j]))
        if IOV[j] != 0 :
            TabConst2[j]=TabConst2[j]*exp(random.normal(0,IOV[j])) 
    return TabConst2

#n est le nombre de fois où on effectue les tests, pour réaliser une étude statistique
def simulation_globale(n):
    resultat_T=[]
    resultat_Res=[]
    resultat_E=[]
    Variables=[]
    for i in range(n):  
        
        print("Pourcentage d'avancement : ", str(i/n*100), " %")
    
          
        #prend le tableau des variables et perturbe les variables suivant une loi normale de variance (iiv+iov)
        TabConst2 = varConst()          
                    
        """ on effectue une simulation pour ces valeurs obtenues"""
        T, Res, E = simul(Yinit, Inj, TabConst2)
        resultat_T.append(T)
        resultat_Res.append(Res)
        resultat_E.append(E)
        Variables.append(TabConst2)
    return resultat_T,resultat_Res,resultat_E,Variables


"""enregistrement donnees"""
#On exécute cette fonction une fois, avec n fixé. Ensuite, on peut accéder à tous les résultats sans la reéxecuter:
#on peut travailler sur le tableau donnees comme ci-dessous.
def enregistrer_donnees(n):
    tableau_T,tableau_Res,tableau_E,Variables=simulation_globale(n)
    tableau_total=[tableau_T,tableau_Res,tableau_E,Variables]
    Fichier=open('tableau_total.txt','wb')
    pickle.dump(tableau_total,Fichier)
    Fichier.close()
    

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
SEUL ENDROIT QUI EXECUTE UNE FONCTION, REPRESENTATIONS GRAPHIQUES EXCLUES. ON L'EXECUTE UNE FOIS PUIS ON LA PASSE EN COMMENTAIRE POUR TRAVAILLER SUR LES DONNEES OBTENUES ET NE PAS LES RECALCULER A CHAQUE FOIS"""

# enregistrer_donnees(3)

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


#donnees de la forme: [tableau_T,tableau_Res,tableau_E]
#avec tableau_T=[T pour la simulation 1, ........., T pour la simulation n]
#avec tableau_Res=[Res pour la simulation 1, ........., Res pour la simulation n]
#avec tableau_E=[E pour la simulation 1, ........., E pour la simulation n]


def recuperer_donnees():
    
    """Recupère les données à partir du fichier"""
    
    fichier=open('tableau_total.txt','rb')
    donnees=pickle.load(fichier)
    fichier.close()
    return donnees

def extraire_tab(donnees):
    
    """Extrait les trois tableaux principaux à partir des données"""
    
    n=len(donnees[0])
    T=donnees[0][0]
    tableau_Res=donnees[1]
    tableau_E=donnees[2]
    Variables=donnees[3]
    return T, tableau_Res, tableau_E,Variables

def extraire_tab_var(donnees, nom_var):
    
    """Permet d'extraire des données le tableau d'évolution d'une des variables parmi
    "c1", "c2", "c3", "c4", "c5", "c6", "ros", "rom" (à mettre en argument 
    nom_var)"""
    
    return [extraire_var(tabRes[i], nom_var) for i in range(len(tabRes))]


def tracer_courbe_stat(T, tab):
    
    """Trace moyenne, mediane, quantiles à 5/95 % pour une variable donnée dont les valeurs
    sont contenues dans tab. T reste le temps."""
    
    #n nombre de simualtions réalisées
    #calcul de l'évolution des differents parametres de l'effet
    moyenne=mean(tab,axis=0)
    mediane=median(tab,axis=0)
    variance=var(tab,axis=0)
    q5Percent=percentile(tab,5,axis=0)
    q95Percent=percentile(tab,95,axis=0)
    #on trace
    plt.plot(T, moyenne, 'r-', linewidth = 0.5, label="moyenne")
    plt.plot(T, mediane, 'b-', linewidth = 0.5, label="mediane")
    plt.plot(T, q5Percent, 'g-', linewidth = 0.5, label="quantiles a 5% et 95%")
    plt.plot(T, q95Percent, 'g-', linewidth = 0.5)
    plt.legend()
    plt.show()
    
    
def analyse_sensibilite_sobol(metrique,samplesize):
    # indices des variables dans le tableau Tabconst 
    indices=[1,3,4,6,7,10,17,22,23]
    problem = {
        'num_vars':9,
        'names':['v','clm','fpe','ka21','ka22','bio','b0','kons','koffs'],
        'bounds':[[TabConst[i]*exp(-3*(IIV[i]+IOV[i])),TabConst[i]*exp(3*(IIV[i]+IOV[i]))] for i in indices]#regle des 3 sigma
    }
    # Generate samples
    param_values = saltelli.sample(problem, samplesize)
    
    # Run model
#     Y= tableau de la moyenne des battements (metrique qu'on peut changer) pour les entrée param_values 
    Y=zeros(len(param_values))
    for i in range (len(param_values)):
        TabSimul=deepcopy(TabConst)
        for j in range (len(indices)):
            TabSimul[indices[j]]=param_values[i][j]
        Y[i]=metrique(Yinit, TabSimul)

    # Perform analysis
    Si = sobol.analyze(problem, Y)
    return(Si)
    
    
def plot_sensibilite_sobol(Si) :
    # plot https://matplotlib.org/examples/lines_bars_and_markers/barh_demo.html
    plt.close()
    plt.rcdefaults()
    fig, ax = plt.subplots()
    
    Variables = ('v','clm','fpe','ka21','ka22','bio','b0','kons','koffs')
    y_pos = arange(len(Variables))
    sobol_index = Si['ST'] 
    error = Si['ST_conf']
    
    
    ax.bar( y_pos , sobol_index,yerr=error, align='center', ecolor='blue' , color='orange')
    ax.set_xticks(y_pos)
    ax.set_xticklabels(Variables)
    ax.set_xlabel('Coefficient')
    ax.set_title('Indice Sensibilité totale par méthode de Sobol')
    
    plt.show()
'''
Pour avoir un plot pour l'analyse de sensibilité executer le script :

Si=analyse_sensibilite_sobol(metrique,100) # vous pouvez changer la taille de l'échantillon 
plot_sensibilite_sobol(Si)

'''    
# enregistrer_donnees(1000)
''''T, tabRes, tabE ,Variables= extraire_tab(recuperer_donnees())
tabC5 = extraire_tab_var(tabRes, "c5")
tracer_courbe_stat(T, log(tabC5))'''


def moyennebat(Yinit, TabConst, dose):#calcule la moyenne des battements entre le jour 15 et le dernier jour d'intégration.tabE est à une dimension, tmax la durée d'intégration, ptsperh le nombre d'intégrations par h
    
    if tmax<=24*19: #Vérification que la durée d'itégration est d'au moins 19 jours pour avoir une moyenne pertinente.
        print("Incorrect integration length. The simulation must be at least 19 days long.")
    else:
        Inj = [(24*i, dose) for i in range(0, tmax//24 + 1)]
        T, Res, tabE = simul(Yinit, Inj, TabConst)
        ind_debut=24*15*ptsperh
        delta_t=int(len(tabE)-ind_debut)
        ind_fin=int(ind_debut+24*ptsperh*floor(delta_t/(24*ptsperh)))#on cherche à moyenner sur un nombre entier de jours.
        return sum(tabE[ind_debut:ind_fin])/(ind_fin-ind_debut)
    

def dose_dichoto(Yinit, TabConst, fbat, eps):
    
    """Renvoie la dose necessaire à eps près pour obtenir des battements à la fréquence fbat 
    pour un patient décrit par les constantes TabConst à partir d'une dichotomie"""
    
    dose1 = 0
    bat1 = moyennebat(Yinit, TabConst, dose1)
    dose2 = 200
    bat2 = moyennebat(Yinit, TabConst, dose2)
    while dose2 - dose1 > eps:
        dosem = 0.5*(dose1 + dose2)
        batm = moyennebat(Yinit, TabConst, dosem)
        if batm > fbat:
            dose1 = dosem
            bat1 = batm
        else:
            dose2 = dosem
            bat2 = batm
    return 0.5*(dose1 + dose2)


def moyenne_stats(n, moyenne_bat, tol):
    
    dose = dose_dichoto(Yinit, TabConst, moyenne_bat, 10**-3)
    M = []
    T = []
    for i in range(n):
        print(i/n*100)
        V = varConst()
        M.append(moyennebat(Yinit, V, dose))
        T.append(V)
    return M, T


def coeff_reg(Yinit, TabConst, nbpoints = 3):
    
    """Renvoie les coefficients de la regression linéaire pour les constantes données"""
    
    Dose = linspace(0, 200, nbpoints)
    E = array([moyennebat(Yinit, TabConst, dose) for dose in Dose])
    E0 = E[0]
    Einf = moyennebat(Yinit, TabConst, 10**9)
    Ecor = (E0 - E)/(E-Einf)
    slope, intercept, rvalue, pvalue, stderr = linregress(Dose, Ecor)
    return slope


def creer_excel(nom, nom_donnees, donnees):
    
    """Cree un Excel du nom précisé contenant les données précisées"""
    
    # création 
    book = Workbook()
    # création de la feuille 1
    feuil1 = book.add_sheet('feuille 1')
    # Remplissage
    for i in range(len(nom_donnees)):
        feuil1.write(0, i, nom_donnees[i])
    for i in range(len(donnees)):
        for j in range(len(donnees[i])):
            feuil1.write(i+1, j, donnees[i][j])
    #Sauvegarde
    book.save(nom + ".xls")
    
    
def creer_donnees_simul(n, nomvars, ivars, nomfichier = "test"):
    
    """Ecrit sur un Excel le resultat de n simulations aleatoires avec les variables
    dont sont spécifiés le nom et l'indice dans les tableaux nomvars et ivars"""
    
    #Indices pour les variables intéressantes
    #v, clm, ka21, ka22, b0, kons, koffs --> 1, 3, 6, 7, 17, 22, 23
    
    nvar = len(ivars)
    donnees = zeros((n, nvar + 1))
    for i in range(n):
        print(i/n*100)
        TabConst2 = deepcopy(TabConst)
        for j in ivars:
            if IIV[j] != 0 :
                TabConst2[j]=TabConst2[j]*exp(random.normal(0,IIV[j]))
            if IOV[j] != 0 :
                TabConst2[j]=TabConst2[j]*exp(random.normal(0,IOV[j]))
        TabConst2 = varConst()
        slope = coeff_reg(Yinit, TabConst2)
        for j in range(nvar):
            donnees[i][j] = TabConst2[ivars[j]]
        donnees[i][nvar] = slope
        
    nomvars.append("pente")
    creer_excel(nomfichier, nomvars, donnees)
    return donnees


def evol_pente(iv, nbpoints, f = lambda x : x, invf = lambda x : x):
    
    """Renvoie l'évolution de la pente en fonction de la variable dont l'indice est donné
    (voir fonction du dessus pour quel indice correspond à quelle variable)"""
    
    a, b = TabConst[iv]*exp(-3*(IIV[iv]+IOV[iv])),TabConst[iv]*exp(3*(IIV[iv]+IOV[iv]))
    fa, fb = min(f(a), f(b)), max(f(a), f(b))
    valvar = linspace(fa, fb, nbpoints)
    TabConst2 = deepcopy(TabConst)
    pente = []
    for i in range(len(valvar)):
        print(i/len(valvar)*100)
        TabConst2[iv] = invf(valvar[i])
        slope = coeff_reg(Yinit, TabConst2)
        pente.append(slope)
    plt.plot(valvar, pente)
    plt.show()
    return array(valvar), array(pente)


def sobol2(samplesize):
    # indices des variables dans le tableau Tabconst 
    indices=[1,10,22,23]
    problem = {
        'num_vars':4,
        'names':['v','bio','kons','koffs'],
        'bounds':[[TabConst[i]*exp(-3*(IIV[i]+IOV[i])),TabConst[i]*exp(3*(IIV[i]+IOV[i]))] for i in indices]#regle des 3 sigma
    }
    # Generate samples
    param_values = saltelli.sample(problem, samplesize)
    
    # Run model
#     Y= tableau de la moyenne des battements (metrique qu'on peut changer) pour les entrée param_values 
    Y=zeros(len(param_values))
    print(len(param_values))
    for i in range (len(param_values)):
        print(i/len(param_values)*100)
        TabSimul=deepcopy(TabConst)
        for j in range (len(indices)):
            TabSimul[indices[j]]=param_values[i][j]
        Y[i] = coeff_reg(Yinit, TabSimul)

    # Perform analysis
    Si = sobol.analyze(problem, Y)
    return(Si)
    
    
def plot_sobol2(Si) :
    # plot https://matplotlib.org/examples/lines_bars_and_markers/barh_demo.html
    plt.close()
    plt.rcdefaults()
    fig, ax = plt.subplots()
    
    Variables = ('v','bio','kons','koffs')
    y_pos = arange(len(Variables))
    sobol_index = Si['ST'] 
    error = Si['ST_conf']
    
    
    ax.bar( y_pos , sobol_index,yerr=error, align='center', ecolor='blue' , color='orange')
    ax.set_xticks(y_pos)
    ax.set_xticklabels(Variables)
    ax.set_xlabel('Coefficient')
    ax.set_title('Indice Sensibilité totale par méthode de Sobol')
    
    plt.show()
    
    
def evol_pente3d(Iv, nbpts):
    
    ivx, ivy = Iv
    
    a, b = TabConst[ivx]*exp(-3*(IIV[ivx]+IOV[ivx])),TabConst[ivx]*exp(3*(IIV[ivx]+IOV[ivx]))
    X = linspace(a, b, nbpts)
    a, b = TabConst[ivy]*exp(-3*(IIV[ivy]+IOV[ivy])),TabConst[ivy]*exp(3*(IIV[ivy]+IOV[ivy]))
    Y = linspace(a, b, nbpts)
    
    def f(x,y):
        TabConst2 = deepcopy(TabConst)
        TabConst2[ivx] = x
        TabConst2[ivy] = y
        return coeff_reg(Yinit, TabConst2)
    
    return plot3d(X, Y, f)
    

def plot3d(X, Y, f):
    
    X3D = array([X for i in range(len(Y))])
    Y3D = array([Y for i in range(len(X))])
    Y3D = transpose(Y3D)
    Z = array([[ f(X[i], Y[j]) for i in range(len(X)) ] for j in range(len(Y))])
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_wireframe(X3D, Y3D, Z)
    return X3D, Y3D, Z


def plothisto(L, a, b, n = -1, h = -1):
    
    if n == -1:
        n = int((b-a)/h)
    h = (b-a)/n
    H = zeros(n)
    for i in range(len(L)):
        index = int(floor((L[i]-a)/h))
        if index >= 0 and index < n:
            H[index] += 1
    X = [a + (0.5+i)*h for i in range(n)]
    plt.plot(X, H)
    return X, H

for i in range(100):
    V = varConst()
    print(V[17])


#creer_donnees_simul(1000, ['v','clm','fpe','ka21','ka22','bio','b0','kons','koffs'], 
#                    [1,3,4,6,7,10,17,22,23])
#Notes
#'v','clm','fpe','ka21','ka22','bio','b0','kons','koffs' -> 1,3,4,6,7,10,17,22,23
#pente linéaire en v : 6e-5
#pente linéaire en 1/koffs : 0.001
#pente linéaire en kons : 0.6
#pente linéaire en bio : 0.01
#pente linéaire en 1/clm
#pente linéaire en fpe
#pente peu influée par ka21
#pour ka22, un peu étrange, d'abord log puis constant...
#b0 sans AUCUNE influence