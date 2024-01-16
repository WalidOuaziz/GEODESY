from math import *
from fonction import *
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

def excentriciteNew(a,b):
    return sqrt(1 - (b**2/(a**2)))

def phi_from_beta(a,b,beta) :
    return atan(a/b * tan(beta))

def phi_from_psi(a,b,psi) :
    return atan((a/b)**2 * tan(psi))

#calcul de beta

def beta_from_phi(a,b,phi) :
    return atan(b/a * tan(phi))

def beta_from_psi(a,b,psi) :
    return atan(a/b * tan(psi))

#calcul de psi

def psi_from_phi(a,b,phi) :
    return atan((b/a)**2 * tan(phi))

def psi_from_beta(a,b,beta) :
    return atan((b/a) * tan(beta))


def tracer_latitude(phi,psi,beta):
    # Création d'une figure
    fig, ax = plt.subplots(subplot_kw={'polar': True})

    # Convertir les angles de degrés en radians
    phi_rad = phi * (pi / 180)  # Conversion en radians
    psi_rad = psi * (pi / 180)
    beta_rad = beta * (pi / 180)

    # Tracez les angles sur le cercle
    ax.plot([0, phi_rad], [0, 1], label='Phi', color='r', marker='o')
    ax.plot([0, psi_rad], [0, 1], label='Psi', color='g', marker='o')
    ax.plot([0, beta_rad], [0, 1], label='Beta', color='b', marker='o')

    ax.set_thetamin(-90)
    ax.set_thetamax(90)

    # Configurer les légendes
    ax.legend(loc='lower right')
    plt.title('Angles Phi, Psi et Beta sur un Cercle')
    # Afficher le graphique
    plt.show()
# tracer_latitude(10, 12,33)
#_____________________________________________________________________________________________________________

def surface(lambda_1,lambda_2,phi_1,phi_2,a,b):
    e = excentriciteNew(a,b)
    def f(phi):
        return abs((lambda_2-lambda_1)/2)*(b**2)*(((1/(2*e))*log((1+e*sin(phi))/(1-e*sin(phi))))+((sin(phi))/(1-(e**2)*(sin(phi))**2)))
        # (((lambda_1 - lambda_2)*b**2/2) * ((1/(2*e**2)*log((1+e*sin(phi))/(1-e*sin(phi)))) + ((sin(phi))/(1-(e**2) * (sin(phi))**2))))
    
    return (f(phi_2) - f(phi_1))

def surface_(l1 , l2 , phi, a ,b) :
    e = excentriciteNew(a,b)
    return (l2 - l1) * (b**2) * (sin(phi) + 2/3 * ((e**2) * (sin(phi)**3)) + 3/5 * ((e**4) * (sin(phi)**5)) + 4/7 * ((e**6) * (sin(phi)**7)))



#___________________________________________________________________Longuer d'arc____________________
def longeur_arc_meridien(phi , a, b):
    e2 =( excentriciteNew(a,b))**2
    Liste = [
        [1, 3 / 4, 45 / 64, 175 / 256, 11025 / 16348, 43659 / 65536],
        [3 / 4, 15 / 16, 525 / 512, 2205 / 2048, 72765 / 65536],
        [15 / 16, 105 / 256, 2205 / 4096, 10395 / 16348],
        [35 / 512, 315 / 2048, 31185 / 131072],
        [315 / 16384, 3465 / 65536],
        [639 / 131072],
    ]

    Le = [
        [1, e2, e2**2, e2**3, e2**4, e2**5],
        [e2, e2**2, e2**3, e2**4, e2**5],
        [e2**2, e2**3, e2**4, e2**5],
        [e2**3, e2**4, e2**5],
        [e2**4, e2**5],
        [e2**5],
    ]

    A_F = []
    for L1, L2 in zip(Liste, Le):
        coef = 0
        for i in range(len(L1)):
            coef += L1[i] * L2[i]
        A_F.append(coef)

    alpha_mue = [a * (1 - e2) * A_F[0]]
    for i in range(1, 6):
        alpha_mue.append(a * (1 - e2) * A_F[i] / (2 * i))

    S = alpha_mue[0] * phi
    for i in range(1, 6):
        S += sin(2 * i * phi) * alpha_mue[i] * (-1) ** i
    return abs(S)

def longeur_arc_parallele(phi , a , b ,lambda_1 , lambda_2) :
    e = excentriciteNew(a,b)
    W = sqrt(1-(e * sin(phi))**2)
    delta_lambda = lambda_2 - lambda_1 

    return a/W * cos(phi) * delta_lambda

def ArcMeridienDeuxParallele(a, b, phi1, phi2):
    Arc1 = longeur_arc_meridien(phi1 , a, b)
    Arc2 = longeur_arc_meridien(phi2 , a, b)
    if phi1 * phi2 > 0:
        if Arc1 > Arc2:
            return Arc1 - Arc2
        else:
            return Arc2 - Arc1
    else:
        return Arc1 + Arc2

#___Calcul des rayons______________________________________________________________________________

def rayon_2_demiAxe(a , b):
    return (2*a+b)/3

def rayon_Surface_sphere_Suraface_Elipse(E):
    # E est la surface d'ellipsoide
    return sqrt(E/(4*pi))


def rayon_Volume_sphere_EQUAL_Volume_Elipse(a , b):
    return (a**2 * b)**(1/3)

def rayon_GAUSS(a,b,phi):
    M = Rayon_courbure(a,b,phi)[1]
    N = Rayon_courbure(a,b,phi)[0]
    return sqrt(M*N)


#____Problème directe Sphère_________________________________________________________________________
def arcoTan(x):
    return atan(1/x)

def cotan(x):
    return 1/tan(x)

def pb_direct_sphere( phi1, lambda_1, D_12, A_12, R ):
    sigma_12 = D_12 / R #transformer la distance en rad

    def calcul_ph2(phi1 , A_12 , sigma_12):
        return asin( sin(phi1) * cos(sigma_12) + cos(phi1) * sin(sigma_12) * cos(A_12))

    def calcul_lambda_2(phi1 , A_12 , sigma_12):
        delta =  arcoTan((1/sin(A_12)) * ((cotan(sigma_12) * sin(pi/2 - phi1)) - (cos(pi/2 - phi1) * cos(A_12))))
        return lambda_1 + delta

    
    def calcul_A21(phi1 , A_12 , sigma_12):
        return  arcoTan((1/sin(A_12)) * ((cos(sigma_12) * cos(A_12)) - (tan(phi1) * sin(sigma_12))))

    phi2 = calcul_ph2(phi1 , A_12 , sigma_12)
    lambda_2 = calcul_lambda_2(phi1 , A_12 , sigma_12)
    A21 = calcul_A21(phi1 , A_12 , sigma_12)
    if A21 <=pi :
        A_21 = A21 + pi
    else :
        A_21 = A21 - pi
        

    return [phi2*180/pi , lambda_2*180/pi, A_21*180/pi]


def  pb_inverse_sphere (phi1 , phi2 , lambda_1 , lambda_2 , R):
    delta = lambda_2 - lambda_1
    def calcul_sigma_12(phi1 , phi2 , delta):
        return acos( sin(phi1)* sin(phi2) + cos(phi1) * cos(phi2) *  cos(delta))

    def calcul_A12(phi1 , phi2 , delta):
        # return atan(tan(phi2) * cos(phi1) / sin(delta) - sin(phi1) * cotan(delta))
        return asin(sin(delta) * cos(phi2)/ sin(calcul_sigma_12(phi1 , phi2 , delta)))

    def calcul_A21(phi1 , phi2 , delta):
        #return atan(-(tan(phi1) * cos(phi2) / sin(delta) - sin(phi2) * cotan(delta)))
        return asin(sin(delta) * cos(phi1)/ sin(calcul_sigma_12(phi1 , phi2 , delta)))

    sigma_12 = calcul_sigma_12(phi1 , phi2 , delta)
    A_12 = calcul_A12(phi1 , phi2 , delta)
    A21 = calcul_A21(phi1 , phi2 , delta)
    if A21 <=pi :
        A_21 = A21 + pi
    else :
        A_21 = A21 - pi

    return [sigma_12*R, A_12*180/pi , A_21*180/pi]



phi = 45 *pi/180  
phi1 = 35.05*pi/180
lambda_1 = -6.8*pi/180
D_12 = 10000
A12 = 30*pi/180
a = 6378249.145
b = 6356514.9658
M = Rayon_courbure(a,b,phi)[1]
N = Rayon_courbure(a,b,phi)[0]
R = rayon_GAUSS(a,b,phi)
phi2 = round(pb_direct_sphere( phi1, lambda_1, D_12, A12, R )[0] , 10)*pi/180
lambda_2 = round(pb_direct_sphere( phi1, lambda_1, D_12, A12, R )[1] , 10)*pi/180
# A21 = round(pb_direct_sphere( phi1, lambda_1, D_12, A12, R )[2] , 10)
# # print([phi2  , lambda2 , A21])
# print(pb_direct_sphere(phi1 , lambda_1, D_12, A12 , R))

# phi1 = 35.05*pi/180
# lambda_1 = -6.8*pi/180

# phi2 = pb_direct_sphere(phi1 , lambda_1, D_12, A12 , R)[0]*pi/180
# lambda_2 = pb_direct_sphere(phi1 , lambda_1, D_12, A12 , R)[1]*pi/180

# print( pb_inverse_sphere(phi1, phi2, lambda_1, lambda_2,R))


#______problème direct du ellipsoide formule puissante riguereuse et simplifier______________________________________________________

# formule Puissant : 

def pb_direct_ellisoide_puissant_general(Phi1, lam1, S, A12,a,b):
    e2 = excentricite(a,b)
    M1 = Rayon_courbure(a,b,Phi1)[1]
    N1 = Rayon_courbure(a,b,Phi1)[0]
    deltaPhiPrime = (
        S * cos(A12) / N1
        - (S / N1) ** 2 / 2 * tan(Phi1) * sin(A12) ** 2
        - (S / N1) ** 3 / 6 * cos(A12) * sin(A12) ** 2 * (1 + 3 * tan(Phi1) ** 2)
    )

    B = 1 / M1
    C = 3 / 2 * e2 * sin(Phi1) * cos(Phi1) / (1 - e2 * sin(Phi1) ** 2)
    D = tan(Phi1) / (2 * M1 * N1)
    E = (1 + 3 * tan(Phi1) ** 2) / (6 * N1**2)
    h = S * cos(A12) / M1
    sigmaPhi = N1 * deltaPhiPrime / M1

    deltaPhi = (
        S * cos(A12) * B
        - S**2 * sin(A12) ** 2 * D
        - h * S**2 * sin(A12) ** 2 * E
        - C * sigmaPhi**2
    )
    Phi2 = Phi1 + deltaPhi


    N2 = Rayon_courbure(a,b,Phi2)[0]
    deltaLam = (
        S
        / N2
        * sin(A12)
        / cos(Phi2)
        * (1 - 6 * (S / N2) ** 2 * (1 - (sin(A12) / cos(Phi2)) ** 2))
    )
    lam2 = lam1 + deltaLam

    PhiMoyenne = (Phi1 + Phi2) / 2
    deltaAzimut = 2 * arcoTan(cos(deltaPhi / 2) / sin(PhiMoyenne) * cotan(deltaLam / 2))

    # A21 = A12 + deltaAzimut + pi
    # if A21 >= 2 * pi:
    #     A21 -= 2 * pi
    # return Phi2 * 180 / pi, lam2 * 180 / pi, A21 * 180 / pi
    A_21 = A12 + deltaAzimut
    if A_21 <=pi :
        A21 = A_21 + pi
    else :
        A21 = A_21 - pi

    return [Phi2*180/pi , lam2*180/pi , A21*180/pi ]


def pb_direct_ellisoide_puissant_simplifier(a , b , S , phi1, lambda_1, A12 ):
    e = excentriciteNew(a,b)
    M1 = Rayon_courbure(a,b,phi1)[1]
    N1 =  Rayon_courbure(a,b,phi1)[0]
    B = 1/M1
    C = 3/2 * e**2 * sin(phi1) * cos(phi1)/(1-e**2 * sin(phi1)**2)
    D = tan(phi1)/(2*M1*N1)

    def delta_phi_prime(N1 , S , phi1 , A12):
        return S/N1 * cos(A12) - 0.5*(S/N1)**2 * tan(phi1) * sin(A12)**2 - (S/N1)**3 * cos(A12) * sin(A12)**2 * ( 1 + 3 * tan(phi1)**2)
    
    deta_phi = N1/M1 * delta_phi_prime(N1 , S , phi1 , A12)

    phi2 = deta_phi + phi1
    M2 = Rayon_courbure(a, b, phi2)[1]
    N2 = Rayon_courbure(a, b, phi2)[0]

    delta_phi = S * cos(A12)*B - S**2 * sin(A12)**2 * D - deta_phi**2 * C
    delta_lambda = S * sin(A12) / (cos(phi2) * N2)

    lambda_2 = delta_lambda + lambda_1

    delta_A = delta_lambda  * sin((phi1+phi2)/2)
    A_21 = (A12 + delta_A ) 
    if A_21 <=pi :
        A21 = A_21 + pi
    else :
        A21 = A_21 - pi
    return [phi2 *180/pi, lambda_2*180/pi , A21*180/pi]


def pb_direct_ellisoide_Gauss(Phi1, lam1, S, A12,a,b,p=10**(-20)):
    e2 = excentricite(a,b)
    deltaLam = []
    # étape 0
    deltaA = [0]
    deltaPhi = [0]
    Mm = Rayon_courbure(a,b,Phi1)[1]
    Nm = Rayon_courbure(a,b,Phi1)[0]
    Phim = Phi1

    deltaLam.append(S * sin(A12 + deltaA[-1] / 2) / (Nm * cos(Phim)))
    deltaA.append(2 * atan(tan(deltaLam[-1] / 2) * sin(Phim) / cos(deltaPhi[-1] / 2)))
    deltaPhi.append(S * cos(A12 + deltaA[-1] / 2) / (Mm * cos(deltaLam[-1] / 2)))

    Phi2 = [Phi1]
    lam2 = [lam1]
    A21 = [A12]

    Phi2.append(deltaPhi[-1] + Phi2[-1])
    lam2.append(deltaLam[-1] + lam2[-1])
    if deltaA[-1] + A21[-1] > pi:
        A21.append(deltaA[-1] + A21[-1] - pi)
    else:
        A21.append(deltaA[-1] + A21[-1] + pi)

    while (abs(Phi2[-2] - Phi2[-1]) > p or abs(lam2[-2] - lam2[-1]) > p or abs(A21[-2] - A21[-1]) > p ):
        Nm = (Rayon_courbure(a,b,Phi2[-1])[0] + Rayon_courbure(a,b,phi1)[0]) / 2
        Phim = (Phi2[-1] + Phi1) / 2
        Mm = (Rayon_courbure(a,b,Phi2[-1])[1] + Rayon_courbure(a,b,phi1)[1]) / 2

        deltaLam.append(S * sin(A12 + deltaA[-1] / 2) / (Nm * cos(Phim)))
        deltaA.append(
            2 * atan(tan(deltaLam[-1] / 2) * sin(Phim) / cos(deltaPhi[-1] / 2))
        )
        deltaPhi.append(S * cos(A12 + deltaA[-1] / 2) / (Mm * cos(deltaLam[-1] / 2)))

        Phi2.append(deltaPhi[-1] + Phi1)
        lam2.append(deltaLam[-1] + lam1)
        A21.append(deltaA[-1] + A12)
        A_21 = A21[-1]
        if A_21 <=pi :
            A_21 = A_21 + pi
        else :
            A_21 = A_21 - pi
        # if (deltaA[-1] + A12) >= pi:
        #     A21.append(deltaA[-1] + A12 - pi)
        # else:
        #     A21.append(deltaA[-1] + A12 + pi)

    return [Phi2[-1]*180/pi , lam2[-1]*180/pi , A_21*180/pi]
#problème inverse cad d'ellipsoide ____________________________________________________________________________________


def  pb_inverse_ellipsoide_Gauss (a, b , phi1 , phi2 , lambda1 , lambda2):

    Mm = (Rayon_courbure(a,b,phi1)[1] + Rayon_courbure(a,b,phi2)[1])/2
    Nm = (Rayon_courbure(a,b,phi1)[0] + Rayon_courbure(a,b,phi2)[0])/2

    phim = (phi1 + phi2)/2

    delta_lambda = lambda2 - lambda1
    delta_phi = phi2 - phi1

    delta_A =  2 * arcoTan(cos(delta_phi / 2) / sin(phim) * cotan(delta_lambda / 2))

    A12 = atan(cos(phim) * tan(delta_lambda / 2) / sin(Mm * delta_phi / (2 * Nm)))- delta_A / 2

    A_21 =  A12 + delta_A
    if A_21 <=pi :
            A_21 = A_21 + pi
    else :
        A_21 = A_21 - pi
    S  = delta_phi * Mm * cos(delta_lambda / 2) / cos(A12 + delta_A / 2)
    
    return [A12 *180/pi,  A_21*180/pi , S ]

# test du programme des problemes directes et inverses__________________________________________________________________________________________

phi1 = 35.05*pi/180
lambda_1 = -6.8*pi/180
S = 10000
A12 = 30*pi/180
a = 6378249.145
b = 6356514.9658
p = 10**(-20)
# print(pb_direct_ellisoide_puissant_general(phi1, lambda_1, S, A12,a,b))
# print(pb_direct_ellisoide_puissant_simplifier(a , b , S , phi1, lambda_1, A12 ))
# print(pb_direct_ellisoide_Gauss(phi1, lambda_1, S, A12 , a , b ,p = 10**(-20)))

# print([rayon_2_demiAxe(a,b) , rayon_GAUSS(a,b,45) , rayon_Surface_sphere_Suraface_Elipse(surface(0,pi,0,pi/2 , a,b)) , rayon_Volume_sphere_EQUAL_Volume_Elipse(a,b)])
# phi2_1 = pb_direct_ellisoide_puissant_general(a , b , S , phi1, lambda_1, A12 )[0] * pi/180
# phi2_2 = pb_direct_ellisoide_puissant_simplifier(a , b , S , phi1, lambda_1, A12 )[0] * pi/180
# phi2_3 = pb_direct_ellisoide_Gauss(a , b , S , phi1, lambda_1, A12 , p)[0] * pi/180

# lambda2_1 =  pb_direct_ellisoide_puissant_general(a , b , S , phi1, lambda_1, A12 )[1] * pi/180
# lambda2_2 = pb_direct_ellisoide_puissant_simplifier(a , b , S , phi1, lambda_1, A12 )[1] * pi/180
# lambda2_3 = pb_direct_ellisoide_Gauss(a , b , S , phi1, lambda_1, A12, p )[1] * pi/180

# print('probleme Direct')
# print(pb_direct_ellisoide_puissant_general(a , b , S , phi1, lambda_1, A12 ))
# print(pb_direct_ellisoide_puissant_simplifier(a , b , S , phi1, lambda_1, A12 ))
# print(pb_direct_ellisoide_Gauss(a , b , S , phi1, lambda_1, A12 , p))

# print('probleme inverse')
# print(pb_inverse_ellipsoide_Gauss(a, b , phi1 , phi2_1 , lambda_1 , lambda2_1))
# print(pb_inverse_ellipsoide_Gauss(a, b , phi1 , phi2_2 , lambda_1 , lambda2_2))
# print(pb_inverse_ellipsoide_Gauss(a, b , phi1 , phi2_3 , lambda_1 , lambda2_3))

def visMer(phi1 , phi2, lam = -6 ) :
    # Création de la carte basée sur un ellipsoïde
    fig = plt.figure(figsize=(8, 4))
    ax = plt.axes(projection=ccrs.PlateCarree())

    # Afficher le fond d'écran de la carte (l'ellipsoïde)
    ax.stock_img()

    # Coordonnées des points à matérialiser (latitude et longitude)
    points = [(phi1,lam), (phi2, lam)]

    # Matérialisation des points sur l'ellipsoïde
    for lat, lon in points:
        ax.plot(lon, lat, 'ro', markersize=2, transform=ccrs.PlateCarree())  # Placer un point rouge sur la carte
    lats, lons = zip(*points)
    ax.plot(lons, lats, 'b-', transform=ccrs.PlateCarree())  # Ligne bleue
    # Titre de la carte
    plt.title('Ellipsoïde avec des points matérialisés')

    # Afficher la carte
    plt.show()

import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QFrame
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

def plot_parallel(phi1, phi2, lam, ax):
    # Création de la carte basée sur un ellipsoïde
    ax.stock_img()

    # Tracer les parallèles de latitude phi1 et phi2
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', linestyle='--')
    
    gl.ylocator = plt.FixedLocator([phi1, phi2])  # Spécifier les parallèles de latitude

    # Ajouter des points rouges sur les parallèles
    ax.plot(0, phi1, 'ro', markersize=5, transform=ccrs.PlateCarree())
    ax.plot(0, phi2, 'ro', markersize=5, transform=ccrs.PlateCarree())

    # Tracer le méridien reliant les deux points
    meridian_points = [(phi1, lam), (phi2, lam)]
    lats, lons = zip(*meridian_points)
    ax.plot([0, 0], lats, 'g-', linewidth=2, transform=ccrs.PlateCarree())  # Ligne verte pour le méridien

    # Titre de la carte
    



def plot_meridian(lam1, lam2, phi, ax):
    # Création de la carte basée sur un ellipsoïde
    ax.stock_img()

    # Tracer les parallèles de latitude phi1 et phi2
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', linestyle='--')
    
    gl.xlocator = plt.FixedLocator([lam1, lam2])  # Spécifier les parallèles de latitude

    # Ajouter des points rouges sur les parallèles
    ax.plot(lam1,phi, 'ro', markersize=5, transform=ccrs.PlateCarree())
    ax.plot(lam2,phi ,'ro', markersize=5, transform=ccrs.PlateCarree())

    # Tracer le méridien reliant les deux points
    meridian_points = [(phi, lam1), (phi, lam2)]
    lats, lons = zip(*meridian_points)
    ax.plot(lons, [phi, phi], 'r-', linewidth=2, transform=ccrs.PlateCarree())  # Ligne red pour le méridien

    # Titre de la carte

import matplotlib.pyplot as plt
import numpy as np
from math import *

def VisAng(Phi , Beta ,Psi):
    # Paramètres de l'ellipse
    a = 3
    b = 1

    # # Angle φ en radians
    # Phi = pi / 3

    # # Calculer les angles β et ψ
    # Beta = beta_from_phi(a, b, Phi)
    # Psi = psi_from_phi(a, b, Phi)

    # Création des points de l'ellipse
    theta = np.linspace(-np.pi / 2, np.pi / 2, 10000)
    x_ellipse = a * np.cos(theta)
    y_ellipse = b * np.sin(theta)

    # Création des points pour l'angle φ
    phi_x = [0, a * np.cos(Phi)]
    phi_y = [0, b * np.sin(Phi)]

    # Création des points pour l'angle β
    beta_x = [0, a * np.cos(Beta)]
    beta_y = [0, b * np.sin(Beta)]

    # Création des points pour l'angle ψ
    psi_x = [0, a * np.cos(Psi)]
    psi_y = [0, b * np.sin(Psi)]

    # Dessiner l'ellipse
    plt.plot(x_ellipse, y_ellipse)

    # Dessiner le cercle (en prenant a comme rayon)
    plt.plot(x_ellipse, a * np.sin(theta))

    # Relier le centre de l'ellipse à l'angle β
    plt.plot([0, a * np.cos(Beta)], [0, b * np.sin(Beta)], color="g", linestyle="--")

    # Relier le centre de l'ellipse à l'angle ψ
    plt.plot([0, a * np.cos(Psi)], [0, b * np.sin(Psi)], color="r", linestyle="--")

    # Dessiner les angles φ, β, ψ
    plt.plot(phi_x, phi_y, color="b", label=f"Phi= {round(np.degrees(Phi) , 0)}°")
    plt.plot(beta_x, beta_y, color="g", label=f"Beta= {round(np.degrees(Beta) , 0)}°")
    plt.plot(psi_x, psi_y, color="r", label=f"Psi= {round(np.degrees(Psi)  , 0)}°")

    # Configurer le graphique
    plt.axhline(0, color="black", linewidth=0.5)
    plt.axvline(0, color="black", linewidth=0.5)
    plt.gca().set_aspect("equal", adjustable="box")
    plt.title("Cercle et Ellipse avec Angles")



    # Retirer les graduations des axes x et y
    plt.xticks([])
    plt.yticks([])

    # Légende
    plt.legend()
    plt.legend(loc='lower right')
    # Afficher le graphique
    plt.show()

