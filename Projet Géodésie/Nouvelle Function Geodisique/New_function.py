from math import *
from fonction import *

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
    return atan((b/a)**2 * tan(beta))

# a = 6378249.145
# b = 6356514.9658
# phi = pi/2
# psi = pi/3


# beta = pi/4


# Test ------
# print('phi a partir de beta est : ', phi_from_beta(a,b,beta))
# print('psi a partir de phi est : ',  psi_from_phi(a,b,phi_from_beta(a,b,beta)))
# print('beta a partir de psi est : ', beta_from_psi(a,b,psi_from_phi(a,b,phi_from_beta(a,b,beta))))
# print(pi/4)

#------------Ajouter la representation des 3 angles



# import matplotlib.pyplot as plt

# # Les angles en radians
# phi = 30  # Angle phi en degrés
# psi = 45  # Angle psi en degrés
# beta = 60  # Angle beta en degrés

# # Création d'une figure
# fig, ax = plt.subplots(subplot_kw={'polar': True})

# # Convertir les angles de degrés en radians
# phi_rad = phi * (3.14159265359 / 180)  # Conversion en radians
# psi_rad = psi * (3.14159265359 / 180)
# beta_rad = beta * (3.14159265359 / 180)

# # Tracez les angles sur le cercle
# ax.plot([0, phi_rad], [0, 1], label='Phi', color='r', marker='o')
# ax.plot([0, psi_rad], [0, 1], label='Psi', color='g', marker='o')
# ax.plot([0, beta_rad], [0, 1], label='Beta', color='b', marker='o')

# ax.set_thetamin(-90)
# ax.set_thetamax(90)

# # Configurer les légendes
# ax.legend(loc='lower right')
# plt.title('Angles Phi, Psi et Beta sur un Cercle')

# # Afficher le graphique
# plt.show()

#_____________________________________________________________________________________________________________


def surface(lambda_1,lambda_2,phi_1,phi_2,a,b):
    e = excentriciteNew(a,b)
    def f(phi):
        return ((1/2*e) * log((1+e*sin(phi))/(1-e*sin(phi)))) + ((sin(phi))/(1-(e**2) * (sin(phi))**2))
    
    return 1/2 * (lambda_2 - lambda_1) * (b**2) * (f(phi_2) - f(phi_1))




def surface_(l1 , l2 , phi, a ,b) :
    e = excentriciteNew(a,b)
    return 1/2 * (l2 - l1) * (b**2) * (sin(phi) + 2/3 * ((e**2) * (sin(phi)**3)) + 3/5 * ((e**4) * (sin(phi)**5)) + 4/7 * ((e**6) * (sin(phi)**7)))


# l1 = 0
# l2 = pi
# phi1 = 0
# phi2 = pi/2
# a = 6378249.145
# b = 6356514.9658


# s1 = surface_2(l1 , l2 , phi2 ,  a ,b)

# s2 = surface(l1 , l2 , phi1 , phi2 , a ,b)

# print(s1)
# print(s2)
# print(s1 - s2)


#___________________________________________________________________Longuer d'arc____________________


def longeur_arc_meridien (phi , a , b) :
    e = excentriciteNew(a,b)
    w = a*(1-e**2)

    A = sum([i*j for i,j in zip([1,3/4 , 45/64 , 175/256 , 11025/16348 , 43659/65536], [1 , e**2 ,e**4,e**6,e**8,e**10 ])])
    B = sum([i*j for i,j in zip([3/4 , 15/16, 525/512 , 2205/2048 , 72765/65536], [ e**2 , e**4, e**6, e**8, e**10  ])])
    C = sum([i*j for i,j in zip([15/16, 105/256, 2205/4096, 10395/16348], [ e**4, e**6, e**8, e**10  ])])
    D = sum([i*j for i,j in zip([35/512 , 315/2048, 31185/131072], [e**6, e**8, e**10 ])])
    E = sum([i*j for i,j in zip([315/16348 , 3465/65536], [e**8, e**10  ])])
    F = 639/131072  *  e**10

    List_1 = [A,B,C,D,E,F]
    List_2 = [w,w/2 , w/4, w/6 , w/8 , w/10]

    List_3 = [i*j for i,j in zip(List_1, List_2)]              #alpha , beta..........

    List_sin = [phi , sin(2*phi),sin(4*phi),sin(6*phi),sin(8*phi),sin(10*phi)]

    L = List_3[0]*List_sin[0] - List_3[1]*List_sin[1] + List_3[2]*List_sin[2] - List_3[3]*List_sin[3] + List_3[4]*List_sin[4]-List_3[5]*List_sin[5]

    return L



def longeur_arc_parallele(phi , a , b ,lambda_1 , lambda_2) :
    e = excentriciteNew(a,b)
    W = sqrt(1-(e * sin(phi))**2)
    delta_lambda = lambda_2 - lambda_1 

    return a/W * cos(phi) * delta_lambda




    
# a = 6378249.145
# b = 6356514.9658
# phi = pi/2

# l1 = 0
# l2 = pi

# print(longeur_arc_meridien(phi , a , b))

# print(longeur_arc_parallele(phi , a , b ,l1 , l2 ))
#___Calcul des rayons______________________________________________________________________________

def rayon_2_demiAxe(a , b):
    return (2*a+b)/3

def rayon_Surface_sphere_Suraface_Elipse(E):
    # E est la surface d'ellipsoide
    return sqrt(E/(4*pi))


def rayon_Volume_sphere_EQUAL_Volume_Elipse(a , b):
    return (a**2 * b)**(1/3)

def rayon_GAUSS(M , N):
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
    A_21 = calcul_A21(phi1 , A_12 , sigma_12)

    return [phi2*180/pi , lambda_2*180/pi, 360-A_21*180/pi]


def  pb_inverse_sphere (phi1 , phi2 , lambda_1 , lambda_2):
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
    A_21 = calcul_A21(phi1 , phi2 , delta)

    return [sigma_12*R, A_12*180/pi , 360 - A_21*180/pi]


#test_________________________________________________

# phi = 45 *pi/180  
# phi1 = 35.05*pi/180
# lambda_1 = -6.8*pi/180
# D_12 = 10000
# A12 = 30*pi/180
# a = 6378249.145
# b = 6356514.9658
# M = Rayon_courbure(a,b,phi)[1]
# N = Rayon_courbure(a,b,phi)[0]
# R = rayon_GAUSS(M,N)

# print(pb_direct_sphere(phi1 , lambda_1, D_12, A12 , R))

# phi1 = 35.05*pi/180
# lambda_1 = -6.8*pi/180

# phi2 = pb_direct_sphere(phi1 , lambda_1, D_12, A12 , R)[0]*pi/180
# lambda_2 = pb_direct_sphere(phi1 , lambda_1, D_12, A12 , R)[1]*pi/180

# print( pb_inverse_sphere(phi1, phi2, lambda_1, lambda_2))


#______problème direct du ellipsoide formule puissante riguereuse et simplifier______________________________________________________

# formule Puissant : 

def pb_direct_ellisoide_puissant_general(a , b , S , phi1, lambda_1, A12 ):
    e = excentriciteNew(a,b)
    M1 = Rayon_courbure(a,b,phi1)[1]
    N1 = Rayon_courbure(a,b,phi1)[0]
    B = 1/M1
    C = 3/2 * e**2 * sin(phi1) * cos(phi1) / (1-e**2 * sin(phi1)**2)
    D = tan(phi1)/(2 * M1 * N1)
    E = (1 + 3*tan(phi1)**2)/(6 * N1**2)
    H = (S/M1) * cos(A12)
    
    def delta_phi_prime(N1 , S , phi1 , A12):
        return S/N1 * cos(A12) - 0.5*(S/N1)**2 * tan(phi1) * sin(A12)**2 - (S/N1)**3 * cos(A12) * sin(A12)**2 * ( 1 + 3 * tan(phi1)**2)

    deta_phi = N1/M1 * delta_phi_prime(N1 , S , phi1 , A12)

    def delta_phi(N1 , S , phi1 , A12):
        return S*B*cos(A12) - S**2 * D *sin(A12)**2 - H * S**2 * sin(A12)**2 * E - C * deta_phi**2

    def calcul_phi2(N1 , S , phi1 , A12):
        return phi1 + delta_phi(N1 , S , phi1 , A12)

    phi2 = calcul_phi2(N1 , S , phi1 , A12)

    M2 = Rayon_courbure(a, b, phi2)[1]
    N2 = Rayon_courbure(a, b, phi2)[0]

    def sec(x):
        return 1/cos(x)

    def calcul_lambda_2(N2 , S , phi2, A12 ):
        delta_lambda =  S/N2 * sin(A12) * sec(phi2) * (1 - (1/6) * (S/N2**2) * (1 - sin(A12)**2 * sec(phi2)**2))
        return delta_lambda + lambda_1

    lambda_2 = calcul_lambda_2(N2 , S , phi2, A12 )

    def calcul_A21(A12 , phi1 , phi2 , lambda_1 , lambda_2):
        phi_moy =( phi1 + phi2)/2
        delta_phi = phi2  -phi1
        delta_lambda = lambda_1 - lambda_2

        def delta_A (phi_moy , delta_phi , delta_lambda):
            return 2 * atan(sin(phi_moy) /(cos(delta_phi) * cotan(delta_lambda/2)))


        return  A12 + delta_A(phi_moy , delta_phi , delta_lambda)

    A21 = calcul_A21(A12 , phi1 , phi2 , lambda_1 , lambda_2)

    return [phi2* 180/pi , lambda_2 * 180/pi , 360-(A21* 180/pi)]

    

def pb_direct_ellisoide_puissant_simplifier(a , b , S , phi1, lambda_1, A12 ):
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
    A21 = (A12 + delta_A ) 
    return [phi2 *180/pi, lambda_2*180/pi , 360-(A21*180/pi) ]


#______problème direct du ellipsoide methode de Gauss______________________________________________________


def pb_direct_ellisoide_Gauss(a , b , S , phi1, lambda_1, A12, p ):
    
    M1 = Rayon_courbure(a,b,phi1)[1]
    N1 =  Rayon_courbure(a,b,phi1)[0]

    Nm = [N1]
    Mm = [M1]
    Phi_m = [phi1]


    delta_A = [0]
    delta_phi = [0]
    

    
    delta_lambda = [(S * sin(A12 + delta_A[-1]/2 )) / (Nm[-1] * cos(Phi_m [-1]))]

    delta_A.append(2 * atan (tan  ( delta_lambda[-1]/2) * sin(Phi_m[-1]) / cos(delta_phi[-1]) ))

    delta_phi.append(( S/Mm[-1]) * cos(A12 + delta_A[-1]/2)/cos(delta_lambda[-1]/2))

    phi2_0 = delta_phi[-1] + phi1
    lambda2_0 = delta_lambda[-1] + lambda_1

    N2 = Rayon_courbure(a,b , phi2_0)[0]
    M2 = Rayon_courbure(a,b , phi2_0)[1]

    Mm.append( (M1 + M2)/2)
    Nm.append( (N1 + N2)/2)
    Phi_m.append((phi1 + phi2_0)/2)

    
    delta_lambda.append((S * sin(A12 + delta_A[-1]/2 )) / (Nm[-1] * cos(Phi_m [-1])))

    list_phi2 = [0 , phi2_0]
    list_lambda = [0 , lambda2_0]
    list_A21 = [0 , delta_A[-1] + A12]

    while (abs(list_A21[-1] - list_A21[-2]) > p or abs(list_phi2[-1] - list_phi2[-2]) > p  or abs(list_lambda[-1] - list_lambda[-2]) > p ):
       
        delta_A.append(2 * atan (tan  ( delta_lambda[-1]) * sin(Phi_m[-1]) / cos(delta_phi[-1]) ))
        delta_phi.append(( S/Mm[-1]) * cos(A12 + delta_A[-1]/2)/cos(delta_lambda[-1]/2))

        list_A21.append(delta_A[-1] + A12)

        phi2_0 = delta_phi[-1] + phi1
        list_phi2.append(phi2_0)
        lambda2_0 = delta_lambda[-1] + lambda_1
        list_lambda.append(lambda2_0)

        N2 = Rayon_courbure(a,b , phi2_0)[0]
        M2 = Rayon_courbure(a,b , phi2_0)[1]

        Mm.append( (M1 + M2)/2)
        Nm.append( (N1 + N2)/2)
        Phi_m.append((phi1 + phi2_0)/2)
        delta_lambda.append((S * sin(A12 + delta_A[-1]/2 )) / (Nm[-1] * cos(Phi_m [-1])))
    return [phi2_0 *180/pi, lambda2_0 *180/pi , 360-(delta_A[-1]+A12) *180/pi ]
    
#problème inverse cad d'ellipsoide ____________________________________________________________________________________


def  pb_inverse_ellipsoide_Gauss (a, b , phi1 , phi2 , lambda1 , lambda2):

    Mm = (Rayon_courbure(a,b,phi1)[1] + Rayon_courbure(a,b,phi2)[1])/2
    Nm = (Rayon_courbure(a,b,phi1)[0] + Rayon_courbure(a,b,phi2)[0])/2

    phim = (phi1 + phi2)/2

    delta_lambda = lambda2 - lambda1
    delta_phi = phi2 - phi1

    delta_A =  2 * arcoTan(cos(delta_phi / 2) / sin(phim) * cotan(delta_lambda / 2))

    A12 = atan(cos(phim) * tan(delta_lambda / 2) / sin(Mm * delta_phi / (2 * Nm)))- delta_A / 2

    A21 = 2*pi - A12 + delta_A 

    S  = delta_phi * Mm * cos(delta_lambda / 2) / cos(A12 + delta_A / 2)
    
    return [A12 *180/pi, A21 *180/pi , S ]

# test du programme des problemes directes et inverses__________________________________________________________________________________________

phi1 = 35.05*pi/180
lambda_1 = -6.8*pi/180
S = 10000
A12 = 30*pi/180
a = 6378249.145
b = 6356514.9658
p = 10**(-20)

phi2_1 = pb_direct_ellisoide_puissant_general(a , b , S , phi1, lambda_1, A12 )[0] * pi/180
phi2_2 = pb_direct_ellisoide_puissant_simplifier(a , b , S , phi1, lambda_1, A12 )[0] * pi/180
phi2_3 = pb_direct_ellisoide_Gauss(a , b , S , phi1, lambda_1, A12 , p)[0] * pi/180

lambda2_1 =  pb_direct_ellisoide_puissant_general(a , b , S , phi1, lambda_1, A12 )[1] * pi/180
lambda2_2 = pb_direct_ellisoide_puissant_simplifier(a , b , S , phi1, lambda_1, A12 )[1] * pi/180
lambda2_3 = pb_direct_ellisoide_Gauss(a , b , S , phi1, lambda_1, A12, p )[1] * pi/180


print(pb_direct_ellisoide_puissant_general(a , b , S , phi1, lambda_1, A12 ))
print(pb_direct_ellisoide_puissant_simplifier(a , b , S , phi1, lambda_1, A12 ))
print(pb_direct_ellisoide_Gauss(a , b , S , phi1, lambda_1, A12 , p))


print(pb_inverse_ellipsoide_Gauss(a, b , phi1 , phi2_1 , lambda_1 , lambda2_1))
print(pb_inverse_ellipsoide_Gauss(a, b , phi1 , phi2_2 , lambda_1 , lambda2_2))
print(pb_inverse_ellipsoide_Gauss(a, b , phi1 , phi2_3 , lambda_1 , lambda2_3))
