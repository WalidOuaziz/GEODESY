from math import*
#calcul_parametres() :
def aplatissement(a,b):
    return (1/(1 - b/a))

def excentricite(a,b):
    return (1 - (b**2/(a**2)))

def grand_axe1(b,inv_f):
    return b/(1-1/inv_f)

def petit_axe1(a,inv_f):
    return (1-1/inv_f)*a

def petit_axe2(a,e2):
    return a*sqrt(1-e2)

def grand_axe2(b,e2):
    return b/sqrt(1-e2)

def dms(angle):
    angle_en_degrees = degrees(angle)
    deg = int(angle_en_degrees)
    decimal_part_degrees, integer_part_minutes = modf(angle_en_degrees)
    minutes = int(decimal_part_degrees * 60)
    decimal_part_minutes, integer_part_seconds = modf(decimal_part_degrees * 60)
    seconds = decimal_part_minutes * 60
    return deg, minutes, round(seconds,8)

  
# Calcule de rayn de courbure :-----------------------------------------
def Rayon_courbure(a,b,phi):
    e2 = excentricite(a,b)
    W = sqrt(1-(e2) * (sin(phi)**2))
    N = a/W
    M = a*(1-e2)/W**3
    R = [N , M]
    return R

#Calcule de la latitude isomètrique :--------------------------------------------------

def U(a,b,phi):
    # phi  = phi * pi/200
    e2 = excentricite(a,b)
    u = log(tan(pi/4 + phi/2))- (((e2)/2) * log((1 + sqrt(e2)*sin(phi))/(1 - sqrt(e2)*sin(phi))))
    return u

#Calcule de phi a partir de U :
def calculer_phi(a,b,u,phi_0, P ):
    e2 = excentricite(a, b)
    phi_0 = phi_0 * pi/200
    E = (e2/2)*(log((1+sqrt(e2)*sin(phi_0))/(1-sqrt(e2)*sin(phi_0))))
    phi = 2*atan((exp(u+E)-1)/(exp(u+E)+1))
    liste_phi = [phi_0 , phi]
    while abs(liste_phi[-1] - liste_phi[-2])>= P:
        E = (e2/2)*(log((1+sqrt(e2)*sin(liste_phi[-1]))/(1-sqrt(e2)*sin(liste_phi[-1]))))
        phi = 2*atan((exp(u+E)-1)/(exp(u+E)+1))
        liste_phi.append(phi)

    return  liste_phi[-1]

# -------------------------------------------------------
#Coord geo =====> coord Cartographique de LAMBERt

def geo_vers_Lambert(phi , L , Zone):
    # parametres de l'ellipsoide de Clark 1880
    phi = phi*pi/180
    L = L * pi/180
    a = 6378249.145
    b = 6356514.966
    u = U(a , b , phi)

    if Zone == 1 :
        X0 = 500000
        Y0 = 300000
        L0 = -6 * pi/200
        phi_0 = 37*pi/200
        K0 = 0.999625769
        u0 = U(a , b , phi_0)
        N0 = Rayon_courbure(a,b,phi_0)[0]
        e0 = K0 * N0 / (tan(phi_0))
        dU = (u - u0)
        dL = (L - L0)
        X = e0 * exp(-dU * sin(phi_0)) * sin(sin(phi_0)*dL) + X0
        Y = -e0 * (exp(-dU * sin(phi_0)) * cos(sin(phi_0)*dL) - 1) + Y0

    if Zone == 2 :
        X0 = 500000
        Y0 = 300000
        L0 = -6 * pi/200
        phi_0 = 33*pi/200
        K0 = 0.999615596
        u0 = U(a , b , phi_0)
        N0 = Rayon_courbure(a,b,phi_0)[0]
        e0 = K0 * N0 / (tan(phi_0))
        dU = (u - u0)
        dL = (L - L0)
        X = e0 * exp(-dU * sin(phi_0)) * sin(sin(phi_0)*dL) + X0
        Y = -e0 * (exp(-dU * sin(phi_0)) * cos(sin(phi_0)*dL) - 1) + Y0

    if Zone == 3 :
        X0 = 1200000
        Y0 = 400000
        L0 = -6 * pi/200
        phi_0 = 29*pi/200
        K0 = 0.999616304
        u0 = U(a , b , phi_0)
        N0 = Rayon_courbure(a,b,phi_0)[0]
        e0 = K0 * N0 / (tan(phi_0))
        dU = (u - u0)
        dL = (L - L0)
        X = e0 * exp(-dU * sin(phi_0)) * sin(sin(phi_0)*dL) + X0
        Y = -e0 * (exp(-dU * sin(phi_0)) * cos(sin(phi_0)*dL) - 1) + Y0

    if Zone == 4 :
        X0 = 1500000
        Y0 = 400000
        L0 = -6 * pi/200
        phi_0 = 25*pi/200
        K0 = 0.999616437
        u0 = U(a , b , phi_0)
        N0 = Rayon_courbure(a,b,phi_0)[0]
        e0 = K0 * N0 / (tan(phi_0))
        dU = (u - u0)
        dL = (L - L0)
        X = e0 * exp(-dU * sin(phi_0)) * sin(sin(phi_0)*dL) + X0
        Y = -e0 * (exp(-dU * sin(phi_0)) * cos(sin(phi_0)*dL) - 1) + Y0
    return round(X,3),round(Y,3)

# transformation de coord Lambret vers coord Geo :--------------------------------------------------------------------------
def coord_Lambert_vers_Geo(X,Y,Zone):
    a = 6378249.2000
    b = 6356515.0000
    e2 = excentricite(a, b)

    if Zone == 1 :
        X = X - 500000
        Y = Y - 300000
        L0 = -6 * pi/200
        phi_0 = 37*pi/200
        K0 = 0.999625769
  

    if Zone == 2 :
        X = (X - 500000)
        Y = Y - 300000
        L0 = -6 * pi/200
        phi_0 = 33*pi/200
        K0 = 0.999615596

    if Zone == 3 :
        X = X - 1200000
        Y = Y - 400000
        L0 = -6 * pi/200
        phi_0 = 29*pi/200
        K0 = 0.999616304

    if Zone == 4 :
        X = X - 1500000
        Y = Y - 400000
        L0 = -6 * pi/200
        phi_0 = 25*pi/200
        K0 =  0.999616437 


    u0 = U(a , b , phi_0)
    N0 = Rayon_courbure(a,b,phi_0)[0]
    e0 = K0 * N0 / tan(phi_0)

    # longitude Lambda
    L_cal = (1/(sin(phi_0))) * (atan(X/(e0 - Y)))
    L = L_cal + L0


    # Latitude phi
    du = (1/(2*sin(phi_0))) * log(e0**2/(X**2 + (e0 - Y)**2))
    u = abs(du + u0)
    # E = log(tan((pi/4) + (phi_0/2))) - u

    E = (e2/2)*(log((1+sqrt(e2)*sin(phi_0))/(1-sqrt(e2)*sin(phi_0))))

    phi = 2*atan((exp(u+E)-1)/(exp(u+E)+1))
    liste_phi = [phi_0 , phi]
    # print(liste_phi[0])
    # print(liste_phi[-1]*200/pi)

    while abs(liste_phi[-1] - liste_phi[-2])>= 10**(-10):
        E = (e2/2)*(log((1+sqrt(e2)*sin(liste_phi[-1]))/(1-sqrt(e2)*sin(liste_phi[-1]))))
        phi = 2*atan((exp(u+E)-1)/(exp(u+E)+1))
        liste_phi.append(phi)

    return  liste_phi[-1], L

# conversion entre les unités des angles
def angle_rad(A):
    A_deg_dec = A*180/pi
    A_dms = dms(A)
    A_gon = A*200/pi
    return round(A_deg_dec , 5) , round(A_gon,5), A_dms 

def angle_deg_dec(A):
    A_rad = A*pi/180
    A_gon = A*200/180
    A_dms = dms(A_rad)
    return A_rad , round(A_gon,5) , A_dms

def angle_gon(A):
    A_rad = A*pi/200
    A_deg_dec = A*180/200
    A_dms = dms(A_rad)
    return A_rad ,round(A_deg_dec,5) , A_dms

def angle_dms(A_deg , A_min ,A_sec):
    A_rad = (A_deg + A_min/60 + A_sec/3600)*pi/180
    A_deg_dec = A_deg + A_min/60 + A_sec/3600
    A_gon = (A_deg + A_min/60 + A_sec/3600) *200/180
    return A_rad , round(A_deg_dec , 5) , round(A_gon,5)


# def G_L(phi , L):
    # # parametres de l'ellipsoide de Clark 1880
    # phi = phi*pi/200
    # L = L * pi/200
    # a = 6378249.145
    # b = 6356514.966
    # u = U(a , b , phi)

    # # if Zone == 1 :
    # X0 = 500000
    # Y0 = 300000
    # L0 = -6 * pi/200
    # phi_0 = 37*pi/200
    # K0 = 0.999625769
    # u0 = U(a , b , phi_0)
    # N0 = Rayon_courbure(a,b,phi_0)[0]
    # # ------
    # e0 = K0 * N0 / (tan(phi_0))
    # dU = (u - u0) 
    # dL = (L - L0)
    # X = e0 * exp(-dU * sin(phi_0)) * sin(sin(phi_0)*dL) + X0
    # Y = -e0 * (exp(-dU * sin(phi_0)) * cos(sin(phi_0)*dL) - 1) + Y0
    # # cood Polaire
    # gama = (L-L0)*sin(phi_0)
    # R = e0*exp(-sin(phi_0)*dU)
    # # Coord Lambert 
    # X = X0 + R*sin(gama)
    # Y = Y0 + e0 - R*cos(gama)
    # return X , Y



    