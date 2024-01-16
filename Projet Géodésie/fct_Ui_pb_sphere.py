import sys
from PyQt5 import QtWidgets, uic
from PyQt5.QtWidgets import *
from fct_0 import *
from New_function import *


# fct_0 ----> GUI.ui  
# fct_UI_1 -----> UI_1
# fct_UI_2 -----> UI_2
class F5(QtWidgets.QMainWindow):
    import math
    def __init__(self , parent=None ):
        super(F5 , self).__init__(parent)
        uic.loadUi('Projet Géodésie/UI_project/UI_pb_sphere.ui' , self)
        self.show()
        self.setWindowTitle("Aproximation d'Ellpsoide par une sphère")
        F5.setFixedSize(self , 660 , 580) #fixer les dimentions


        self.calculer_direct.clicked.connect(self.caluler_Dr)
        self.clear_direct.clicked.connect(self.clear_Dr)

        self.calculer_inverse.clicked.connect(self.caluler_inv)
        self.clear_inverse.clicked.connect(self.clear_inv)
    
    def caluler_Dr(self):
        try : 
            a = float(str(eval(self.axeA.text())))
            b = float(str(eval(self.axeB.text())))
            phi1 = float(str(eval(self.phi1_direct.text())))* pi/180
            lambda1 = float(str(eval(self.lambda1_direct.text())))* pi/180
            D12 = float(str(eval(self.D_12_direct.text())))
            A12 = float(str(eval(self.A_12_direct.text())))*pi/180

            if (self.rayon_moyen.isChecked()== True):
                R = rayon_2_demiAxe(a,b)
            elif (self.rayon_surface.isChecked()== True):
                E = surface(0,pi,0,pi/2,a,b) * 4
                R = rayon_Surface_sphere_Suraface_Elipse(E)
            elif (self.rayo_Gauss.isChecked()== True):
                phi = float(str(eval(self.phi_gauss.text())))*pi/180
                R = rayon_GAUSS(a , b , phi)
            elif (self.rayon_volume.isChecked()== True):
                R = rayon_Volume_sphere_EQUAL_Volume_Elipse(a,b)
            
            phi2 = round(pb_direct_sphere( phi1, lambda1, D12, A12, R )[0] , 10)
            lambda2 = round(pb_direct_sphere( phi1, lambda1, D12, A12, R )[1] , 10)
            A21 = round(pb_direct_sphere( phi1, lambda1, D12, A12, R )[2] , 10)

            self.phi2_result_direct.setText(str(phi2))
            self.lambda2_result_direct.setText(str(lambda2))
            self.azimut_retour_direct.setText(str(A21))
        except :
            from PyQt5.QtWidgets import QMessageBox
            msgBox = QMessageBox()
            msgBox.setIcon(QMessageBox.Critical)
            msgBox.setWindowTitle("Erreur")
            msgBox.setText("Vérifier les paramètres saisis !")
            msgBox.exec_()



    def caluler_inv(self):
        try : 
            a = float(str(eval(self.axeA.text())))
            b = float(str(eval(self.axeB.text())))
            phi1 = float(str(eval(self.phi1_inverse.text())))* pi/180
            lambda1 = float(str(eval(self.lambda1_inverse.text())))* pi/180
            phi2 = float(str(eval(self.phi2_inverse.text())))* pi/180
            lambda2 = float(str(eval(self.lambda2_inverse.text())))* pi/180
            
            if (self.rayon_moyen.isChecked()== True):
                R = rayon_2_demiAxe(a,b)
            elif (self.rayon_surface.isChecked()== True):
                E = surface_(0,pi,pi/2,a,b) * 4
                R = rayon_Surface_sphere_Suraface_Elipse(E)
                print(E)
            elif (self.rayo_Gauss.isChecked()== True):
                phi = float(str(eval(self.phi_gauss.text())))*pi/180
                R = rayon_GAUSS(a , b , phi)
            elif (self.rayon_volume.isChecked()== True):
                R = rayon_Volume_sphere_EQUAL_Volume_Elipse(a,b)

            D12 = round(pb_inverse_sphere (phi1 , phi2 , lambda1 , lambda2 , R)[0], 3)
            A12 = round(pb_inverse_sphere (phi1 , phi2 , lambda1 , lambda2 , R)[1] , 9)
            A21 = round(pb_inverse_sphere (phi1 , phi2 , lambda1 , lambda2 , R)[2], 9)

            self.azimut_aller_inverse.setText(str(A12))
            self.azimut_retour_inverse.setText(str(A21))
            self.D_inverse.setText(str(D12))
        except:
            from PyQt5.QtWidgets import QMessageBox
            msgBox = QMessageBox()
            msgBox.setIcon(QMessageBox.Critical)
            msgBox.setWindowTitle("Erreur")
            msgBox.setText("Vérifier les paramètres saisis !")
            msgBox.exec_()


    def clear_Dr(self):
        self.phi1_direct.clear()
        self.lambda1_direct.clear()
        self.D_12_direct.clear()
        self.A_12_direct.clear()
        self.phi2_result_direct.clear()
        self.lambda2_result_direct.clear()
        self.azimut_retour_direct.clear()

    def clear_inv(self):
        self.phi1_inverse.clear()
        self.lambda1_inverse.clear()
        self.phi2_inverse.clear()
        self.lambda2_inverse.clear()
        self.D_inverse.clear()
        self.azimut_aller_inverse.clear()
        self.azimut_retour_inverse.clear()

    

        
if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window  = F5()
    sys.exit(app.exec_())