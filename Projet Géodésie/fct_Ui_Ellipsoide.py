import sys
from PyQt5 import QtWidgets, uic
from PyQt5.QtWidgets import *
from fct_0 import *
from New_function import *

# fct_0 ----> GUI.ui  
# fct_UI_1 -----> UI_1
# fct_UI_2 -----> UI_2
class F6(QtWidgets.QMainWindow):
    def __init__(self ,parent=None):
        super(F6 , self).__init__(parent)
        uic.loadUi('Projet Géodésie/UI_project/Ui_pb_Ellipsoide.ui' , self)
        self.show()
        self.setWindowTitle("Problème sur Ellipsoide")
        F6.setFixedSize(self , 700 , 570) #fixer les dimentions

        self.calculer_direct.clicked.connect(self.cal_ellipsoide_direct)
        self.calculer_inverse.clicked.connect(self.cal_ellipsoide_inverse)
        self.clear_direct.clicked.connect(self.clear_dir)
        self.clear_inverse.clicked.connect(self.clear_inv)

    def cal_ellipsoide_direct(self):
        try : 
            a = float(str(eval(self.AxeA.text())))
            b = float(str(eval(self.AxeB.text())))
            phi1 = float(str(eval(self.phi1_direct.text())))* pi/180
            lambda1 = float(str(eval(self.lambda1_direct.text())))* pi/180
            D12 = float(str(eval(self.D_12_direct.text())))
            A12 = float(str(eval(self.A_12_direct.text())))*pi/180

            if self.puissant_simplifier.isChecked() == True:
                phi2 = round(pb_direct_ellisoide_puissant_simplifier( a, b , D12 , phi1, lambda1,  A12)[0] , 10)
                lambda2 = round(pb_direct_ellisoide_puissant_simplifier( a, b , D12 , phi1, lambda1,  A12)[1] , 10)
                A21 = round(pb_direct_ellisoide_puissant_simplifier( a, b , D12 , phi1, lambda1,  A12)[2] , 10)

            if self.Gauss.isChecked() == True:
                
                phi2 = round(pb_direct_ellisoide_Gauss(phi1, lambda1,D12,A12,a,b,p = 10**(-20))[0] , 10)
                lambda2 = round(pb_direct_ellisoide_Gauss(phi1, lambda1,D12,A12,a,b,p = 10**(-20))[1] , 10)
                A21 = round(pb_direct_ellisoide_Gauss(phi1, lambda1,D12,A12,a,b,p = 10**(-20))[2] , 10)
                print(phi2 , lambda2 , A21)
            
            if self.puissan_general.isChecked() == True:
                phi2 = round(pb_direct_ellisoide_puissant_general( phi1, lambda1,D12,A12,a,b)[0] , 10)
                lambda2 = round(pb_direct_ellisoide_puissant_general(  phi1, lambda1,D12,A12,a,b)[1] , 10)
                A21 = round(pb_direct_ellisoide_puissant_general(  phi1, lambda1,D12,A12,a,b)[2] , 10)

            if self.DvLimite.isChecked() == True:
                phi2 = round(pb_direct_ellisoide_puissant_general(  phi1, lambda1,D12,A12,a,b)[0] , 10)
                lambda2 = round(pb_direct_ellisoide_puissant_general(  phi1, lambda1,D12,A12,a,b)[1] , 10)
                A21 = round(pb_direct_ellisoide_puissant_general(  phi1, lambda1,D12,A12,a,b)[2] , 10)

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
    def cal_ellipsoide_inverse(self):
        try : 
            a = float(str(eval(self.AxeA.text())))
            b = float(str(eval(self.AxeB.text())))
            phi1 = float(str(eval(self.phi1_inverse.text())))* pi/180
            lambda1 = float(str(eval(self.lambda1_inverse.text())))* pi/180
            phi2 = float(str(eval(self.phi2_inverse.text())))* pi/180
            lambda2 = float(str(eval(self.lambda2_inverse.text())))* pi/180
    

            
            A12 = round(pb_inverse_ellipsoide_Gauss (a,b,phi1 , phi2 , lambda1 , lambda2)[0] , 10)
            A21 = round(pb_inverse_ellipsoide_Gauss (a,b,phi1 , phi2 , lambda1 , lambda2)[1] , 10)
            D12 = round(pb_inverse_ellipsoide_Gauss (a,b, phi1 , phi2 , lambda1 , lambda2)[2] , 10)

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




    
    def clear_dir(self):
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
    window  = F6()
    sys.exit(app.exec_())