import sys
from PyQt5 import QtWidgets, uic
from PyQt5.QtWidgets import *
from New_function import *


# fct_0 ----> GUI.ui  
# fct_UI_1 -----> UI_1
# fct_UI_2 -----> UI_2
class F4(QtWidgets.QMainWindow):
    import math
    def __init__(self , parent=None):
        import math
        super(F4 , self).__init__(parent)
        uic.loadUi('Projet Géodésie/UI_project/UI_surface.ui' , self)
        self.show()
        self.setWindowTitle("Calcul des surfaces")
        F4.setFixedSize(self , 720 , 200) #fixer les dimentions

        self.calculSurface.clicked.connect(self.calcul_Surface)
        self.clear2.clicked.connect(self.remove)

    def calcul_Surface(self):
        try : 
            a = float(str(eval(self.axeA.text())))
            b = float(str(eval(self.axeB.text())))
            phi1 = float(str(eval(self.phi1.text())))* pi/180
            phi2 = float(str(eval(self.phi2.text())))*pi/180
            lambda1 = float(str(eval(self.lambda1.text())))* pi/180
            lambda2 = float(str(eval(self.lambda2.text())))*pi/180

            S = surface(lambda1 , lambda2 , phi1, phi2 ,a,b)

            self.Rusult_surface.setText(str(S))
        except:
            from PyQt5.QtWidgets import QMessageBox
            msgBox = QMessageBox()
            msgBox.setIcon(QMessageBox.Critical)
            msgBox.setWindowTitle("Erreur")
            msgBox.setText("Vérifier les paramètres saisis !")
            msgBox.exec_()
    
    def remove(self) :
        self.axeA.clear()
        self.axeB.clear()
        self.phi1.clear()
        self.phi2.clear()
        self.Rusult_surface.clear()
        self.lambda1.clear()
        self.lambda2.clear()
        
if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window  = F4()
    sys.exit(app.exec_())