import sys
from PyQt5 import QtWidgets, uic
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from PyQt5.QtWidgets import *
from New_function import *


from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QFrame
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt

# fct_0 ----> GUI.ui  
# fct_UI_1 -----> UI_1
# fct_UI_2 -----> UI_2
class F3(QtWidgets.QMainWindow):
    import math
    def __init__(self , parent=None):
        import math
        super(F3 , self).__init__(parent)
        uic.loadUi('Projet Géodésie/UI_project/UI_longeurArc.ui' , self)
        self.show()
        self.setWindowTitle("Calcul des Longueurs")
        F3.setFixedSize(self , 580 , 320) #fixer les dimentions

        self.calcule_arcMer.clicked.connect(self.long_Mer)
        self.calculArcPara.clicked.connect(self.calcul_ArcPara)
        self.clear1.clicked.connect(self.remove)

    def long_Mer(self):
        try :
            a = float(str(eval(self.axeA.text())))
            b = float(str(eval(self.axeB.text())))
            phi1 = float(str(eval(self.phi1_mer.text())))* pi/180
            phi2 = float(str(eval(self.phi2_mer.text())))*pi/180
            longeurMer = round(ArcMeridienDeuxParallele(a,b,phi1, phi2) , 3)
            self.Rusult_ArcMeridian.setText(str(longeurMer))
            F3.setFixedSize(self , 376 , 120)
            # Créer un graphe avec Matplotlib
            fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
            plot_parallel(phi1 * 180/pi, phi2* 180/pi, -5.4, ax)
            F3.setFixedSize(self , 580 , 580)
            # Placer le graphe dans le QFrame
            canvas = FigureCanvas(fig)
            layout = QVBoxLayout(self.frame_4)
            layout.addWidget(canvas)
            self.frame_4.setLayout(layout)
        

        except :
            from PyQt5.QtWidgets import QMessageBox
            msgBox = QMessageBox()
            msgBox.setIcon(QMessageBox.Critical)
            msgBox.setWindowTitle("Erreur")
            msgBox.setText("Vérifier les paramètres saisis !")
            msgBox.exec_()


    def calcul_ArcPara(self):
        try : 
            a = float(str(eval(self.axeA.text())))
            b = float(str(eval(self.axeB.text())))
            phi = float(str(eval(self.phi_para.text())))* pi/180
            lambda_1 = float(str(eval(self.lambda1_para.text())))* pi/180
            lambda_2 = float(str(eval(self.lambda2_para.text())))* pi/180

            longeur = round(longeur_arc_parallele(phi , a , b ,lambda_1 , lambda_2 ) , 3)
            self.Result_ArcParallele.setText(str(longeur))
            
            # Créer un graphe avec Matplotlib
            fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
            plot_meridian(lambda_1 * 180/pi, lambda_2* 180/pi, phi* 180/pi, ax)
            F3.setFixedSize(self , 580 , 880)
            # Placer le graphe dans le QFrame
            canvas = FigureCanvas(fig)
            layout = QVBoxLayout(self.frame_5)
            layout.addWidget(canvas)
            self.frame_5.setLayout(layout)
        except : 
            from PyQt5.QtWidgets import QMessageBox
            msgBox = QMessageBox()
            msgBox.setIcon(QMessageBox.Critical)
            msgBox.setWindowTitle("Erreur")
            msgBox.setText("Vérifier les paramètres saisis !")
            msgBox.exec_()

    def remove(self):
        try: 
            F3.setFixedSize(self , 580 , 320)
            layout = self.frame_4.layout().deleteLater()
        except:
            from PyQt5.QtWidgets import QMessageBox
            msgBox = QMessageBox()
            msgBox.setIcon(QMessageBox.Critical)
            msgBox.setWindowTitle("Erreur")
            msgBox.setText("Saisis les angles !")
            msgBox.exec_()

        try: 
            F3.setFixedSize(self , 580 , 320)
            layout = self.frame_5.layout().deleteLater()
        except:
            from PyQt5.QtWidgets import QMessageBox
            msgBox = QMessageBox()
            msgBox.setIcon(QMessageBox.Critical)
            msgBox.setWindowTitle("Erreur")
            msgBox.setText("Saisis les angles !")
            msgBox.exec_()

        self.phi1_mer.clear()
        self.phi2_mer.clear()
        self.Rusult_ArcMeridian.clear()
        self.phi_para.clear()
        self.lambda1_para.clear()
        self.lambda2_para.clear()
        self.Result_ArcParallele.clear()

    
        
        
if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window  = F3()
    sys.exit(app.exec_())