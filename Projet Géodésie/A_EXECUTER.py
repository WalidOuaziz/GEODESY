import sys
from PyQt5 import QtWidgets, uic
from PyQt5.QtWidgets import *
from fct_0 import *
from fct_UI_Latitudes import *
from ui_fct_langueurArc import *
from ui_fct_surface import *
from fct_Ui_pb_sphere import *
from fct_Ui_Ellipsoide import *
import icone

# fct_0 ----> GUI.ui  
# fct_UI_1 -----> UI_1
# fct_UI_2 -----> UI_2
class F1(QtWidgets.QMainWindow):
    def __init__(self ):
        super(F1 , self).__init__()
        uic.loadUi('Projet Géodésie/UI_project/UI_1.ui' , self)
        self.show()
        self.setWindowTitle("Géodésie")
        F1.setFixedSize(self , 1280, 800) #fixer les dimentions

        # Entre to 'parametre ellipsoide':
        self.Ellipsoide_paramters.clicked.connect(self.switch_fct_0)
        self.ui2 = None
        # Entre to 'calcul des latitdes':
        self.Angle.clicked.connect(self.switch_fct_1)
        self.ui3 = None

        # Entre to 'calcul des latitdes':
        self.Longueur_Arc_2.clicked.connect(self.switch_Longeur)
        self.ui4= None

        # Entre to 'calcul des surfaces':
        self.Surface.clicked.connect(self.switch_surface)
        self.ui5= None

        # Entre to 'pb directe sue sphère':
        self.pb_sphere.clicked.connect(self.switch_pb_sphere)
        self.ui6= None

        # Entre to 'pb directe sue sphère':
        self.pb_elliipsoide.clicked.connect(self.switch_pb_ellipsoide)
        self.ui7= None


    # Entre to 'parametre ellipsoide':
    def switch_fct_0(self):
        if self.ui2 is None :
            self.ui2 = Calculs_Geodisiques(self)
        self.ui2.show() # switch to UI2
        # self.hide() #permet de fermer 1er UI
    
    # Entre to 'calcul des angles':
    def switch_fct_1(self):
        if self.ui3 is None :
            self.ui3 = F2(self)
        self.ui3.show()
    # Entre to 'calcul des latitdes':
    def switch_Longeur(self):
        if self.ui4 is None :
            self.ui4 = F3(self)
        self.ui4.show() # switch to UI3
        # self.hide() #permet de fermer 1er UI

    # Entre to 'calcul des surfaces':
    def switch_surface(self):
        if self.ui5 is None :
            self.ui5 = F4(self)
        self.ui5.show()

    # Entre to 'pb directe sue sphère':
    def switch_pb_sphere(self):
        if self.ui6 is None :
            self.ui6 = F5(self)
        self.ui6.show()

        
    # Entre to 'pb directe sue ellipsoide':
    def switch_pb_ellipsoide(self):
        if self.ui7 is None :
            self.ui7 = F6(self)
        self.ui7.show()
        
if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window  = F1()
    sys.exit(app.exec_())