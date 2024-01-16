import sys
from PyQt5 import QtWidgets, uic
from PyQt5.QtWidgets import *
from New_function import *
from fonction import *
import math
class F2(QtWidgets.QMainWindow):
    def __init__(self , parent = None):
        super(F2 , self).__init__(parent)
        uic.loadUi('Projet Géodésie/UI_project/UI_2.ui' , self)
        self.show()
        self.setWindowTitle("Latitudes")
        F2.setFixedSize(self , 615 , 310) #fixer les dimentions


        self.calculer.clicked.connect(self.calcul)
        self.clear.clicked.connect(self.remove)
    def remove(self):
            self.grandAxeA.clear()
            self.petitAxeB.clear()
            self.phi_deg.clear()
            self.phi_min.clear()
            self.phi_sec.clear()
            self.beta_deg.clear()
            self.beta_min.clear()
            self.beta_sec.clear()
            self.psy_deg.clear()
            self.psy_min.clear()
            self.psy_sec.clear()

    def calcul(self):
        # a = 6378249.145
        # b = 6356514.9658
        try:
            a = float(str(eval(self.grandAxeA.text())))
            b = float(str(eval(self.petitAxeB.text())))
            try :
                phi_degValue = float(str(eval(self.phi_deg.text())))
                phi_minValue = float(str(eval(self.phi_min.text())))
                phi_secValue = float(str(eval(self.phi_sec.text())))
                phi_Value = (phi_degValue + phi_minValue/60 + phi_secValue/3600) * math.pi/180
                
                beta_Value = beta_from_phi(a , b , phi_Value)
                psy_Value = psi_from_phi(a , b , phi_Value)

                beta_degValue = dms(beta_Value)[0]
                beta_minValue = dms(beta_Value)[1]
                beta_secValue = dms(beta_Value)[2]
                self.beta_deg.setText(str(beta_degValue))
                self.beta_min.setText(str(beta_minValue))
                self.beta_sec.setText(str(beta_secValue))

                psy_degValue = dms(psy_Value)[0]
                psy_minValue = dms(psy_Value)[1]
                psy_secValue = dms(psy_Value)[2]
                self.psy_deg.setText(str(psy_degValue))
                self.psy_min.setText(str(psy_minValue))
                self.psy_sec.setText(str(psy_secValue))


                tracer_latitude(phi_Value*180/pi, psy_Value *180/pi,beta_Value*180/pi)

                # self.graph.canvas(visAngle(phi_Value *180/pi, beta_Value *180/pi, psy_Value*180/pi))

                # VisAng(phi_Value *180/pi, beta_Value *180/pi, psy_Value*180/pi)
                # fig, ax = visAngle(phi_Value *180/pi, beta_Value *180/pi, psy_Value*180/pi)
                
                # canvas = FigureCanvas(fig)
                # layout = QVBoxLayout(self.graph)
                # layout.addWidget(canvas)
                # self.graph.setLayout(layout)
                
                
            
            except : 
                try:
                    beta_degValue = float(str(eval(self.beta_deg.text())))
                    beta_minValue = float(str(eval(self.beta_min.text())))
                    beta_secValue = float(str(eval(self.beta_sec.text())))
                    beta_Value = (beta_degValue + beta_minValue/60 + beta_secValue/3600) * math.pi/180

                    phi_Value = phi_from_beta(a , b , beta_Value)
                    psy_Value = psi_from_beta(a , b , beta_Value)

                    phi_degValue = dms(phi_Value)[0]
                    phi_minValue = dms(phi_Value)[1]
                    phi_secValue = dms(phi_Value)[2]
                    self.phi_deg.setText(str(phi_degValue))
                    self.phi_min.setText(str(phi_minValue))
                    self.phi_sec.setText(str(phi_secValue))

                    psy_degValue = dms(psy_Value)[0]
                    psy_minValue = dms(psy_Value)[1]
                    psy_secValue = dms(psy_Value)[2]
                    self.psy_deg.setText(str(psy_degValue))
                    self.psy_min.setText(str(psy_minValue))
                    self.psy_sec.setText(str(psy_secValue))

                    tracer_latitude(phi_Value*180/pi, psy_Value *180/pi,beta_Value*180/pi)

                    
                    
                except:
                    try:
                        psy_degValue = float(str(eval(self.psy_deg.text())))
                        psy_minValue = float(str(eval(self.psy_min.text())))
                        psy_secValue = float(str(eval(self.psy_sec.text())))
                        psy_Value = (psy_degValue + psy_minValue/60 + psy_secValue/3600) * math.pi/180
                        
                        phi_Value = phi_from_psi(a , b , psy_Value)
                        beta_Value = beta_from_psi(a , b , psy_Value)
                        

                        beta_degValue = dms(beta_Value)[0]
                        beta_minValue = dms(beta_Value)[1]
                        beta_secValue = dms(beta_Value)[2]
                        self.beta_deg.setText(str(beta_degValue))
                        self.beta_min.setText(str(beta_minValue))
                        self.beta_sec.setText(str(beta_secValue))

                        phi_degValue = dms(phi_Value)[0]
                        phi_minValue = dms(phi_Value)[1]
                        phi_secValue = dms(phi_Value)[2]
                        self.phi_deg.setText(str(phi_degValue))
                        self.phi_min.setText(str(phi_minValue))
                        self.phi_sec.setText(str(phi_secValue))
                        tracer_latitude(phi_Value*180/pi, psy_Value *180/pi,beta_Value*180/pi)

                        
                    except :
                        print('error')
                        from PyQt5.QtWidgets import QMessageBox
                        msgBox = QMessageBox()
                        msgBox.setIcon(QMessageBox.Critical)
                        msgBox.setWindowTitle("Erreur")
                        msgBox.setText("Vérifier les paramètres saisis !")
                        msgBox.exec_()

                    
        except:
            print('error')
            from PyQt5.QtWidgets import QMessageBox
            msgBox = QMessageBox()
            msgBox.setIcon(QMessageBox.Critical)
            msgBox.setWindowTitle("Erreur")
            msgBox.setText("Vérifier les paramètres saisis !")
            msgBox.exec_()


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window  = F2()
    sys.exit(app.exec_())