import sys
from PyQt5 import QtWidgets, uic
from PyQt5.QtWidgets import *
from fonction import *
import icone

class Calculs_Geodisiques(QtWidgets.QMainWindow):
    def __init__(self , parent = None ):
        super(Calculs_Geodisiques , self).__init__(parent)
        uic.loadUi('Projet Géodésie\GUI.ui' , self)
        self.show()
        self.setWindowTitle("Calculs Géodésiques")
        Calculs_Geodisiques.setFixedSize(self , 741, 381)

        self.pushButton.clicked.connect(self.inPut)
        self.pushButton_2.clicked.connect(self.Clear)
        self.pushButton_6.clicked.connect(self.Clear1)
        self.pushButton_5.clicked.connect(self.Clear2)
        self.pushButton_32.clicked.connect(self.Clear3)
        self.pushButton_14.clicked.connect(self.Clear4)
        
        self.pushButton_3.clicked.connect(self.convertir)
        self.pushButton_4.clicked.connect(self.convertir1)
        self.pushButton_31.clicked.connect(self.convertir2)
        self.pushButton_13.clicked.connect(self.convertir3)
        

        self.comboBox.currentIndexChanged.connect(self.combox_select) # permetd de creer une fonction convert select qui stocke la valeur choisie par utilisateur
        self.comboBox_2.currentIndexChanged.connect(self.combox_select1)
        
    # recuperer la valeur choisie au ComBox
    def combox_select(self):
        Zone = self.comboBox.currentText()
        return int(Zone)

    def combox_select1(self):
        Zone = self.comboBox_2.currentText()
        return int(Zone)

    def Clear4(self):
        self.lineEdit_22.clear()
        self.lineEdit_23.clear()
        self.lineEdit_24.clear()
        self.lineEdit_16.clear()
        self.lineEdit_18.clear()
        self.lineEdit_17.clear()

    def convertir3(self):
        try :
            A_rad = float(str(eval(self.lineEdit_22.text())))
            # A_deg_dec = float(str(eval(self.lineEdit_23.text())))
            # A_gon = float(str(eval(self.lineEdit_24.text())))

            # A_dms_deg = float(str(eval(self.lineEdit_16.text())))
            # A_dms_min = float(str(eval(self.lineEdit_18.text())))
            # A_dms_sec = float(str(eval(self.lineEdit_17.text())))

            # if str(A_rad)!="":
            A_deg_dec = angle_rad(A_rad)[0]
            A_gon = angle_rad(A_rad)[1]
            A_dms_deg = angle_rad(A_rad)[2][0]
            A_dms_min = angle_rad(A_rad)[2][1]
            A_dms_sec = angle_rad(A_rad)[2][2]
            self.lineEdit_23.setText(str(A_deg_dec))
            self.lineEdit_24.setText(str(A_gon))
            self.lineEdit_16.setText(str(A_dms_deg))
            self.lineEdit_18.setText(str(A_dms_min))
            self.lineEdit_17.setText(str(A_dms_sec))
        except :
            try :
                A_deg_dec = float(str(eval(self.lineEdit_23.text())))
                A_rad = angle_deg_dec(A_deg_dec)[0]
                A_gon = angle_deg_dec(A_deg_dec)[1]
                A_dms_deg = angle_deg_dec(A_deg_dec)[2][0]
                A_dms_min = angle_deg_dec(A_deg_dec)[2][1]
                A_dms_sec = angle_deg_dec(A_deg_dec)[2][2]
                self.lineEdit_22.setText(str(A_rad))
                self.lineEdit_24.setText(str(A_gon))
                self.lineEdit_16.setText(str(A_dms_deg))
                self.lineEdit_18.setText(str(A_dms_min))
                self.lineEdit_17.setText(str(A_dms_sec))
            except :
                try :
                    A_gon = float(str(eval(self.lineEdit_24.text())))
                    A_rad = angle_gon(A_gon)[0]
                    A_deg_dec = angle_gon(A_gon)[1]
                    A_dms_deg = angle_gon(A_gon)[2][0]
                    A_dms_min = angle_gon(A_gon)[2][1]
                    A_dms_sec = angle_gon(A_gon)[2][2]
                    self.lineEdit_22.setText(str(A_rad))
                    self.lineEdit_23.setText(str(A_deg_dec))
                    self.lineEdit_16.setText(str(A_dms_deg))
                    self.lineEdit_18.setText(str(A_dms_min))
                    self.lineEdit_17.setText(str(A_dms_sec))
                except:
                    try :
                        A_dms_deg = float(str(eval(self.lineEdit_16.text())))
                        A_dms_min = float(str(eval(self.lineEdit_18.text())))
                        A_dms_sec= float(str(eval(self.lineEdit_17.text())))

                        A_rad = angle_dms(A_dms_deg, A_dms_min, A_dms_sec)[0]
                        A_deg_dec = angle_dms(A_dms_deg, A_dms_min, A_dms_sec)[1]
                        A_gon = angle_dms(A_dms_deg, A_dms_min, A_dms_sec)[2]

                        self.lineEdit_22.setText(str(A_rad))
                        self.lineEdit_23.setText(str(A_deg_dec))
                        self.lineEdit_24.setText(str(A_gon))
                    except :
                        from PyQt5.QtWidgets import QMessageBox

                        # create a message box with a critical icon
                        msgBox = QMessageBox()
                        msgBox.setIcon(QMessageBox.Critical)

                        # set the window title and message text
                        msgBox.setWindowTitle("Erreur")
                        msgBox.setText("Vérifier les paramètres saisis !")

                        # display the message box
                        msgBox.exec_()

    def convertir2(self):
        try :
            a = float(str(eval(self.lineEdit_27.text())))
            b = float(str(eval(self.lineEdit_28.text())))
            phi_0_deg = float(str(eval(self.lineEdit_29.text())))
            phi_0_min = float(str(eval(self.lineEdit_30.text())))/60
            phi_0_sec = float(str(eval(self.lineEdit_31.text())))/3600

            phi_0 = phi_0_deg + phi_0_min + phi_0_sec

            U = float(str(eval(self.lineEdit_68.text())))
            P = float(str(eval(self.lineEdit_88.text())))

            phi = calculer_phi(a, b, U, phi_0 , P)
            phi_deg = dms(phi)[0]
            phi_min = dms(phi)[1]
            phi_sec = dms(phi)[2]

            self.label_268.setText(str(phi_deg) + "°")
            self.label_269.setText(str(phi_min) + "'")
            self.label_266.setText(str(phi_sec) + "\"")
        except :

            from PyQt5.QtWidgets import QMessageBox

            # create a message box with a critical icon
            msgBox = QMessageBox()
            msgBox.setIcon(QMessageBox.Critical)

            # set the window title and message text
            msgBox.setWindowTitle("Erreur")
            msgBox.setText("Une erreur s'est produite.")

            # display the message box
            msgBox.exec_()

    def Clear3(self):
        self.lineEdit_27.clear()
        self.lineEdit_28.clear()
        self.lineEdit_29.clear()
        self.lineEdit_30.clear()
        self.lineEdit_31.clear()
        self.lineEdit_68.clear()
        self.lineEdit_88.clear()
        self.label_268.clear()
        self.label_269.clear()
        self.label_266.clear()

    def convertir1(self):
        try :
            def combox_select1(self):
                Zone1 = self.comboBox_2.currentText()
                return int(Zone1)

            X = float(str(eval(self.lineEdit_12.text())))
            Y = float(str(eval(self.lineEdit_13.text())))
            Zone1 = combox_select1(self)
            phi = coord_Lambert_vers_Geo(X, Y, Zone1)[0]
            L = coord_Lambert_vers_Geo(X, Y, Zone1)[1]
            phi_deg = dms(phi)[0]
            phi_min = dms(phi)[1]
            phi_sec = dms(phi)[2]

            L_deg = abs(dms(L)[0])
            L_min = abs(dms(L)[1])
            L_sec = abs(dms(L)[2])
            self.label_45.setText(str(phi_deg) + "°")
            self.label_48.setText(str(phi_min) + "'")
            self.label_49.setText(str(phi_sec) + "\"")

            self.label_50.setText(str(L_deg) + "°")
            self.label_51.setText(str(L_min) + "'")
            self.label_52.setText(str(L_sec) + "\"")
        except :
            from PyQt5.QtWidgets import QMessageBox

            # create a message box with a critical icon
            msgBox = QMessageBox()
            msgBox.setIcon(QMessageBox.Critical)

            # set the window title and message text
            msgBox.setWindowTitle("Erreur")
            msgBox.setText("Une erreur s'est produite.")

            # display the message box
            msgBox.exec_()

    def convertir(self):
        try :
            def combox_select(self):
                Zone = self.comboBox.currentText()
                return int(Zone)

            phi_deg = float(str(eval(self.lineEdit_6.text())))
            phi_min = float(str(eval(self.lineEdit_8.text())))/60
            phi_sec = float(str(eval(self.lineEdit_9.text())))/3600
            phi = phi_deg + phi_min +phi_sec
          
      
            L_deg = -float(str(eval(self.lineEdit_7.text())))
            L_min = -float(str(eval(self.lineEdit_11.text())))/60
            L_sec = -float(str(eval(self.lineEdit_10.text())))/3600
            L = L_deg + L_min + L_sec
        
            Zone = combox_select(self)
            X = geo_vers_Lambert(phi, L, Zone)[0]
            Y = geo_vers_Lambert(phi, L, Zone)[1]
            self.label_5.setText(str(X))
            self.label_18.setText(str(Y))
        except :
            from PyQt5.QtWidgets import QMessageBox

            # create a message box with a critical icon
            msgBox = QMessageBox()
            msgBox.setIcon(QMessageBox.Critical)

            # set the window title and message text
            msgBox.setWindowTitle("Erreur")
            msgBox.setText("Une erreur s'est produite.")

            # display the message box
            msgBox.exec_()
      
    def Clear(self):
        self.lineEdit_1.clear()
        self.lineEdit_2.clear()
        self.lineEdit_3.clear()
        self.lineEdit_4.clear()
        self.lineEdit_5.clear()
        self.lineEdit_14.clear()
        self.lineEdit_15.clear()
        self.label_2.clear()
        self.label_4.clear()
        self.label_28.clear()
    def Clear2(self):
        self.label_5.clear()
        self.label_18.clear()
        self.lineEdit_6.clear()
        self.lineEdit_7.clear()
        self.lineEdit_8.clear()
        self.lineEdit_9.clear()
        self.lineEdit_10.clear()
        self.lineEdit_11.clear()
        
    def Clear1(self):
        self.label_45.clear()
        self.label_48.clear()
        self.label_49.clear()
        self.label_50.clear()
        self.label_51.clear()
        self.label_52.clear()
        self.lineEdit_12.clear()
        self.lineEdit_13.clear()
       
    def inPut(self):
        try :
            a = float(str(eval(self.lineEdit_1.text()))) #a
            b = float(str(eval(self.lineEdit_2.text()))) #b
            f = aplatissement(a, b)
            e2 = excentricite(a, b)
            self.lineEdit_3.setText(str(f))
            self.lineEdit_4.setText(str(e2))
            if self.lineEdit_5.text()!="" :

                phi_deg = float(str(eval(self.lineEdit_5.text())))
                phi_min = float(str(eval(self.lineEdit_14.text())))/60
                phi_sec = float(str(eval(self.lineEdit_15.text())))/3600
                phi = phi_deg + phi_min +phi_sec
                phi = phi *pi/180

                M = Rayon_courbure(a, b, phi)[1]
                N = Rayon_courbure(a, b, phi)[0]
                u = U(a, b, phi)
                self.label_28.setText(str(u))
                self.label_2.setText(str(M))
                self.label_4.setText(str(N))
        except :
            try :
                a = float(str(eval(self.lineEdit_1.text()))) #a
                f = float(str(eval(self.lineEdit_3.text()))) #f
                b = petit_axe1(a, f)
                e2 = excentricite(a, b)
                self.lineEdit_2.setText(str(b))
                self.lineEdit_4.setText(str(e2))
                if self.lineEdit_5.text()!="":
                    phi_deg = float(str(eval(self.lineEdit_5.text())))
                    phi_min = float(str(eval(self.lineEdit_14.text())))/60
                    phi_sec = float(str(eval(self.lineEdit_15.text())))/3600
                    phi = phi_deg + phi_min +phi_sec
                    phi = phi *pi/180
                    M = Rayon_courbure(a, b, phi)[1]
                    N = Rayon_courbure(a, b, phi)[0]
                    u = U(a, b, phi)
                    self.label_28.setText(str(u))
                    self.label_2.setText(str(M))
                    self.label_4.setText(str(N))
                
            except :
                try :
                    a = float(str(eval(self.lineEdit_1.text()))) #a
                    e2 = float(str(eval(self.lineEdit_4.text()))) #e2
                    b = petit_axe2(a, e2)
                    f = aplatissement(a, b)
                    self.lineEdit_2.setText(str(b))
                    self.lineEdit_3.setText(str(f))
                    if self.lineEdit_5.text()!="":
                        phi_deg = float(str(eval(self.lineEdit_5.text())))
                        phi_min = float(str(eval(self.lineEdit_14.text())))/60
                        phi_sec = float(str(eval(self.lineEdit_15.text())))/3600
                        phi = phi_deg + phi_min +phi_sec
                        phi = phi *pi/180
                        M = Rayon_courbure(a, b, phi)[1]
                        N = Rayon_courbure(a, b, phi)[0]
                        u = U(a, b, phi)
                        self.label_28.setText(str(u))
                        self.label_2.setText(str(M))
                        self.label_4.setText(str(N))
                except :
                    try :
                        b = float(str(eval(self.lineEdit_2.text()))) 
                        e2 = float(str(eval(self.lineEdit_4.text()))) 
                        a = grand_axe2(b, e2)
                        f = aplatissement(a, b)
                        self.lineEdit_1.setText(str(a))
                        self.lineEdit_3.setText(str(f))
                        if self.lineEdit_5.text()!="":
                            phi_deg = float(str(eval(self.lineEdit_5.text())))
                            phi_min = float(str(eval(self.lineEdit_14.text())))/60
                            phi_sec = float(str(eval(self.lineEdit_15.text())))/3600
                            phi = phi_deg + phi_min +phi_sec
                            phi = phi *pi/180
                            M = Rayon_courbure(a, b, phi)[1]
                            N = Rayon_courbure(a, b, phi)[0]
                            u = U(a, b, phi)
                            self.label_28.setText(str(u))
                            self.label_2.setText(str(M))
                            self.label_4.setText(str(N))
                    except :
                        try :
                            b = float(str(eval(self.lineEdit_2.text()))) 
                            f = float(str(eval(self.lineEdit_3.text()))) 
                            a = grand_axe1(b, f)
                            e2 = excentricite(a, b)
                            self.lineEdit_1.setText(str(a))
                            self.lineEdit_4.setText(str(e2))
                            if self.lineEdit_5.text()!="":
                                phi_deg = float(str(eval(self.lineEdit_5.text())))
                                phi_min = float(str(eval(self.lineEdit_14.text())))/60
                                phi_sec = float(str(eval(self.lineEdit_15.text())))/3600
                                phi = phi_deg + phi_min +phi_sec
                                phi = phi *pi/180
                                M = Rayon_courbure(a, b, phi)[1]
                                N = Rayon_courbure(a, b, phi)[0]
                                u = U(a, b, phi)
                                self.label_28.setText(str(u))
                                self.label_2.setText(str(M))
                                self.label_4.setText(str(N))
                        except :
                            from PyQt5.QtWidgets import QMessageBox

                            # create a message box with a critical icon
                            msgBox = QMessageBox()
                            msgBox.setIcon(QMessageBox.Critical)

                            # set the window title and message text
                            msgBox.setWindowTitle("Erreur")
                            msgBox.setText("Une erreur s'est produite.")

                            # display the message box
                            msgBox.exec_()

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window  = Calculs_Geodisiques()
    sys.exit(app.exec_())