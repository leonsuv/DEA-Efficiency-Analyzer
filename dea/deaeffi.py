import numpy as np
from dealib.dea.core import dea
from dealib.dea.utils.options import RTS, Orientation
import scipy.optimize as opt
import warnings


class DEAEfficiency:
    """
    Data Envelopment Analysis (DEA) Effizienzanalyse Klasse
    
    Diese Klasse enthält Methoden für verschiedene DEA-Analysen wie:
    - CCR (konstante Skaleneffekte)
    - BCC (variable Skaleneffekte)
    - Skaleneffizienz
    - Kreuzeffizienz
    
    """
    
    def __init__(self, x, y, orientation=Orientation.input, instant_calculation=False):
        """
        Initialisierung der DEA Effizienzanalyse
        
        Parameter:
            x (2D-Array): Input-Matrix, wobei Zeilen die DMUs und Spalten die Inputs sind
            y (2D-Array): Output-Matrix, wobei Zeilen die DMUs und Spalten die Outputs sind
            orientation: Input- oder Output-Orientierung (Standard: Input)
        
        """
        assert len(x) == len(y), "Input und Output müssen gleich viele DMUs haben"
        assert len(x) >= 6, "Mindestens 6 DMUs in x benötigt"
        assert len(x[0]) > 0, "Mindestens eine Input-Variable nötig"
        assert len(y) >= 6, "Mindestens 6 DMUs in y benötigt"
        assert len(y[0]) > 0, "Mindestens eine Output-Variable nötig"
        #assert (len(x[0]) + len(y[0])) * 3 >= len(x), "Mindestens (u+y)*3 DMUs benötigt"
        
        std_len_x = len(x[0])
        if not all(len(xs) == std_len_x for xs in x):
            raise ValueError('Nicht alle Listen haben die gleiche Länge!')
        
        std_len_y = len(y[0])
        if not all(len(ys) == std_len_y for ys in y):
            raise ValueError('Nicht alle Listen haben die gleiche Länge!')

        self.x = np.array(x)
        self.y = np.array(y)
        self.orientation = orientation
        
        self.n_dmus = len(self.x)
        self.m = self.x.shape[1]  # Anzahl der Inputs
        self.s = self.y.shape[1]  # Anzahl der Outputs
        
        # Speicher für Ergebnisse
        self._ccr_result = None
        self._bcc_result = None
        self._scale_efficiency = None
        self._cross_eff_matrix = None
        # Direkte berechnung der werte bei der Initialisierung eines Objekts
        if instant_calculation:
            self.ccr_result = self.get_ccr()
            self.bcc_result = self.get_bcc()
            self.scale_efficiency = self.get_scale_efficiency()
            self.cross_eff_matrix = self.calculate_cross_efficiency()
    
    def get_ccr(self):
        """
        Berechnet CCR (Constant Returns to Scale) Effizienzwerte
        
        Rückgabe:
            numpy.ndarray: CCR Effizienzwerte für jede DMU
        """
        if self._ccr_result is None:
            self._ccr_result = dea(self.x, self.y, rts=RTS.crs, orientation=self.orientation)
        
        return self._ccr_result.eff
    
    def get_bcc(self):
        """
        Berechnet BCC (Variable Returns to Scale) Effizienzwerte
        
        Rückgabe:
            numpy.ndarray: BCC Effizienzwerte für jede DMU
        """
        if self._bcc_result is None:
            self._bcc_result = dea(self.x, self.y, rts=RTS.vrs, orientation=self.orientation)
        
        return self._bcc_result.eff
    
    def get_scale_efficiency(self):
        """
        Berechnet die Skaleneffizienz (CCR/BCC)
        
        Rückgabe:
            numpy.ndarray: Skaleneffizienzwerte für jede DMU
        """
        if self._scale_efficiency is None:
            ccr_scores = self.get_ccr()
            bcc_scores = self.get_bcc()
            self._scale_efficiency = ccr_scores / bcc_scores
        
        return self._scale_efficiency
    
    
    def calculate_cross_efficiency(self, float_values=True):
        """
        Berechnet die Kreuzeffizienzmatrix
        
        Parameter:
            float_values: Wenn True, werden Fließkommawerte zurückgegeben, sonst gerundete Ganzzahlen (als Prozent)
            
        Rückgabe:
            2D-Array: Kreuzeffizienz-Matrix
        
        Hinweis: Der Parameter 'rts' wird in der Doku erwähnt aber nicht verwendet!
        """
        # Standard DEA Effizienzwerte berechnen
        eff_result = dea(self.x, self.y, rts=RTS.crs, orientation=self.orientation)
        ux = eff_result.ux  # Input-Gewichtungen
        vy = eff_result.vy  # Output-Gewichtungen
        
        # Kreuzeffizienz-Matrix initialisieren
        cross_eff_matrix = np.zeros((self.n_dmus, self.n_dmus))
        
        # Für jede DMU optimale Gewichtungen finden und anwenden
        for i in range(self.n_dmus):
            u = ux[i]
            v = vy[i]
            for j in range(self.n_dmus):
                x = 0
                y = 0

                for n in range(len(u)):
                    x += u[n] * self.x[j][n]
                
                for m in range(len(v)):
                    y += v[m] * self.y[j][m]
                
                
                if (float_values):
                    cross_eff_matrix[i][j] = y/x
                else:
                    cross_eff_matrix[i][j] = int(round((y/x)*100))

        return cross_eff_matrix


    def calculate_column_laplace(cross_mat):
        """
        Berechnet die Lapplace Werte
        
        Parameter:
            cross_mat von calculate_cross_efficiency
        Rückgabe:
            Array mit Laplace Werten pro Spalte
        """
        temp_laplace = 0
        laplace_result = np.zeros(len(cross_mat))
        for j in range(len(cross_mat[0])):
            for i in range(len(cross_mat)):
                temp_laplace += cross_mat[i][j]
            laplace_result[j] = temp_laplace/len(cross_mat)
        return laplace_result

    def calculate_column_maxmin(cross_mat):
        """
        Berechnet die Maxmin Werte
        
        Parameter:
            cross_mat von calculate_cross_efficiency
        Rückgabe:
            Array mit Maxmin Werten pro Spalte
        """
        temp_min = 0
        maxmin_result = np.zeros(len(cross_mat))
        for j in range(len(cross_mat[0])):
            for i in range(len(cross_mat)):
                if cross_mat[i][j] < temp_min:
                    temp_min = cross_mat[i][j]
            maxmin_result[j] = temp_min
        return maxmin_result

    def calculate_column_maxmax(cross_mat):
        """
        Berechnet die Maxmax Werte
        
        Parameter:
            cross_mat von calculate_cross_efficiency
        Rückgabe:
            Array mit Maxmax Werten pro Spalte
        """
        temp_max = 0
        maxmax_result = np.zeros(len(cross_mat))
        for j in range(len(cross_mat[0])):
            for i in range(len(cross_mat)):
                if cross_mat[i][j] > temp_max:
                    temp_max = cross_mat[i][j]
            maxmax_result[j] = temp_max
        return maxmax_result


