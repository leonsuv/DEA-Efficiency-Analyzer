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
    
    def __init__(self, x, y, orientation=Orientation.input, instant_calculation=False, const_or_variable=RTS.crs):
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
        
        assert const_or_variable in (RTS.crs, RTS.vrs), "Nur CRS (Konstant) oder VRS (Variabel) auswählbar"

        self.x = np.array(x)
        self.y = np.array(y)
        self.orientation = orientation
        
        self.n_dmus = len(self.x)
        self.m = self.x.shape[1]  # Anzahl der Inputs
        self.s = self.y.shape[1]  # Anzahl der Outputs

        self.rts = const_or_variable
        
        # Speicher für Ergebnisse
        self._ccr_result = None
        self._bcc_result = None
        self._scale_efficiency = None
        self._cross_eff_matrix = None
        self._laplace_values = None
        self._maxmin_values = None
        self._maxmax_values = None
        self._maverick_values = None
        
        # Direkte Berechnung bei Initialisierung
        if instant_calculation:
            self.ccr_result = self.get_ccr()
            self.bcc_result = self.get_bcc()
            self.scale_efficiency = self.get_scale_efficiency()
            self.cross_eff_matrix = self.calculate_cross_efficiency()
            self.laplace_values = self.calculate_column_laplace(self.cross_eff_matrix)
            self.maxmin_values = self.calculate_column_maxmin(self.cross_eff_matrix)
            self.maxmax_values = self.calculate_column_maxmax(self.cross_eff_matrix)
            self.maverick_values = self.calculate_column_maverick(self.cross_eff_matrix)
    
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
        Berechnet die Kreuzeffizienzmatrix mit Berücksichtigung der Orientierung.
        
        Falls self.rts == RTS.crs wird CCR 2nd Phase verwendet:
        effₖₗ* = (uₖᵀ yₗ)/(vₖᵀ xₗ)
        
        Falls self.rts == RTS.vrs wird BCC 2nd Phase verwendet:
        effₖₗ* = (Uₖᵀ yₗ + uₖ)/(Vₖᵀ xₗ)
        
        Die Orientierung bestimmt die Normalisierung im LP:
        - Input-Orientierung: vᵀ xₖ = 1, max uᵀ yₖ
        - Output-Orientierung: uᵀ yₖ = 1, min vᵀ xₖ
        """
        cross_eff_matrix = np.zeros((self.n_dmus, self.n_dmus))
        
        if self.rts == RTS.crs:
            # Second-phase CCR using Model 6
            for i in range(self.n_dmus):
                u, v = self._solve_second_phase_ccr(i)
                print(f"{i}: {u}, {v}")
                if u is None or v is None:
                    warnings.warn(f"LP für DMU {i} konvergierte nicht; überspringe DMU i.")
                    continue
                for j in range(self.n_dmus):
                    denom = np.dot(v, self.x[j])
                    # Avoid division by zero
                    if denom == 0:
                        eff = 0
                    else:
                        eff = np.dot(u, self.y[j]) / denom
                    if float_values:
                        cross_eff_matrix[i][j] = eff
                    else:
                        cross_eff_matrix[i][j] = int(round(eff * 100))
        
        else:
            # Second-phase BCC using Model 5
            for i in range(self.n_dmus):
                U, V, u_free = self._solve_second_phase_bcc(i)
                if U is None or V is None:
                    warnings.warn(f"BCC LP für DMU {i} konvergierte nicht; überspringe DMU i.")
                    continue
                for j in range(self.n_dmus):
                    denom = np.dot(V, self.x[j])
                    if denom == 0:
                        eff = 0
                    else:
                        eff = (np.dot(U, self.y[j]) + u_free) / denom
                    if float_values:
                        cross_eff_matrix[i][j] = eff
                    else:
                        cross_eff_matrix[i][j] = int(round(eff * 100))
                        
        return cross_eff_matrix


    def _solve_second_phase_ccr(self, dmu_index):
        """
        LP für CCR (CRS) Second Phase mit Orientierung gemäß den Formeln:
        
        Output-Orientierung:
        min  sum(s_j)
        s.t. U_k^T y_k = 1
            -U_k^T y_j + V_k^T x_j ≥ 0,  ∀ j
            V_k^T x_k = (eff_k*)^(-1)
            -U_k^T y_j + V_k^T x_j - s_j ≤ 0,  ∀ j
            U_k, V_k ≥ 0, s_j ≥ 0
            
        Input-Orientierung:
        min  sum(s_j)
        s.t. V_k^T x_k = 1
            -V_k^T x_j + U_k^T y_j ≤ 0,  ∀ j
            U_k^T y_k = eff_k*
            -V_k^T x_j + U_k^T y_j + s_j ≥ 0,  ∀ j
            V_k, U_k ≥ 0, s_j ≥ 0
        """
        s_dim = self.s  # Anzahl Output-Variablen
        m_dim = self.m  # Anzahl Input-Variablen
        n = self.n_dmus
        
        # First, get the efficiency score for the DMU being evaluated
        if self.orientation == Orientation.output:
            eff_result = dea(self.x, self.y, rts=RTS.crs, orientation=Orientation.output)
            eff_k = eff_result.eff[dmu_index]
            
            # Decision variables: [U (s_dim), V (m_dim), s (n_dmus)]
            num_variables = s_dim + m_dim + 1 + n
            
            # Objective: minimize sum of slacks
            c = [0] * (s_dim + m_dim + 1) + [1] * n
            
            # Constraint 1: U_k^T y_k = 1
            A_eq = [[0] * num_variables]  # Will fill below
            for i in range(s_dim):
                A_eq[0][i] = self.y[dmu_index][i]
            b_eq = [1]
            
            # Constraint 2: V_k^T x_k = (eff_k*)^(-1)
            A_eq.append([0] * num_variables)
            for i in range(m_dim):
                A_eq[1][s_dim + i] = self.x[dmu_index][i]
            b_eq.append(1.0 / eff_k)
            
            # Constraints 3 & 4: For each DMU
            A_ub = []
            b_ub = []
            
            # Constraint 3: -U_k^T y_j + V_k^T x_j ≥ 0
            for j in range(n):
                row = [0] * num_variables
                # -U_k^T y_j
                for i in range(s_dim):
                    row[i] = -self.y[j][i]
                # +V_k^T x_j
                for i in range(m_dim):
                    row[s_dim + i] = self.x[j][i]
                # +u_k
                row[s_dim + m_dim] = 1
                A_ub.append([-val for val in row])  # Flip sign for ≥ to ≤
                b_ub.append(0)
            
            # Constraint 4: -U_k^T y_j + V_k^T x_j - s_j ≤ 0
            for j in range(n):
                row = [0] * num_variables
                # -U_k^T y_j
                for i in range(s_dim):
                    row[i] = -self.y[j][i]
                # +V_k^T x_j
                for i in range(m_dim):
                    row[s_dim + i] = self.x[j][i]
                # +u_k
                row[s_dim + m_dim] = 1
                # -s_j
                row[s_dim + m_dim + j] = -1
                A_ub.append(row)
                b_ub.append(0)
            
        else:  # Input orientation
            # For input orientation, we need to calculate CCR score
            eff_result = dea(self.x, self.y, rts=RTS.crs, orientation=Orientation.input)
            eff_k = eff_result.eff[dmu_index]
            
            # Decision variables: [U (s_dim), V (m_dim), s (n_dmus)]
            num_variables = s_dim + m_dim + 1 + n
            
            # Objective: minimize sum of slacks
            c = [0] * (s_dim + m_dim + 1) + [1] * n
            
            # Constraint 1: V_k^T x_k = 1
            A_eq = [[0] * num_variables]  # Will fill below
            for i in range(m_dim):
                A_eq[0][s_dim + i] = self.x[dmu_index][i]
            b_eq = [1]
            
            # Constraint 2: U_k^T y_k = eff_k*
            A_eq.append([0] * num_variables)
            for i in range(s_dim):
                A_eq[1][i] = self.y[dmu_index][i]
            A_eq[1][s_dim + m_dim] = 1  # u_k coefficient
            b_eq.append(eff_k)
            
            # Constraints 3 & 4: For each DMU
            A_ub = []
            b_ub = []
            
            # Constraint 3: -V_k^T x_j + U_k^T y_j ≤ 0
            for j in range(n):
                row = [0] * num_variables
                # -V_k^T x_j
                for i in range(m_dim):
                    row[s_dim + i] = -self.x[j][i]
                # +U_k^T y_j
                for i in range(s_dim):
                    row[i] = self.y[j][i]
                # +u_k
                row[s_dim + m_dim] = 1
                A_ub.append(row)
                b_ub.append(0)
            
            # Constraint 4: -V_k^T x_j + U_k^T y_j + s_j ≥ 0
            for j in range(n):
                row = [0] * num_variables
                # -V_k^T x_j
                for i in range(m_dim):
                    row[s_dim + i] = -self.x[j][i]
                # +U_k^T y_j
                for i in range(s_dim):
                    row[i] = self.y[j][i]
                # +u_k
                row[s_dim + m_dim] = 1
                # +s_j
                row[s_dim + m_dim + j] = 1
                A_ub.append([-val for val in row])  # Flip sign for ≥ to ≤
                b_ub.append(0)
        
        # Bounds: U, V, s ≥ 0
        bounds = [(0, None)] * s_dim + [(0, None)] * m_dim + [(0, 0)] + [(0, None)] * n
        
        # Solve the LP
        res = opt.linprog(c=c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                        bounds=bounds, method='highs')
        
        print(f"c={c}\nA_ub={A_ub}\nb_ub={b_ub}\nA_eq={A_eq}\nb_eq={b_eq}\nbounds={bounds}")
        print(f"res: {res}")
        if res.success:
            sol = res.x
            u = sol[:s_dim]
            v = sol[s_dim:s_dim+m_dim]
            return u, v
        else:
            return None, None


    def _solve_second_phase_bcc(self, dmu_index):
        """
        LP für BCC (VRS) Second Phase mit Orientierung.
        
        Input-Orientierung:
        min  sum(s_j)
        s.t. V_k^T x_k = 1
             -V_k^T x_j + U_k^T y_j + u_k ≤ 0  ∀ j
             U_k^T y_k + u_k = eff_k*
             -V_k^T x_j + U_k^T y_j + u_k + s_j ≥ 0  ∀ j
             U_k, V_k ≥ 0, s_j ≥ 0, u_k free
             
        Output-Orientierung:
        min  sum(s_j)
        s.t. U_k^T y_k = 1
             -U_k^T y_j + V_k^T x_j + u_k ≥ 0  ∀ j
             V_k^T x_k = (eff_k*)^(-1)
             -U_k^T y_j + V_k^T x_j + u_k - s_j ≤ 0  ∀ j
             U_k, V_k ≥ 0, s_j ≥ 0, u_k free
        """
        s_dim = self.s
        m_dim = self.m
        n = self.n_dmus
        
        # First, get the efficiency score for the DMU being evaluated
        if self.orientation == Orientation.output:
            eff_result = dea(self.x, self.y, rts=RTS.vrs, orientation=Orientation.output)
            eff_k = eff_result.eff[dmu_index]
            
            # Decision variables: [U (s_dim), V (m_dim), u0 (1), s (n_dmus)]
            num_variables = s_dim + m_dim + 1 + n
            
            # Objective: minimize sum of slacks
            c = [0] * (s_dim + m_dim + 1) + [1] * n
            
            # Constraint 1: U_k^T y_k = 1
            A_eq = [[0] * num_variables]
            for i in range(s_dim):
                A_eq[0][i] = self.y[dmu_index][i]
            b_eq = [1]
            
            # Constraint 2: V_k^T x_k = (eff_k*)^(-1)
            A_eq.append([0] * num_variables)
            for i in range(m_dim):
                A_eq[1][s_dim + i] = self.x[dmu_index][i]
            b_eq.append(1.0 / eff_k)
            
            # Constraints for each DMU
            A_ub = []
            b_ub = []
            
            # Constraint: -U_k^T y_j + V_k^T x_j + u_k ≥ 0
            for j in range(n):
                row = [0] * num_variables
                # -U_k^T y_j
                for i in range(s_dim):
                    row[i] = self.y[j][i]
                # +V_k^T x_j
                for i in range(m_dim):
                    row[s_dim + i] = -self.x[j][i]
                # +u_k
                row[s_dim + m_dim] = -1
                A_ub.append(row)
                b_ub.append(0)
            
            # Constraint: -U_k^T y_j + V_k^T x_j + u_k - s_j ≤ 0
            for j in range(n):
                row = [0] * num_variables
                # -U_k^T y_j
                for i in range(s_dim):
                    row[i] = -self.y[j][i]
                # +V_k^T x_j
                for i in range(m_dim):
                    row[s_dim + i] = self.x[j][i]
                # +u_k
                row[s_dim + m_dim] = 1
                # -s_j
                row[s_dim + m_dim + 1 + j] = -1
                A_ub.append(row)
                b_ub.append(0)
                
        else:  # Input orientation
            eff_result = dea(self.x, self.y, rts=RTS.vrs, orientation=Orientation.input)
            eff_k = eff_result.eff[dmu_index]
            
            # Decision variables: [U (s_dim), V (m_dim), u0 (1), s (n_dmus)]
            num_variables = s_dim + m_dim + 1 + n
            
            # Objective: minimize sum of slacks
            c = [0] * (s_dim + m_dim + 1) + [1] * n
            
            # Constraint 1: V_k^T x_k = 1
            A_eq = [[0] * num_variables]
            for i in range(m_dim):
                A_eq[0][s_dim + i] = self.x[dmu_index][i]
            b_eq = [1]
            
            # Constraint 2: U_k^T y_k + u_k = eff_k*
            A_eq.append([0] * num_variables)
            for i in range(s_dim):
                A_eq[1][i] = self.y[dmu_index][i]
            A_eq[1][s_dim + m_dim] = 1  # u_k coefficient
            b_eq.append(eff_k)
            
            # Constraints for each DMU
            A_ub = []
            b_ub = []
            
            # Constraint: -V_k^T x_j + U_k^T y_j + u_k ≤ 0
            for j in range(n):
                row = [0] * num_variables
                # -V_k^T x_j
                for i in range(m_dim):
                    row[s_dim + i] = -self.x[j][i]
                # +U_k^T y_j
                for i in range(s_dim):
                    row[i] = self.y[j][i]
                # +u_k
                row[s_dim + m_dim] = 1
                A_ub.append(row)
                b_ub.append(0)
            
            # Constraint: -V_k^T x_j + U_k^T y_j + u_k + s_j ≥ 0
            for j in range(n):
                row = [0] * num_variables
                # -V_k^T x_j
                for i in range(m_dim):
                    row[s_dim + i] = self.x[j][i]
                # +U_k^T y_j
                for i in range(s_dim):
                    row[i] = -self.y[j][i]
                # +u_k
                row[s_dim + m_dim] = -1
                # +s_j
                row[s_dim + m_dim + 1 + j] = -1
                A_ub.append(row)
                b_ub.append(0)
        
        # Bounds: U, V, s ≥ 0, u_k free
        bounds = [(0, None)] * s_dim + [(0, None)] * m_dim + [(None, None)] + [(0, None)] * n
        
        # Solve the LP
        res = opt.linprog(c=c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                        bounds=bounds, method='highs')
        if res.success:
            sol = res.x
            u = sol[:s_dim]
            v = sol[s_dim:s_dim+m_dim]
            u0 = sol[s_dim+m_dim]
            return u, v, u0
        else:
            return None, None, None


    def calculate_column_laplace(self, cross_mat):
        """
        Berechnet die Laplace Werte (Durchschnitt je Spalte)
        
        Parameter:
            cross_mat von calculate_cross_efficiency
        Rückgabe:
            Array mit Laplace Werten pro Spalte
        """
        laplace_result = np.zeros(len(cross_mat))
        for j in range(len(cross_mat[0])):
            temp_laplace = 0
            for i in range(len(cross_mat)):
                temp_laplace += cross_mat[j][i]
            laplace_result[j] = temp_laplace/len(cross_mat)
        return laplace_result

    def calculate_column_maxmin(self, cross_mat):
        """
        Berechnet die Maxmin Werte (Minimum je Spalte)
        
        Parameter:
            cross_mat von calculate_cross_efficiency
        Rückgabe:
            Array mit Maxmin Werten pro Spalte
        """
        maxmin_result = np.zeros(len(cross_mat))
        for j in range(len(cross_mat[0])):
            temp_min = float('inf')  # Start with infinity for finding minimum
            for i in range(len(cross_mat)):
                if cross_mat[i][j] < temp_min:
                    temp_min = cross_mat[i][j]
            maxmin_result[j] = temp_min
        return maxmin_result

    def calculate_column_maxmax(self, cross_mat):
        """
        Berechnet die Maxmax Werte (Maximum je Spalte)
        
        Parameter:
            cross_mat von calculate_cross_efficiency
        Rückgabe:
            Array mit Maxmax Werten pro Spalte
        """
        maxmax_result = np.zeros(len(cross_mat))
        for j in range(len(cross_mat[0])):
            temp_max = float('-inf')  # Start with negative infinity for finding maximum
            for i in range(len(cross_mat)):
                if i != j:
                    if cross_mat[i][j] > temp_max:
                        temp_max = cross_mat[i][j]
            maxmax_result[j] = temp_max
        return maxmax_result


    def calculate_column_maverick(self, cross_mat):
        """
        Berechnet die maverick Werte
        
        Formel: (self_efficiency - avg_peer_efficiency) / avg_peer_efficiency
        wobei avg_peer_efficiency der Durchschnitt aller Effizienzwerte für eine DMU
        außer der Selbstevaluation ist.
       
        Parameter:
            cross_mat von calculate_cross_efficiency
        Rückgabe:
            Array mit maverick Werten pro Spalte
        """
        maverick_result = np.zeros(len(cross_mat))
        dmus = len(cross_mat)
        for j in range(len(cross_mat[0])):
            # Reset temp_sum for each DMU j
            temp_sum = 0
            
            # Calculate sum of peer appraisals for DMU j
            for i in range(len(cross_mat)):
                if i != j:  # Exclude self-evaluation
                    temp_sum += cross_mat[j][i]
            
            # Calculate average peer-appraisal and maverick index
            average_peer_appraisal = temp_sum / (dmus - 1)
            maverick_result[j] = (cross_mat[j][j] - average_peer_appraisal) / average_peer_appraisal
        
        return maverick_result