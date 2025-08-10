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
        TOL = 1e-12
        cross_eff_matrix = np.zeros((self.n_dmus, self.n_dmus))
        for i in range(self.n_dmus):
            if self.rts == RTS.crs:
                u, v = self._solve_second_phase_ccr(i)
                if u is None or v is None:
                    warnings.warn(f"CCR LP for DMU {i} did not converge; row left as zeros.")
                    continue
                for j in range(self.n_dmus):
                    denom = float(np.dot(v, self.x[j]))
                    eff = 0.0 if denom <= TOL else float(np.dot(u, self.y[j])) / denom
                    # cap tiny numerical overshoot
                    if eff > 1.0 and eff < 1.0 + 1e-9:
                        eff = 1.0
                    cross_eff_matrix[i, j] = eff if float_values else int(round(eff * 100))
            else:
                # BCC (VRS)
                u, v, u0 = self._solve_second_phase_bcc(i)
                if u is None or v is None:
                    warnings.warn(f"BCC LP for DMU {i} did not converge; row left as zeros.")
                    continue
                for j in range(self.n_dmus):
                    if self.orientation == Orientation.output:
                        # Use denominator shift to avoid negative numerators
                        denom = float(np.dot(v, self.x[j])) - float(u0)
                        num = float(np.dot(u, self.y[j]))
                    else:
                        # Input‑oriented BCC (original formula is fine)
                        denom = float(np.dot(v, self.x[j]))
                        num = float(np.dot(u, self.y[j])) + float(u0)
                    eff = 0.0 if denom <= TOL else num / denom
                    # confine to [0,1] with small tolerance
                    if eff < 0.0 and eff > -1e-9:
                        eff = 0.0
                    if eff > 1.0 and eff < 1.0 + 1e-9:
                        eff = 1.0
                    cross_eff_matrix[i, j] = eff if float_values else int(round(eff * 100))
        return cross_eff_matrix


    def _solve_second_phase_ccr(self, dmu_index):
        """
        Standard CCR multiplier model (no slacks):
        - Output-oriented:
            min v^T x_k
            s.t. u^T y_k = 1
                u^T y_j - v^T x_j ≤ 0, ∀j
                u ≥ 0, v ≥ 0
        - Input-oriented:
            max u^T y_k  <=>  min -u^T y_k
            s.t. v^T x_k = 1
                u^T y_j - v^T x_j ≤ 0, ∀j
                u ≥ 0, v ≥ 0
        """
        s_dim, m_dim, n = self.s, self.m, self.n_dmus
        num_variables = s_dim + m_dim

        # Small positive lower bounds to stabilize solution and avoid zero denominators
        eps = 1e-9
        bounds = [(eps, None)] * s_dim + [(eps, None)] * m_dim

        # Common inequality constraints: u^T y_j - v^T x_j ≤ 0
        A_ub = []
        b_ub = []
        for j in range(n):
            row = [0.0] * num_variables
            # u part
            for i in range(s_dim):
                row[i] = float(self.y[j, i])
            # v part (-x)
            for i in range(m_dim):
                row[s_dim + i] = -float(self.x[j, i])
            A_ub.append(row)
            b_ub.append(0.0)

        if self.orientation == Orientation.output:
            # min v^T x_k
            c = [0.0] * s_dim + [float(self.x[dmu_index, i]) for i in range(m_dim)]
            # u^T y_k = 1
            A_eq = [[0.0] * num_variables]
            for i in range(s_dim):
                A_eq[0][i] = float(self.y[dmu_index, i])
            b_eq = [1.0]
        else:
            # min -u^T y_k  (i.e., maximize u^T y_k)
            c = [-float(self.y[dmu_index, i]) for i in range(s_dim)] + [0.0] * m_dim
            # v^T x_k = 1
            A_eq = [[0.0] * num_variables]
            for i in range(m_dim):
                A_eq[0][s_dim + i] = float(self.x[dmu_index, i])
            b_eq = [1.0]

        res = opt.linprog(c=c, A_ub=A_ub, b_ub=b_ub,
                        A_eq=A_eq, b_eq=b_eq, bounds=bounds,
                        method="highs")
        if res.success:
            sol = res.x
            u = sol[:s_dim]
            v = sol[s_dim:s_dim + m_dim]
            return u, v
        # Fallback: try highs-ipm if highs fails
        res = opt.linprog(c=c, A_ub=A_ub, b_ub=b_ub,
                        A_eq=A_eq, b_eq=b_eq, bounds=bounds,
                        method="highs-ipm")
        if res.success:
            sol = res.x
            u = sol[:s_dim]
            v = sol[s_dim:s_dim + m_dim]
            return u, v
        return None, None


    def _solve_second_phase_bcc(self, dmu_index):
        """
        Standard BCC (VRS) multiplier model (u0 free):
        Constraints (all orientations):
            u^T y_j - v^T x_j + u0 ≤ 0, ∀j
            u ≥ 0, v ≥ 0, u0 free

        - Output-oriented:
            min v^T x_k - u0
            s.t. u^T y_k = 1

        - Input-oriented:
            max u^T y_k + u0 <=> min -(u^T y_k + u0)
            s.t. v^T x_k = 1
        """
        s_dim, m_dim, n = self.s, self.m, self.n_dmus
        num_variables = s_dim + m_dim + 1  # u, v, u0

        eps = 1e-9
        # u ≥ eps, v ≥ eps, u0 free
        bounds = [(eps, None)] * s_dim + [(eps, None)] * m_dim + [(None, None)]

        # Inequalities: u^T y_j - v^T x_j + u0 ≤ 0
        A_ub = []
        b_ub = []
        for j in range(n):
            row = [0.0] * num_variables
            for i in range(s_dim):
                row[i] = float(self.y[j, i])
            for i in range(m_dim):
                row[s_dim + i] = -float(self.x[j, i])
            row[s_dim + m_dim] = 1.0  # + u0
            A_ub.append(row)
            b_ub.append(0.0)

        if self.orientation == Orientation.output:
            # min v^T x_k - u0
            c = [0.0] * s_dim + [float(self.x[dmu_index, i]) for i in range(m_dim)] + [-1.0]
            # u^T y_k = 1
            A_eq = [[0.0] * num_variables]
            for i in range(s_dim):
                A_eq[0][i] = float(self.y[dmu_index, i])
            b_eq = [1.0]
        else:
            # min -(u^T y_k + u0)
            c = [-float(self.y[dmu_index, i]) for i in range(s_dim)] + [0.0] * m_dim + [-1.0]
            # v^T x_k = 1
            A_eq = [[0.0] * num_variables]
            for i in range(m_dim):
                A_eq[0][s_dim + i] = float(self.x[dmu_index, i])
            b_eq = [1.0]

        res = opt.linprog(c=c, A_ub=A_ub, b_ub=b_ub,
                        A_eq=A_eq, b_eq=b_eq, bounds=bounds,
                        method="highs")
        if res.success:
            sol = res.x
            u = sol[:s_dim]
            v = sol[s_dim:s_dim + m_dim]
            u0 = sol[s_dim + m_dim]
            return u, v, u0
        # Fallback solver
        res = opt.linprog(c=c, A_ub=A_ub, b_ub=b_ub,
                        A_eq=A_eq, b_eq=b_eq, bounds=bounds,
                        method="highs-ipm")
        if res.success:
            sol = res.x
            u = sol[:s_dim]
            v = sol[s_dim:s_dim + m_dim]
            u0 = sol[s_dim + m_dim]
            return u, v, u0
        return None, None, None



    def calculate_column_laplace(self, cross_mat):
        """
        Column averages (Laplace) for each DMU (i.e., average appraisal received).
        """
        n = len(cross_mat)
        laplace_result = np.zeros(n)
        for j in range(n):
            col_sum = 0.0
            for i in range(n):
                col_sum += cross_mat[i][j]  # correct: across rows i, fixed column j
            laplace_result[j] = col_sum / n
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
        Maverick index:
        (self_evaluation - average_peer_appraisal) / average_peer_appraisal
        where average_peer_appraisal is the column average excluding the diagonal.
        """
        n = len(cross_mat)
        maverick_result = np.zeros(n)
        for j in range(n):
            peer_sum = 0.0
            for i in range(n):
                if i != j:
                    peer_sum += cross_mat[i][j]  # correct: peer appraisals about DMU j
            avg_peer = peer_sum / (n - 1)
            self_eval = cross_mat[j][j]
            # avoid division by zero
            if avg_peer == 0:
                maverick_result[j] = np.nan
            else:
                maverick_result[j] = (self_eval - avg_peer) / avg_peer
        return maverick_result