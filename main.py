import numpy as np
from dea.deaeffi import DEAEfficiency

# Beispieldaten aus der Vorlesung (große Datenmatrix)
x = np.array([
    [1.97, 0.83, 0.4, 1.73],
    [6.73, 4.33, 2.2, 9.66],
    [5.16, 5.64, 4.35, 6.46],
    [4.46, 5.59, 4.96, 5.05],
    [3.42, 4.17, 4.75, 3.71],
    [1.88, 0.85, 0.4, 1.85],
    [7.31, 4.35, 2.2, 10.81],
    [5.94, 6.59, 4.42, 7.58],
    [5.12, 6.16, 5.5, 5.84],
    [5.45, 6.84, 5.34, 6.35],
    [3.92, 3.43, 6.15, 2.75],
    [2.93, 3.32, 5.79, 3.68],
    [2.63, 3.46, 4.32, 2.5],
    [1.71, 1.01, 3.57, 1.1],
    [1.19, 0.79, 0.4, 0.81],
    [6.09, 2.63, 2.2, 5.38],
    [5.32, 5.61, 4.46, 4.52],
    [5.1, 5.34, 5.54, 3.49],
    [4.88, 6.99, 5.42, 3.74],
    [5.87, 7.04, 8.49, 4.53],
    [5.45, 5.64, 9.39, 3.55],
    [3.68, 5.27, 4.88, 2.49],
    [3.77, 4.12, 4.88, 2.41],
])

y = np.array([
    [1.89, 1.02, 2.89],
    [2.7, 6.59, 11.11],
    [4.25, 6.93, 5.76],
    [7.08, 6.39, 5.13],
    [2.23, 1.66, 2.82],
    [4.32, 1.31, 2.8],
    [4.65, 9.34, 12.74],
    [4.92, 12.02, 7.48],
    [5.06, 11.64, 6.07],
    [3.78, 7.89, 5.4],
    [3.98, 3.06, 2.0],
    [1.52, 3.71, 3.18],
    [1.74, 2.62, 2.37],
    [1.66, 1.18, 0.99],
    [4.11, 0.3, 1.28],
    [4.99, 4.52, 6.26],
    [5.73, 6.52, 4.24],
    [9.44, 5.45, 3.59],
    [6.0, 2.48, 3.01],
    [7.96, 4.81, 3.29],
    [3.71, 0.0, 3.06],
    [4.25, 0.0, 2.37],
    [4.05, 0.56, 2.17],
])

# Kleineres Beispiel für einfachere Tests (aus dem Übungsblatt 3)
x2 = np.array([[100], [75], [150], [250], [125], [125]])
y2 = np.array([[100], [50], [125], [100], [75], [75]])

# DEA-Effizienz-Objekt erstellen
dea_analysis = DEAEfficiency(x, y, instant_calculation=True)

# BCC-Werte holen (variable Skaleneffekte)
print("\nBCC:", dea_analysis.bcc_result)

# CCR-Werte holen (konstante Skaleneffekte)
print("\nCCR:", dea_analysis.ccr_result)

# Skaleneffizienz berechnen
print("\nSkaleneffizienz:", dea_analysis.scale_efficiency)

# Kreuzeffizienzanalyse durchführen
print("\nKreuzeffizienz-Matrix:", dea_analysis.cross_eff_matrix)

# Selbstbewertung jeder DMU anzeigen (Diagonale der Kreuzeffizienz-Matrix)
print("\nSelbstbewertung jeder DMU:", dea_analysis.cross_eff_matrix.diagonal())