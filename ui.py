from nicegui import ui
import numpy as np
from dea.deaeffi import DEAEfficiency

# Funktion zum Laden und Anzeigen der DEA-Daten
def daten_laden():
    # Beispieldaten aus dem Übungsblatt (wie in main.py)
    x = np.array([[100], [75], [150], [250], [125], [125]])
    y = np.array([[100], [50], [125], [100], [75], [75]])
    
    # DEA-Analyse erstellen
    dea_analyse = DEAEfficiency(x, y)
    
    # CCR und BCC Werte holen
    ccr = dea_analyse.get_ccr()
    bcc = dea_analyse.get_bcc()
    
    # Skaleneffizienz berechnen (CCR/BCC)
    scale = ccr / bcc
    
    # Tabellendaten erstellen
    tabellen_daten = [
        {
            'DMU': f'DMU {i+1}', 
            'CCR': f'{c:.4f}', 
            'BCC': f'{b:.4f}',
            'Scale': f'{s:.4f}'
        } 
        for i, (c, b, s) in enumerate(zip(ccr, bcc, scale))
    ]
    
    # Vorherige Tabelle löschen, falls vorhanden
    if hasattr(daten_laden, 'tabelle'):
        daten_laden.tabelle.clear()
    
    # Überschrift
    ui.label('DEA-Analyse').classes('text-h5 text-weight-bold')
    
    # Neue Tabelle anzeigen
    daten_laden.tabelle = ui.table(columns=[
        {'name': 'dmu', 'label': 'DMU', 'field': 'DMU', 'align': 'left'},
        {'name': 'ccr', 'label': 'CCR', 'field': 'CCR', 'align': 'right'},
        {'name': 'bcc', 'label': 'BCC', 'field': 'BCC', 'align': 'right'},
        {'name': 'scale', 'label': 'Skaleneffizienz', 'field': 'Scale', 'align': 'right'}
    ], rows=tabellen_daten).classes('w-full')
    
    # Erklärung zur DEA
    with ui.card().classes('w-full'):
        ui.label('Hinweis zur Interpretation:').classes('text-weight-bold')
        ui.label('• CCR: Constant Returns to Scale')
        ui.label('• BCC: Variable Returns to Scale')
        ui.label('• Skaleneffizienz: CCR/BCC - Gibt an, ob eine DMU in optimaler Größe operiert')

# Hauptlayout erstellen
with ui.card().classes('w-full'):
    ui.label('DEA-Rechner').classes('text-h4')
    ui.label('DDE-Teamaufgabe').classes('text-subtitle1')
    ui.separator()
    ui.button('Beispieldaten laden', on_click=daten_laden).classes('text-white bg-primary')

ui.run(native=True, reload=False, window_size=[1000, 900])