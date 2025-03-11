#
#
# Beispiel "Skigebiete" wird entsprechend aufbereitet.
#  - Ca. 30 DMUs, ca. 4 Inputs, 2 Outputs
# Implementierung 2. Phase
# konstante + variable skalenerträge
#  - (uk= 0)               (uk frei)
#  - Zwischen beiden switchen können
# Zwischen Output und Input orientierung switchen können (Check)
# Implementierung der Auswahlregeln (Check)
#  - mittlere spaltensumme
#  - Maxi-min u. Maxi-Max (ggf weiter mit spaltensumme)
#  - Maverick-Index (effjj - SUM(efjj)) / SUM(efjj)
#

from nicegui import ui, app
from views.tab_manager import TabManager

# Configure NiceGUI native settings
app.native.window_args['title'] = 'DEA Efficiency Analyzer'
#app.native.window_args['size'] = (1280, 800)
#app.native.window_args['min_size'] = (800, 600)
app.native.settings['ALLOW_DOWNLOADS'] = True


def main():
    """Main application entry point"""
    tab_manager = TabManager()
    tab_manager.initialize()


if __name__ in {"__main__", "__mp_main__"}:
    main()
    # Run the UI in native mode
    ui.run(native=True, reload=False)
    # Run the UI in cloud mode
    # ui.run(on_air="lTod9aB5LoTYQPKu")