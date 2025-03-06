#
#
# Beispiel "SKigebiete" wird entsprechend aufbereitet.
#  - Ca. 30 DMUs, ca. 4 Inputs, 2 Outputs
# Implementierung 2. Phase
# konstante + variable skalenerträge
#  - (uk= 0)               (uk frei)
#  - Zwischen beiden switchen können
# Zwischen Output und Input orientierung switchen können
# Implementierung der Auswahlregeln
#  - mittlere spaltensumme
#  - Maxi-min u. Maxi-Max (ggf weiter mit spaltensumme)
#  - Maverick-Index (effjj - SUM(efjj)) / SUM(efjj)
#

from nicegui import ui, app

from controllers.app_controller import AppController
from views.main_view import MainView
from models.dea_model import DEAModel

# Configure NiceGUI native settings
app.native.window_args['title'] = 'DEA Efficiency Analyzer'
#app.native.window_args['size'] = (1280, 800)
#app.native.window_args['min_size'] = (800, 600)
app.native.settings['ALLOW_DOWNLOADS'] = True


def main():
    """Main application entry point"""
    # Create model instance
    model = DEAModel()
    
    # Create controller instance with model
    controller = AppController(model=model)
    
    # Create main view with controller
    view = MainView(controller)
    
    # Build the UI
    view.build()
    
    # Initialize controller with view references
    controller.init_views(view)

if __name__ in {"__main__", "__mp_main__"}:
    main()
    # Run the UI in native mode
    ui.run(native=True, reload=False)
    # Run the UI in cloud mode
    # ui.run(on_air="lTod9aB5LoTYQPKu")