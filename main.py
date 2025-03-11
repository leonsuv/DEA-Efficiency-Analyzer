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
from controllers.app_controller import AppController
from views.main_view import MainView
from models.dea_model import DEAModel

# Configure NiceGUI native settings
app.native.window_args['title'] = 'DEA Efficiency Analyzer'
#app.native.window_args['size'] = (1280, 800)
#app.native.window_args['min_size'] = (800, 600)
app.native.settings['ALLOW_DOWNLOADS'] = True


class TabManager:
    """Manages multiple analysis tabs"""
    
    def __init__(self):
        self.tabs = []
        self.next_tab_number = 1   # Always increasing, never decreases
        self.tab_panel = None
        self.tab_headers = None
    
    def initialize(self):
        """Initialize the tab UI structure"""
        ui.add_head_html("""
            <style>
                /* Hide close button unless the tab header has the active class */
                .tab-header:not(.q-tab--active) .close-btn {
                    display: none !important;
                }
            </style>
        """)

        with ui.header().classes('bg-blue-800 text-white'):
            ui.label('DEA Efficiency Analysis Tool').classes('text-xl font-bold')
        
        # Create the tab management bar
        with ui.row().classes('w-full items-center'):
            self.tab_headers = ui.tabs().classes('grow')
            with ui.button('New', on_click=self.add_tab).props('icon=add'):
                ui.tooltip('Create new Analysis')
        
        # Create the tab panel container
        self.tab_panel = ui.tab_panels(self.tab_headers).classes('w-full')
        
        # Add the first tab by default
        self.add_tab()
    
    def close_tab(self, tab_id, tab_header, tab_content):
        if tab_id in self.tabs:
            tab_index = self.tabs.index(tab_id)
            # Remove the tab header and content.
            self.tab_headers.remove(tab_header)
            self.tab_panel.remove(tab_content)
            self.tabs.pop(tab_index)
            
            if self.tabs:
                new_active = self.tabs[tab_index - 1] if tab_index > 0 else self.tabs[0]
                self.tab_headers.set_value(new_active)

    
    def add_tab(self):
        """Add a new analysis tab"""
        tab_id = f"Analysis {self.next_tab_number}"
        self.next_tab_number += 1
        
        # A mutable container to hold tab_content
        tab_content_holder = {}
        
        # Create a new tab header with an embedded close button
        with self.tab_headers:
            # Add custom 'tab-header' class so our CSS can work
            with ui.tab(tab_id).classes('tab-header') as tab_header:
                # For tabs beyond the first, add a close button inside the tab header
                def close_handler(_, tid=tab_id, th=tab_header, tch=tab_content_holder):
                    return self.close_tab(tid, th, tch['content'])
                
                # Add 'close-btn' class to the close button
                with ui.button(icon='close', on_click=close_handler).classes('ml-2 text-xs close-btn').style('width: 10px; height: 10px'):
                    ui.tooltip(f'Close {tab_id}')
        
        # Create tab content with a new instance of the analysis
        with self.tab_panel:
            tab_content = ui.tab_panel(tab_id).classes('w-full')
            tab_content_holder['content'] = tab_content
        
        # Store the tab ID and set as active tab
        self.tabs.append(tab_id)
        self.tab_headers.set_value(tab_id)
        
        # Create a new model-view-controller instance for this tab
        with tab_content:
            model = DEAModel()
            controller = AppController(model=model)
            view = MainView(controller)
            view.build()
            controller.init_views(view)


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