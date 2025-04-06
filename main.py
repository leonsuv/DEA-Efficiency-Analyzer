from nicegui import ui, app
from views.tab_manager import TabManager

app.native.window_args['title'] = 'DEA Efficiency Analyzer'
app.native.settings['ALLOW_DOWNLOADS'] = True

def main():
    """Main application entry point"""
    tab_manager = TabManager()
    tab_manager.initialize()

if __name__ in {"__main__", "__mp_main__"}:
    main()
    
    # Set this to True to run in container mode
    # Impportant: When set to True, use the "run.sh" script to start the app
    run_in_container = False
    
    if run_in_container:
        # Run in web mode when in container
        ui.run(host='0.0.0.0', port=8080, native=False, reload=False)
    else:
        # Run the UI in native mode
        ui.run(native=True, reload=False)