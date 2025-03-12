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
    # Run the UI in native mode
    ui.run(native=True, reload=False)