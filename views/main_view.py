from nicegui import ui
from views.import_view import ImportView
from views.results_view import ResultsView
from views.visualization_view import VisualizationView

class MainView:
    """Main application view"""
    
    def __init__(self, controller):
        self.controller = controller
        self.import_view = None
        self.results_view = None
        self.visualization_view = None
    
    def build(self):
        """Build the main UI layout"""
        
        with ui.tabs().classes('w-full') as tabs:
            ui.tab('Data Import', icon='upload')
            ui.tab('Analysis Results', icon='analytics')
            ui.tab('Visualization', icon='bar_chart')
            ui.tab('About', icon='info')
        
        with ui.tab_panels(tabs, value='Data Import').classes('w-full'):
            with ui.tab_panel('Data Import'):
                self.import_view = ImportView(self.controller)
                self.import_view.build()
            
            with ui.tab_panel('Analysis Results'):
                self.results_view = ResultsView(self.controller)
                self.results_view.build()
            
            with ui.tab_panel('Visualization'):
                self.visualization_view = VisualizationView(self.controller)
                self.visualization_view.build()
            
            with ui.tab_panel('About'):
                self.build_about_panel()
    
    def build_about_panel(self):
        """Build the about panel"""
        with ui.card().classes('w-full p-4'):
            ui.label('About DEA Efficiency Analysis Tool').classes('text-xl font-bold')
            ui.separator()
            
            ui.markdown('''
            ### Data Envelopment Analysis (DEA) Tool
            
            This application allows you to perform DEA efficiency analysis with the following methods:
            
            - **CCR Model** (Constant Returns to Scale)
            - **BCC Model** (Variable Returns to Scale)
            - **Scale Efficiency** Analysis
            - **Cross-Efficiency** Analysis
            
            #### How to use this tool:
            
            1. Import your data through CSV upload or use the sample dataset
            2. View efficiency scores and cross-efficiency results in the Analysis Results tab
            3. Explore visualizations in the Visualization tab
            4. Export your results to CSV for further analysis
            
            #### References:
            
            - Charnes, A., Cooper, W.W., Rhodes, E. (1978). Measuring the efficiency of decision making units. European Journal of Operational Research, 2(6), 429-444.
            - Banker, R.D., Charnes, A., Cooper, W.W. (1984). Some models for estimating technical and scale inefficiencies in data envelopment analysis. Management Science, 30(9), 1078-1092.
            ''')