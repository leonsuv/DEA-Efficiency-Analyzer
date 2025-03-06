from nicegui import ui

class ResultsView:
    """View for results panel"""
    
    def __init__(self, controller):
        self.controller = controller
    
    def build(self):
        """Build the results panel UI"""
        ui.label('DEA Analysis Results').classes('text-lg font-bold')
        
        with ui.card().classes('w-full p-4'):
            self.controller.analysis_controller.results_message = ui.label('Please import data first to view results')
            self.controller.analysis_controller.results_tabs = ui.tabs().classes('w-full')
            
            with self.controller.analysis_controller.results_tabs:
                ui.tab('Efficiency Scores')
                ui.tab('Cross-Efficiency Matrix')
                ui.tab('Rankings')
            
            with ui.tab_panels(self.controller.analysis_controller.results_tabs, value='Efficiency Scores').classes('w-full'):
                with ui.tab_panel('Efficiency Scores'):
                    self.controller.analysis_controller.efficiency_table = ui.table(
                        columns=[
                            {'name': 'dmu', 'label': 'DMU', 'field': 'dmu', 'sortable': True},
                            {'name': 'ccr', 'label': 'CCR Score', 'field': 'ccr', 'sortable': True},
                            {'name': 'bcc', 'label': 'BCC Score', 'field': 'bcc', 'sortable': True},
                            {'name': 'scale', 'label': 'Scale Efficiency', 'field': 'scale', 'sortable': True}
                        ],
                        rows=[],
                        row_key='dmu'
                    ).classes('w-full')
                
                with ui.tab_panel('Cross-Efficiency Matrix'):
                    self.controller.analysis_controller.cross_eff_table = ui.table(
                        columns=[],
                        rows=[],
                        row_key='dmu'
                    ).classes('w-full')
                
                with ui.tab_panel('Rankings'):
                    self.controller.analysis_controller.rankings_table = ui.table(
                        columns=[
                            {'name': 'rank', 'label': 'Rank', 'field': 'rank', 'sortable': True},
                            {'name': 'dmu', 'label': 'DMU', 'field': 'dmu', 'sortable': True},
                            {'name': 'avg_ce', 'label': 'Avg Cross-Efficiency', 'field': 'avg_ce', 'sortable': True}
                        ],
                        rows=[],
                        row_key='rank'
                    ).classes('w-full')
            
            with ui.row():
                self.controller.analysis_controller.export_button = ui.button('Export Results to CSV', 
                                                on_click=self.controller.analysis_controller.export_results).props('icon=download')
                self.controller.analysis_controller.export_button.disable()