from nicegui import ui

class VisualizationView:
    """View for visualization panel"""
    
    def __init__(self, controller):
        self.controller = controller
    
    def build(self):
        """Build the visualization panel UI"""
        ui.label('DEA Visualization').classes('text-lg font-bold')
        
        with ui.card().classes('w-full p-4'):
            self.controller.analysis_controller.vis_message = ui.label('Please import data first to view visualizations')
            
            with ui.row():
                ui.button('Efficiency Comparison Chart', 
                         on_click=lambda: self.controller.analysis_controller.display_chart('efficiency'))
                ui.button('Cross-Efficiency Heatmap', 
                         on_click=lambda: self.controller.analysis_controller.display_chart('heatmap'))
                ui.button('Rankings Chart', 
                         on_click=lambda: self.controller.analysis_controller.display_chart('rankings'))
            
            self.controller.analysis_controller.chart_container = ui.element('div').classes('w-full h-128 flex items-center justify-center')