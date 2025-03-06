from controllers.import_controller import ImportController
from controllers.analysis_controller import AnalysisController

class AppController:
    """Main controller that coordinates between specialized controllers"""
    
    def __init__(self, model):
        """Initialize the controller with model reference"""
        self.model = model
        self.main_view = None
        
        # Create specialized controllers
        self.import_controller = ImportController(model)
        self.analysis_controller = AnalysisController(model)
    
    def init_views(self, main_view):
        """Initialize view references after they've been created"""
        self.main_view = main_view
        
        # Pass appropriate view references to specialized controllers
        if hasattr(main_view, 'import_view'):
            self.import_controller.init_view(main_view.import_view)
        
        if hasattr(main_view, 'results_view') and hasattr(main_view, 'visualization_view'):
            self.analysis_controller.init_views(
                main_view.results_view, 
                main_view.visualization_view
            )
    
    # Delegate methods to specialized controllers
    def handle_csv_upload(self, e):
        """Delegate to import controller"""
        success = self.import_controller.handle_csv_upload(e)
        if success:
            self.analysis_controller.update_results_ui()
    
    def load_sample_data(self):
        """Delegate to import controller"""
        success = self.import_controller.load_sample_data()
        if success:
            self.analysis_controller.update_results_ui()
    
    def change_orientation(self, e):
        """Delegate to analysis controller"""
        self.analysis_controller.change_orientation(e)
    
    def update_results_ui(self):
        """Delegate to analysis controller"""
        self.analysis_controller.update_results_ui()
    
    def display_chart(self, chart_type):
        """Delegate to analysis controller"""
        self.analysis_controller.display_chart(chart_type)
    
    def export_results(self):
        """Delegate to analysis controller"""
        self.analysis_controller.export_results()