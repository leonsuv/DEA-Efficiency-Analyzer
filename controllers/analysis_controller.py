import numpy as np
from nicegui import ui
from dealib.dea.utils.options import Orientation
import base64
import io
import os
import pandas as pd
from utils.export_utils import ExportManager
from utils.chart_utils import ChartGenerator

class AnalysisController:
    """Controller for DEA analysis operations"""
    
    def __init__(self, model):
        """Initialize the analysis controller with model reference"""
        self.model = model
        self.results_view = None
        self.visualization_view = None
        
        # UI element references that will be set by views
        self.efficiency_table = None
        self.cross_eff_table = None
        self.rankings_table = None
        self.vis_message = None
        self.chart_container = None
        self.results_message = None
        self.results_tabs = None
        self.export_button = None
    
    def init_views(self, results_view, visualization_view):
        """Initialize view references after they've been created"""
        self.results_view = results_view
        self.visualization_view = visualization_view
        
        # Get references to UI elements
        if hasattr(self.results_view, 'efficiency_table'):
            self.efficiency_table = self.results_view.efficiency_table
        if hasattr(self.results_view, 'cross_eff_table'):
            self.cross_eff_table = self.results_view.cross_eff_table
        if hasattr(self.results_view, 'rankings_table'):
            self.rankings_table = self.results_view.rankings_table
            
        # Get references to visualization UI elements
        if hasattr(self.visualization_view, 'vis_message'):
            self.vis_message = self.visualization_view.vis_message
        if hasattr(self.visualization_view, 'chart_container'):
            self.chart_container = self.visualization_view.chart_container
    
    def change_orientation(self, e):
        """Change DEA orientation"""
        if e.value == 'Output-Oriented':
            self.model.orientation = Orientation.output
        else:
            self.model.orientation = Orientation.input
            
        if self.model.data_loaded:
            ui.notify('Recalculating with new orientation...', type='info')
            # Re-run analysis with new orientation
            success = self.model.run_analysis()
            if success:
                # Update UI with results
                self.update_results_ui()
    
    def update_results_ui(self):
        """Update UI with latest results from model"""
        # Update efficiency scores table
        efficiency_data = self.model.get_efficiency_data()
        
        if efficiency_data and self.efficiency_table:
            self.efficiency_table.rows = efficiency_data
        
        # Update cross-efficiency matrix table
        cross_eff_columns, cross_eff_rows = self.model.get_cross_efficiency_data()

        if self.export_button:
            self.export_button.enable()
        
        if cross_eff_columns and cross_eff_rows and self.cross_eff_table:
            self.cross_eff_table.columns = cross_eff_columns
            self.cross_eff_table.rows = cross_eff_rows
            
        
        # Update rankings table
        rankings_data = self.model.get_rankings_data()
        
        if rankings_data and self.rankings_table:
            self.rankings_table.rows = rankings_data
        
        # Update visualization messages
        if self.vis_message:
            self.vis_message.text = 'Click on a button to view visualizations'
    
    def display_chart(self, chart_type):
        """Display visualization chart"""
        if not self.model.dea_results:
            ui.notify('Please run analysis first', type='warning')
            return
        
        # Clear previous chart
        if self.chart_container:
            self.chart_container.clear()
        
        if chart_type == 'efficiency':
            self._create_efficiency_chart()
        elif chart_type == 'heatmap':
            self._create_heatmap()
        elif chart_type == 'rankings':
            self._create_rankings_chart()
    
    def _create_efficiency_chart(self):
        """Create efficiency comparison chart"""
        if not self.model.dea_results:
            return
        
        # Use ChartGenerator to create the chart
        img_str = ChartGenerator.create_efficiency_chart(
            self.model.dmu_names, 
            self.model.dea_results.ccr_result, 
            self.model.dea_results.bcc_result, 
            self.model.dea_results.scale_efficiency
        )
        
        # Display the chart
        self._display_chart_image(img_str)

    def _create_heatmap(self):
        """Create cross-efficiency heatmap"""
        if not self.model.dea_results:
            return
        
        # Use ChartGenerator to create the chart
        img_str = ChartGenerator.create_heatmap(
            self.model.dmu_names,
            self.model.dea_results.cross_eff_matrix
        )
        
        # Display the chart
        self._display_chart_image(img_str)

    def _create_rankings_chart(self):
        """Create rankings bar chart"""
        if not self.model.dea_results:
            return
        
        # Use ChartGenerator to create the chart
        img_str = ChartGenerator.create_rankings_chart(
            self.model.dmu_names,
            self.model.dea_results.cross_eff_matrix
        )
        
        # Display the chart
        self._display_chart_image(img_str)

    def _display_chart_image(self, img_str):
        """Helper method to display a base64 image in the chart container"""
        # Clear the chart container first
        if self.chart_container:
            self.chart_container.clear()
            
            # Add the HTML to the container
            with self.chart_container:
                ui.html(f'<img src="data:image/png;base64,{img_str}" style="width:100%;">')
    
    def export_results(self):
        """Export results to CSV files"""
        if not self.model.dea_results:
            ui.notify('No results to export', type='warning')
            return
        
        # Use the ExportManager to handle the export
        success = ExportManager.export_results(self.model)
        
        if success:
            ui.notify(f'Results exported to dea_results folder', type='positive')
        else:
            ui.notify(f'Error exporting results', type='negative')