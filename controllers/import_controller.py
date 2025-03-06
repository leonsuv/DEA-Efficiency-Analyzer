import io
import pandas as pd
import numpy as np
from nicegui import ui

class ImportController:
    """Controller for data import operations"""
    
    def __init__(self, model):
        """Initialize the import controller with model reference"""
        self.model = model
        self.import_view = None
        self.csv_preview = None
        
    def init_view(self, import_view):
        """Initialize view reference after it's been created"""
        self.import_view = import_view
    
    def handle_csv_upload(self, e):
        """Handle CSV file upload"""
        try:
            content = io.StringIO(e.content.read().decode('utf-8'))
            df = pd.read_csv(content)
            
            # Validate the CSV format
            if len(df.columns) < 3:  # At least DMU, 1 input, 1 output
                ui.notify('CSV must have at least 3 columns (DMU, input, output)', type='negative')
                return
            
            # Extract DMU names, inputs, and outputs
            dmu_names = df.iloc[:, 0].values.tolist()
            
            # Extract inputs and outputs based on column prefixes
            input_cols = [col for col in df.columns if col.startswith('input_')]
            output_cols = [col for col in df.columns if col.startswith('output_')]
            
            if not input_cols or not output_cols:
                ui.notify('CSV columns must be prefixed with "input_" or "output_"', type='negative')
                return
            
            # Convert to numpy arrays
            input_data = df[input_cols].values
            output_data = df[output_cols].values
            
            # Set data in model and run analysis
            self.model.set_data(input_data, output_data, dmu_names)
            self.model.set_orientation(self.model.orientation)
            success = self.model.run_analysis()

            # Set Preview of imported data
            preview_df = df.head(10)
            self.csv_preview.columns = [{'name': col, 'label': col, 'field': col} for col in preview_df.columns]
            self.csv_preview.rows = preview_df.to_dict('records')
            
            if success:
                ui.notify('Data loaded successfully!', type='positive')
                # Signal that data is loaded and analysis is complete
                return True
            else:
                ui.notify('Error loading data into model', type='negative')
                return False
        
        except Exception as ex:
            ui.notify(f'Error loading CSV: {str(ex)}', type='negative')
            print(f"Error: {str(ex)}")
            return False
    
    def load_sample_data(self):
        """Load sample dataset"""
        # Sample data inputs (simplified version)
        input_data = np.array([
            [1.97, 0.83, 0.4, 1.73],
            [6.73, 4.33, 2.2, 9.66],
            [5.16, 5.64, 4.35, 6.46],
            [4.46, 5.59, 4.96, 5.05],
            [3.42, 4.17, 4.75, 3.71],
            [1.88, 0.85, 0.4, 1.85],
            [7.31, 4.35, 2.2, 10.81],
            [5.94, 6.59, 4.42, 7.58],
            [5.12, 6.16, 5.5, 5.84],
            [5.45, 6.84, 5.34, 6.35],
            [3.92, 3.43, 6.15, 2.75],
            [2.93, 3.32, 5.79, 3.68],
            [2.63, 3.46, 4.32, 2.5],
            [1.71, 1.01, 3.57, 1.1],
            [1.19, 0.79, 0.4, 0.81],
        ])
        
        # Sample data outputs
        output_data = np.array([
            [1.89, 1.02, 2.89],
            [2.7, 6.59, 11.11],
            [4.25, 6.93, 5.76],
            [7.08, 6.39, 5.13],
            [2.23, 1.66, 2.82],
            [4.32, 1.31, 2.8],
            [4.65, 9.34, 12.74],
            [4.92, 12.02, 7.48],
            [5.06, 11.64, 6.07],
            [3.78, 7.89, 5.4],
            [3.98, 3.06, 2.0],
            [1.52, 3.71, 3.18],
            [1.74, 2.62, 2.37],
            [1.66, 1.18, 0.99],
            [4.11, 0.3, 1.28],
        ])
        
        dmu_names = [f"DMU {i+1}" for i in range(len(input_data))]
    
        self.model.set_data(input_data, output_data, dmu_names)
        success = self.model.run_analysis()

        # Set Preview Data
        # Create preview dataframe
        preview_df = pd.DataFrame()
        preview_df['DMU'] = dmu_names
        
        for i in range(input_data.shape[1]):
            preview_df[f'input_{i+1}'] = input_data[:, i]
        
        for i in range(output_data.shape[1]):
            preview_df[f'output_{i+1}'] = output_data[:, i]
        
        # Update preview table
        self.csv_preview.columns = [{'name': col, 'label': col, 'field': col} for col in preview_df.columns]
        self.csv_preview.rows = preview_df.to_dict('records')
        
        if success:
            ui.notify('Sample data loaded successfully!', type='positive')
            return True
        else:
            ui.notify('Error loading sample data', type='negative')
            return False