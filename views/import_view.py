from nicegui import ui
from dealib.dea.utils.options import RTS

class ImportView:
    """View for data import panel"""
    
    def __init__(self, controller):
        self.controller = controller
        self.returns_toggle = None
    
    def build(self):
        """Build the import panel UI"""
        ui.label('Import Data for DEA Analysis').classes('text-lg font-bold')
        
        with ui.card().classes('w-full'):
            ui.label('Import Options').classes('text-md font-bold')
            
            with ui.row():
                ui.radio(['Input-Oriented', 'Output-Oriented'], value='Input-Oriented',
                         on_change=self.controller.change_orientation).props('inline')
            
            with ui.row():
                ui.radio(['Constant Returns (CCR)', 'Variable Returns (BCC)'], value='Constant Returns (CCR)',
                         on_change=self.controller.change_returns_type).props('inline')
            
            ui.separator()
            
            with ui.row():
                with ui.column().classes('w-1/2'):
                    ui.label('Import CSV File with Inputs and Outputs')
                    ui.upload(label='Upload CSV file', on_upload=self.controller.handle_csv_upload).props(
                        'accept=.csv multiple outlined')
                
                with ui.column().classes('w-1/2'):
                    ui.label('Or Use Sample Data')
                    ui.button('Load Sample Dataset', on_click=self.controller.load_sample_data)
            
            ui.separator()
            
            with ui.row().classes('w-full'):
                ui.label('CSV Format Requirements:')
                with ui.column().classes('text-xs'):
                    ui.label('- First column should contain DMU names/IDs')
                    ui.label('- Input variables should be marked with "input_" prefix')
                    ui.label('- Output variables should be marked with "output_" prefix')
                    ui.label('- Example: DMU,input_labor,input_capital,output_revenue,output_profit')
            
            ui.separator()
            
            with ui.expansion('CSV Preview', icon='table_view'):
                self.controller.import_controller.csv_preview = ui.table(
                    columns=[],
                    rows=[],
                    row_key='DMU'
                ).classes('w-full')