import numpy as np
from dealib.dea.utils.options import Orientation
from dea.deaeffi import DEAEfficiency

class DEAModel:
    """Data model for DEA analysis"""
    
    def __init__(self):
        self.input_data = None
        self.output_data = None
        self.dmu_names = None
        self.orientation = Orientation.input
        self.dea_results = None
    
    def set_data(self, input_data, output_data, dmu_names=None):
        """Set input and output data for analysis"""
        self.input_data = np.array(input_data)
        self.output_data = np.array(output_data)
        
        if dmu_names is None:
            self.dmu_names = [f"DMU {i+1}" for i in range(len(input_data))]
        else:
            self.dmu_names = dmu_names
    
    def set_orientation(self, orientation):
        """Set DEA orientation"""
        self.orientation = orientation
    
    def run_analysis(self):
        """Run DEA analysis using current data and settings"""
        if self.input_data is None or self.output_data is None:
            raise ValueError("Input and output data must be set before analysis")
            
        self.dea_results = DEAEfficiency(
            self.input_data, 
            self.output_data, 
            orientation=self.orientation,
            instant_calculation=True
        )
        
        return self.dea_results
    
    def get_efficiency_data(self):
        """Get efficiency scores as dict for UI presentation"""
        if self.dea_results is None:
            return []
            
        efficiency_data = []
        for i in range(len(self.dmu_names)):
            efficiency_data.append({
                'dmu': self.dmu_names[i],
                'ccr': round(self.dea_results.ccr_result[i], 4),
                'bcc': round(self.dea_results.bcc_result[i], 4),
                'scale': round(self.dea_results.scale_efficiency[i], 4)
            })
            
        return efficiency_data
    
    def get_cross_efficiency_data(self):
        """Get cross-efficiency matrix data for UI presentation"""
        if self.dea_results is None:
            return [], []
            
        cross_eff_columns = [{'name': 'dmu', 'label': 'DMU', 'field': 'dmu', 'sortable': True}]
        for i in range(len(self.dmu_names)):
            cross_eff_columns.append({
                'name': f'dmu_{i+1}',
                'label': self.dmu_names[i],
                'field': f'dmu_{i+1}',
                'sortable': True
            })
        
        cross_eff_rows = []
        for i in range(len(self.dmu_names)):
            row_data = {'dmu': self.dmu_names[i]}
            for j in range(len(self.dmu_names)):
                row_data[f'dmu_{j+1}'] = round(self.dea_results.cross_eff_matrix[i, j], 4)
            cross_eff_rows.append(row_data)
            
        return cross_eff_columns, cross_eff_rows
    
    def get_rankings_data(self):
        """Get rankings data for UI presentation"""
        if self.dea_results is None:
            return []
            
        # Calculate average cross-efficiency (excluding self-evaluation)
        avg_cross_eff = np.zeros(len(self.dmu_names))
        for i in range(len(self.dmu_names)):
            cross_eff_values = np.delete(self.dea_results.cross_eff_matrix[:, i], i)
            avg_cross_eff[i] = np.mean(cross_eff_values)
        
        # Get rankings
        rankings = np.argsort(-avg_cross_eff)
        
        ranking_data = []
        for i, idx in enumerate(rankings):
            ranking_data.append({
                'rank': i + 1,
                'dmu': self.dmu_names[idx],
                'avg_ce': round(avg_cross_eff[idx], 4),
                'laplace': round(self.dea_results.laplace_values[idx], 4),
                'maxmin': round(self.dea_results.maxmin_values[idx], 4),
                'maxmax': round(self.dea_results.maxmax_values[idx], 4),
                'maverick': round(self.dea_results.maverick_values[idx], 4)
            })
            
        return ranking_data