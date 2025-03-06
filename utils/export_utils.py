import os
import pandas as pd
import numpy as np
from nicegui import ui

class ExportManager:
    """Utility class to handle exporting DEA results to files"""
    
    @staticmethod
    def export_results(model, output_dir="dea_results"):
        """
        Export DEA results to CSV files
        
        Args:
            model: The data model containing DEA results
            output_dir: Directory to save the exported files
        
        Returns:
            bool: True if export was successful, False otherwise
        """
        if not model.dea_results:
            return False
        
        try:
            # Create output directory if it doesn't exist
            os.makedirs(output_dir, exist_ok=True)
            
            # Export efficiency scores
            ExportManager._export_efficiency_scores(model, output_dir)
            
            # Export cross-efficiency matrix
            ExportManager._export_cross_efficiency(model, output_dir)
            
            # Export rankings
            ExportManager._export_rankings(model, output_dir)
            
            return True
        except Exception as e:
            print(f"Export error: {str(e)}")
            return False
    
    @staticmethod
    def _export_efficiency_scores(model, output_dir):
        """Export efficiency scores to CSV"""
        efficiency_df = pd.DataFrame({
            'DMU': model.dmu_names,
            'CCR_Score': model.dea_results.ccr_result,
            'BCC_Score': model.dea_results.bcc_result,
            'Scale_Efficiency': model.dea_results.scale_efficiency
        })
        efficiency_df.to_csv(f"{output_dir}/efficiency_scores.csv", index=False)
    
    @staticmethod
    def _export_cross_efficiency(model, output_dir):
        """Export cross-efficiency matrix to CSV"""
        cross_eff_df = pd.DataFrame(
            model.dea_results.cross_eff_matrix,
            index=model.dmu_names,
            columns=model.dmu_names
        )
        cross_eff_df.to_csv(f"{output_dir}/cross_efficiency_matrix.csv")
    
    @staticmethod
    def _export_rankings(model, output_dir):
        """Export rankings based on cross-efficiency to CSV"""
        avg_cross_eff = np.zeros(len(model.dmu_names))
        for i in range(len(model.dmu_names)):
            cross_eff_values = np.delete(model.dea_results.cross_eff_matrix[:, i], i)
            avg_cross_eff[i] = np.mean(cross_eff_values)
        
        rankings = np.argsort(-avg_cross_eff)
        
        ranking_df = pd.DataFrame({
            'Rank': np.arange(1, len(model.dmu_names) + 1),
            'DMU': [model.dmu_names[i] for i in rankings],
            'Avg_Cross_Efficiency': avg_cross_eff[rankings]
        })
        ranking_df.to_csv(f"{output_dir}/rankings.csv", index=False)