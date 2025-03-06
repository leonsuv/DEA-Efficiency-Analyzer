import io
import base64
import numpy as np
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

class ChartGenerator:
    """Generate visualization charts for DEA results"""
    
    @staticmethod
    def create_efficiency_chart(dmu_names, ccr_scores, bcc_scores, scale_efficiency):
        """Create efficiency comparison chart"""
        fig = Figure(figsize=(10, 6))
        ax = fig.subplots()
        
        x = np.arange(len(dmu_names))
        width = 0.25
        
        ax.bar(x - width, ccr_scores, width, label='CCR')
        ax.bar(x, bcc_scores, width, label='BCC')
        ax.bar(x + width, scale_efficiency, width, label='Scale Efficiency')
        
        ax.set_ylabel('Efficiency Score')
        ax.set_title('DEA Efficiency Comparison')
        ax.set_xticks(x)
        ax.set_xticklabels(dmu_names, rotation=45, ha='right')
        ax.set_ylim(0, 1.1)
        ax.legend()
        
        fig.tight_layout()
        
        # Convert to base64 image
        return ChartGenerator._fig_to_base64(fig)
    
    @staticmethod
    def create_heatmap(dmu_names, cross_eff_matrix):
        """Create cross-efficiency heatmap"""
        fig = Figure(figsize=(10, 8))
        ax = fig.subplots()
        
        im = ax.imshow(cross_eff_matrix, cmap='viridis')
        
        # Add colorbar
        cbar = ax.figure.colorbar(im, ax=ax)
        
        # Show all ticks and label them
        ax.set_xticks(np.arange(len(dmu_names)))
        ax.set_yticks(np.arange(len(dmu_names)))
        ax.set_xticklabels(dmu_names, rotation=45, ha='right')
        ax.set_yticklabels(dmu_names)
        
        ax.set_title("Cross-Efficiency Matrix Heatmap")
        
        # Loop over data dimensions and create text annotations
        for i in range(len(dmu_names)):
            for j in range(len(dmu_names)):
                text = ax.text(j, i, f"{cross_eff_matrix[i, j]:.2f}",
                              ha="center", va="center", 
                              color="w" if cross_eff_matrix[i, j] < 0.6 else "k")
                
                # Add a special border or highlight for diagonal elements (self-evaluations)
                if i == j:
                    # Draw a rectangle around diagonal elements with red color
                    rect = plt.Rectangle((j-0.5, i-0.5), 1, 1, fill=False, edgecolor='black', linewidth=1)
                    ax.add_patch(rect)
        
        fig.tight_layout()
        
        # Convert to base64 image
        return ChartGenerator._fig_to_base64(fig)
    
    @staticmethod
    def create_rankings_chart(dmu_names, cross_eff_matrix):
        """Create rankings bar chart"""
        # Calculate average cross-efficiency (excluding self)
        avg_cross_eff = np.zeros(len(dmu_names))
        for i in range(len(dmu_names)):
            cross_eff_values = np.delete(cross_eff_matrix[:, i], i)
            avg_cross_eff[i] = np.mean(cross_eff_values)
        
        # Get rankings
        rankings = np.argsort(-avg_cross_eff)
        ranked_dmus = [dmu_names[i] for i in rankings]
        ranked_scores = avg_cross_eff[rankings]
        
        fig = Figure(figsize=(10, 6))
        ax = fig.subplots()
        
        y_pos = np.arange(len(ranked_dmus))
        
        ax.barh(y_pos, ranked_scores, align='center')
        ax.set_yticks(y_pos)
        ax.set_yticklabels(ranked_dmus)
        ax.invert_yaxis()  # labels read top-to-bottom
        ax.set_xlabel('Average Cross-Efficiency Score')
        ax.set_title('DMUs Ranked by Cross-Efficiency')
        
        fig.tight_layout()
        
        # Convert to base64 image
        return ChartGenerator._fig_to_base64(fig)
    
    @staticmethod
    def _fig_to_base64(fig):
        """Convert a matplotlib figure to base64 string"""
        with io.BytesIO() as buf:
            fig.savefig(buf, format='png')
            buf.seek(0)
            img_str = base64.b64encode(buf.read()).decode('utf-8')
            
        return img_str