import os
import shutil
import subprocess
from datetime import datetime
import logging

class DiagramGenerator:
    def __init__(self, paths):
        self.paths = paths
        self.fig_input_dir = os.path.join(self.paths['output_dir'], 'MEATpy_Figures_Input')
        self.final_figures_dir = os.path.join(self.paths['output_dir'], 'MEATpy_Figures')

    def run(self):
        self.run_biogeochemical_r_script()
        self.draw_sequential_reaction_diagrams()
        self.draw_energy_flow_diagrams()
        
    def run_biogeochemical_r_script(self):
        """
        Run the R script to draw biogeochemical cycle diagrams and organize output.
        """
        r_script = os.path.join(self.paths['r_scripts'], 'draw_biogeochemical_cycles.R')
        input_dir = os.path.join(self.fig_input_dir, 'Nutrient_Cycling_Diagram_Input')
        tmp_output_dir = os.path.join(self.paths['output_dir'], 'tmp')
        final_output_dir = os.path.join(self.final_figures_dir, 'Nutrient_Cycling_Diagrams')

        # Run the R script
        try:
            subprocess.run(
                ["Rscript", r_script, input_dir, tmp_output_dir, "TRUE"],
                check=True,
                stderr=subprocess.DEVNULL  # like 2> /dev/null
            )
        except subprocess.CalledProcessError:
            raise RuntimeError("R script for drawing biogeochemical cycles failed.")

        # Move the generated figures to the final location
        os.makedirs(os.path.dirname(final_output_dir), exist_ok=True)
        if os.path.exists(final_output_dir):
            shutil.rmtree(final_output_dir)
        shutil.move(os.path.join(tmp_output_dir, "draw_biogeochem_cycles"), final_output_dir)

        # Remove temporary output directory
        shutil.rmtree(tmp_output_dir, ignore_errors=True)
        
    def draw_sequential_reaction_diagrams(self):
        """
        Run R script to draw sequential reaction diagrams and move outputs to the final directory.
        
        Args:
            meatpy_dir (str): Path to the MEATpy directory.
            output_dir (str): Path to the output base directory.
            input1_txt (str): Path to Sequential_Transformation_input_1.txt.
            input2_txt (str): Path to Sequential_Transformation_input_2.txt.
            r_mh_tsv (str): Path to draw_sequential_reaction_diagram.R input TSV file.
            r_order_input_1 (str): Path to order input 1 file.
            r_order_input_2 (str): Path to order input 2 file.
        """
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        logging.info("Drawing metabolic handoff diagrams...")

        temp_output_dir = os.path.join(self.paths['output_dir'], 'tmp')
        os.makedirs(temp_output_dir, exist_ok=True)
        os.makedirs(self.final_figures_dir, exist_ok=True)
        
        # Construct Rscript command
        rscript_cmd = [
            "Rscript",
            os.path.join(self.paths['r_scripts'], "draw_sequential_reaction_diagram.R"),
            os.path.join(self.fig_input_dir, 'Sequential_Transformation_input_1.txt'),
            os.path.join(self.fig_input_dir, 'Sequential_Transformation_input_2.txt'),
            os.path.join(self.paths['templates'], 'Sequential_transformations.tsv'),
            os.path.join(self.paths['templates'], 'order_of_input_01.txt'),
            os.path.join(self.paths['templates'], 'order_of_input_02.txt'),
            temp_output_dir
        ]
        subprocess.run(rscript_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        
        # Move output PDFs
        os.rename(
            os.path.join(temp_output_dir, "Bar_plot", "bar_plot_input_1.pdf"),
            os.path.join(self.final_figures_dir, "Sequential_transformation_01.pdf")
        )
        os.rename(
            os.path.join(temp_output_dir, "Bar_plot", "bar_plot_input_2.pdf"),
            os.path.join(self.final_figures_dir, "Sequential_transformation_02.pdf")
        )

        # Cleanup
        shutil.rmtree(temp_output_dir, ignore_errors=True)
        logging.info("Drawing metabolic handoff diagrams finished")
        
    def draw_energy_flow_diagrams(self):
        """
        Runs R scripts to generate metabolic Sankey and functional network diagrams.

        Args:
            scripts_dir (str): Path to the MEATpy directory containing R scripts.
            output_dir (str): Output directory for placing results.
        """

        sankey_input = os.path.join(self.fig_input_dir, 'Metabolic_Sankey_diagram_input.txt')
        sankey_outdir = os.path.join(self.paths['output_dir'], 'Output_energy_flow')
        sankey_final_pdf = os.path.join(self.final_figures_dir, 'Metabolic_Sankey_diagram.pdf')

        fn_input = os.path.join(self.fig_input_dir, 'Functional_network_input.txt')
        fn_outdir = os.path.join(self.paths['output_dir'], 'OutputFolder_Energy')
        fn_final_dir = os.path.join(self.final_figures_dir, 'Functional_network_figures')

        os.makedirs(os.path.dirname(sankey_final_pdf), exist_ok=True)

        subprocess.run([
            "Rscript",
            os.path.join(self.paths['r_scripts'], "draw_metabolic_Sankey_diagram.R"),
            sankey_input,
            sankey_outdir
        ], check=True)

        subprocess.run([
            "mv",
            os.path.join(sankey_outdir, "Energy_plot", "network.plot.pdf"),
            sankey_final_pdf
        ], check=True)

        shutil.rmtree(sankey_outdir, ignore_errors=True)
        
        subprocess.run([
            "Rscript",
            os.path.join(self.paths['r_scripts'], "draw_functional_network_diagram.R"),
            fn_input,
            fn_outdir
        ], check=True)

        # Move and overwrite the output directory
        if os.path.exists(fn_final_dir):
            shutil.rmtree(fn_final_dir)
        shutil.move(os.path.join(fn_outdir, "network_plot"), fn_final_dir)

        # Clean up
        shutil.rmtree(fn_outdir, ignore_errors=True)
        