import argparse
from modules.MEATpy_runner import MEATpyRunner
from version import __version__

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="MEATpy: Microbial Ecosystem Annotation Toolkit (Python)"
    )
    parser.add_argument('-v', '--version', action='version', version=f'MEATpy: Microbial Ecosystem Annotation Toolkit (Python) {__version__} — © 2025 BM Flatgard')
    parser.add_argument('-t', '--threads', type=int, default=20, help='Number of CPU threads to use')
    parser.add_argument('-i', '--input_folder', type=str, help='Input folder containing a "fasta" folder with genome sequences')
    parser.add_argument('-csv', '--input_csv', type=str, help='Input CSV file containing read pairs')
    parser.add_argument('-o', '--output_folder', type=str, default='./MEATpy_out', help='Output directory')
    #parser.add_argument('-n', '--project_name', type=str, default='MEATpy', help='Project name for output files. WARNING: Using the same name will overwrite previous results in using the same output directory!')
    parser.add_argument('-ko_db', '--kofam_database_size', type=str, default='full', choices=['small', 'full'], help='Size of the KOfam database to use')
    parser.add_argument('-m', '--module_cutoff', type=float, default=0.75, help='Cutoff value for module presence')
    parser.add_argument('-pm', '--prodigal_method', type=str, default='meta', choices=['meta', 'single'], help='Prodigal method to annotate ORFs')
    parser.add_argument('-rt', '--reads_type', type=str, default='metaG', choices=['metaG', 'metaT'], help='Type of omic reads')
    parser.add_argument('-st', '--sequencing_type', type=str, default='illumina', choices=['illumina', 'pacbio', 'pacbio_hifi', 'pacbio_asm20', 'nanopore'], help='Sequencing type of the reads')
    parser.add_argument('-tax', '--taxonomy', type=str, default='phylum', choices=['phylum', 'class', 'order', 'family', 'genus', 'species', 'bin'], help='Taxonomic level for MW-score calculation')
    parser.add_argument('-c', '--correction', type=str, default='GTDB', choices=['GTDB', 'NCBI'], help='Correction type for GTDB-Tk to NCBI')
    parser.add_argument("--skip_to_diagrams", action="store_true", help="Start from diagram input step.")
    parser.add_argument("--rename", action="store_true", help="Rename input MAGs from nf-core/mag using GTDBtk/NCBI Classification")
    parser.add_argument("--revert_names", action="store_true", help="Revert names of input MAGs to original names")
    parser.add_argument("--oops", action="store_true", help="For those times you forget to add --rename and you did'n notice until the end...")

    return parser.parse_args()

def main():
    args = parse_arguments()
    runner = MEATpyRunner(args)
    runner.run()

if __name__ == "__main__":
    main()
