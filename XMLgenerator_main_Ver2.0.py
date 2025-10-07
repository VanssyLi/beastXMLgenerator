import pandas as pd
import module1_Taxa as mTaxa
import module2_Seq as mSeq
import module3_treeModel as mTreeModel
import module4_popModel as mPop
import module5_clockModel as mClock
import module6_subsModel as mSubs
import module7_likelihood as mLikeli
import module8_operators as mOperator
import module9_priors as mPriors
import module10_logs as mLog
import functions as f
import math
import xml.etree.ElementTree as ET
import os
import re
import argparse
import sys

# ============================================================================
# DIRECT RUN CONFIGURATION - 直接运行配置
# ============================================================================
# Set your config file path here for direct execution
# DEFAULT_CONFIG_FILE = '/Users/hmj759/Documents/PhD/Scripts/TipDating/Belugas/196modern+26ancient+1narwal_belugas_joint.config'
# DEFAULT_FASTA_FILE = "/Users/hmj759/Documents/PhD/Scripts/TipDating/Belugas/196modern_26ancient_1narwal_belugas_mitogenomes.joint.fasta"

# Set to True for direct run mode, False for command line mode
DIRECT_RUN_MODE = False
# ============================================================================


def parse_config(config_file):
    """Parse configuration file and extract key-value pairs."""
    with open(config_file, 'r') as f:
        content = f.read()
    
    config = {}
    # Single regex pattern for all key-value pairs
    pattern = r'^(?!#)\s*(\w+)\s*=\s*["\']?([^"\'#\n]+?)["\']?\s*(?:#.*)?$'
    
    for match in re.finditer(pattern, content, re.MULTILINE):
        key, value = match.groups()
        value = value.strip()
        
        # Type conversion
        if key == 'split_partition':
            config[key] = value.lower() == 'true'
        elif key in ['root_height_mean', 'root_height_stdev', 'root_height_offset', 'chainLength', 'log_every']:
            config[key] = float(value)
        else:
            config[key] = value
    
    return config


def read_metadata(metadata_file, fasta):
    """Extract taxon sets and root height lower bound from metadata."""
    df = pd.read_excel(metadata_file)
    
    fasta_taxa = set(f.get_taxa_name(fasta))
    
    # Get unique taxon groups and filter by fasta content
    taxon_set = []
    for group in df['Group-By'].unique():
        has_sequence = any(str(group) in taxon for taxon in fasta_taxa)
        if has_sequence:
            taxon_set.append([group])
        else:
            print(f"Warning: Skipping group '{group}' - no sequences found in fasta file")
    
    # Get max calibrated age (numeric only)
    root_height_lower = pd.to_numeric(df['Calibrated_yBP'], errors='coerce').max()
    if pd.isna(root_height_lower):
        root_height_lower = 0.0
    
    return taxon_set, root_height_lower


def create_partition_file(config, output_dir):
    """Generate partition file based on split_partition setting."""
    if config['split_partition']:
        # Parse multi-partition model: "tRNA:HKY+G+I, rRNA:HKY+G+I, ..."
        partitions = [p.split(':') for p in config['subs_model'].split(', ')]
        data = [{'Partition': name, 'Exclude': False, 'Every3': name == 'CDS', 'SubstitutionModel': model}
                for name, model in partitions]
    else:
        # Single partition
        data = [{'Partition': 'all', 'Exclude': False, 'Every3': False, 'SubstitutionModel': config['subs_model']}]
    
    partition_file = os.path.join(output_dir, 'temp_partition.xlsx')
    pd.DataFrame(data).to_excel(partition_file, index=False)
    return partition_file


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate BEAST XML file from configuration",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -c config.config -f alignment.fasta
  %(prog)s -c config.config -f alignment.fasta -o output_dir/

The config file should contain all necessary parameters for XML generation.
The fasta alignment file must be provided via -f flag.
        """
    )
    
    parser.add_argument('-c', '--config', required=True,
                       help='Path to configuration file')
    parser.add_argument('-f', '--fasta', required=True,
                       help='Path to fasta alignment file')
    parser.add_argument('-o', '--output-dir',
                       help='Output directory (default: same as config file)')
    
    return parser.parse_args()


def get_config_path():
    """Get config file path either from command line args or direct configuration."""
    if DIRECT_RUN_MODE:
        # Direct run mode - use configured path
        class DirectArgs:
            def __init__(self, config, fasta=None, output_dir=None):
                self.config = config
                self.fasta = fasta
                self.output_dir = output_dir
        
        # return DirectArgs(DEFAULT_CONFIG_FILE, DEFAULT_FASTA_FILE)
    else:
        # Command line mode - parse arguments
        return parse_arguments()


def main():
    """Main function."""
    args = get_config_path()
    
    print("BEAST XML Generator Ver 2.0")
    if DIRECT_RUN_MODE:
        print("Running in DIRECT MODE")
    print(f"Config file: {args.config}")
    print()
    
    # Validate config file
    if not os.path.exists(args.config):
        print(f"Error: Config file not found: {args.config}")
        sys.exit(1)
    
    # Parse configuration
    config = parse_config(args.config)
    
    # Set output directory
    if hasattr(args, 'output_dir') and args.output_dir:
        output_dir = args.output_dir
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
    else:
        output_dir = os.path.dirname(args.config)
    
    # Get fasta from command line/direct config (NOT from config file)
    fasta = args.fasta
    if not fasta:
        print("Error: Fasta file must be provided via -f flag or DEFAULT_FASTA_FILE")
        sys.exit(1)
    
    # Extract fasta prefix (filename without extension)
    fasta_basename = os.path.basename(fasta)
    fasta_prefix = os.path.splitext(fasta_basename)[0]
    
    # Get output directory from config and ensure it ends with /
    outdir = config.get('outdir', '')
    if outdir and not outdir.endswith('/'):
        outdir += '/'
    
    # Generate log_name: outdir + fasta_prefix
    log_name = fasta_prefix
    
    priors_table = config['tip_priors']
    gff = config.get('annotation', '')
    split_partition = config['split_partition']
    
    # Validate required files
    required_files = [fasta, priors_table, config['metadata']]
    for file_path in required_files:
        if not os.path.exists(file_path):
            print(f"Error: Required file not found: {file_path}")
            sys.exit(1)
    
    print(f"Processing alignment: {fasta}")
    print(f"Using metadata: {config['metadata']}")
    print(f"Output prefix: {log_name}")
    
    # Read metadata and create partition file - NOW PASSING FASTA
    taxon_set, root_height_lower = read_metadata(config['metadata'], fasta)
    partition_file = create_partition_file(config, output_dir)
    
    # BEAST settings
    beast_setting = f.Xml(
        units='years',
        population_model=config['population_model'],
        root_height_mean=config['root_height_mean'],
        root_height_stdev=config['root_height_stdev'],
        root_height_offset=config['root_height_offset'],
        root_height_lower=root_height_lower
    )
    
    # MCMC parameters
    chainLength = int(config['chainLength'])
    log_every = int(config['log_every'])
    
    # Calculate skygrid cutoff
    magnitude = 10 ** (len(str(int(beast_setting.root_height_mean))) - 1)
    skygrid_cutoff = math.ceil(beast_setting.root_height_mean / magnitude) * magnitude
    
    # Read substitution model from partition file
    with open(partition_file, 'rb') as p:
        subs_model = pd.read_excel(p, header=0)['SubstitutionModel'].iloc[0]
    
    print("Building XML structure...")
    
    # Module 1 -- Introduce all the taxa and taxon sets
    beast = ET.Element('beast')
    tree = ET.ElementTree(beast)
    mTaxa.intro_taxa(beast, fasta, priors_table)
    mTaxa.intro_taxonSet(beast, taxon_set, fasta)
    
    # Module 2 -- Introduce aligned sequences with/without partitions
    mSeq.build_partitionSeq(beast, partition_file, fasta, gff, split_partition)
    mSeq.build_patterns(beast, partition_file, split_partition)
    
    # Module 3 -- Introduce the tree model
    mTreeModel.build_initialTree(beast, units=beast_setting.units)
    mTreeModel.build_treeModel(beast, fasta)
    
    mTreeModel.calc_tip_branchLen(beast)
    mTreeModel.calc_taxonSet_branchLen(beast, taxon_set)
    
    # Module 4 -- Define the population model
    mPop.skygrid(beast, mean=beast_setting.root_height_mean)
    
    # Module 5 -- Define the clock model
    mClock.strict_clock(beast, partition_file, split_partition)
    
    # Module 6 -- Define the substitution model
    if 'HKY' in subs_model:
        mSubs.HKY(beast, partition_file, split_partition)
    elif 'GTR' in subs_model:
        mSubs.GTR(beast, partition_file, split_partition)
    
    # Module 7 -- Likelihoods
    mLikeli.calc_likelihood(beast, partition_file, split_partition)
    
    # Module 8 -- Operators
    mOperator.build_operator(beast, partition_file, fasta, split_partition)
    
    # Define MCMC
    MCMC = ET.SubElement(beast, 'mcmc', attrib={'id': 'mcmc', 'chainLength': str(chainLength), 'autoOptimize': 'true'})
    joint = ET.SubElement(MCMC, 'joint', attrib={'id': 'joint'})
    jointPrior = ET.SubElement(joint, 'prior', attrib={'id': 'prior'})
    f.make_comment(beast, ' Define MCMC ')
    
    # Module 9 -- Priors
    mPriors.define_operators(jointPrior, joint, MCMC, fasta, priors_table, partition_file, split_partition,
                             rootHeight_mean=beast_setting.root_height_mean,
                             rootHeight_stdev=beast_setting.root_height_stdev,
                             rootHeight_offset=beast_setting.root_height_offset)
    
    # Module 10 -- Write Logs
    mLog.to_screen(MCMC, partition_file, log_every, split_partition)
    mLog.to_file(MCMC, log_every, log_name, taxon_set, partition_file, fasta, split_partition)
    mLog.trees_to_file(MCMC, log_every, log_name, partition_file, split_partition)
    mLog.report(beast)
    
    # Test and store the XML file
    f.pretty_xml(beast)
    
    # Generate output XML file path
    if output_dir != os.path.dirname(args.config):
        xml_output = os.path.join(output_dir, f'{log_name}.xml')
    else:
        xml_output = f'{log_name}.xml'
    
    tree.write(xml_output, encoding='utf-8', xml_declaration=True)
    
    # Clean up fake partition file
    os.remove(partition_file)
    
    print(f"\nXML file generated successfully: {xml_output}")
    print("Processing complete!")


if __name__ == "__main__":
    main()




