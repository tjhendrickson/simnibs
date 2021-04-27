#!/usr/bin/env python3
import argparse
import os
import subprocess
import nibabel
import numpy
from glob import glob
from simulation import SimNIBS_simulation

__version__ = open(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                'version')).read()

def run(command, env={}):
    merged_env = os.environ
    merged_env.update(env)
    process = subprocess.Popen(command, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, shell=True,
                               env=merged_env)
    while True:
        line = process.stdout.readline()
        line = str(line, 'utf-8')[:-1]
        print(line)
        if line == '' and process.poll() != None:
            break
    if process.returncode != 0:
        raise Exception("Non zero return code: %d"%process.returncode)

parser = argparse.ArgumentParser(description='Example BIDS App entrypoint script.')
parser.add_argument('bids_dir', help='The directory with the input dataset '
                    'formatted according to the BIDS standard.')
parser.add_argument('output_dir', help='The directory where the output files '
                    'should be stored. If you are running group level analysis '
                    'this folder should be prepopulated with the results of the'
                    'participant level analysis.')
parser.add_argument('analysis_level', help='Level of the analysis that will be performed. '
                    'Multiple participant level analyses can be run independently '
                    '(in parallel) using the same output_dir.',
                    choices=['participant', 'group'],default='participant')
parser.add_argument('--coil_name',help='What coil model to use based on coil file.', 
                    choices=['Magstim_70mm_fig8', 'MagVenture_MC_B70', 'No20_Numerical_Helmholtz',
                        'No37_Magstim_double_cone', 'No21_Three_Layer_Double_Coil', 'No49_Numerical_Maxwell'
                        'No25_Magstim_Figure8_25mm', 'No4_Magstim_circular_70mm',
                        'No29_MagVenture_C-B60_Fig8', 'No50_Numerical_Golay', 'No30_MagVenture_B70_Litz',
                        'No5_Magstim_circular_90mm', 'No31_Magstim_70mm_Fig8', 'No6_MST_animal'
                        'No34_Neuronetics', 'No7_MST_human_circular', 'No36_MST_twin_coil_100', 'No9_H1'],
                    required=False) # TODO: will need to add these additional coil names into container
parser.add_argument('--low_mem',action='store_true', help='Run with less RAM (could take longer).'                    
parser.add_argument('--optimization_target',action='store_true',help='Whether to perform optimization (tDCS or TMS) targeting or not. '
                    'If using this feature please cite for TMS: '
                    'Weise, K., Numssen, O., Thielscher, A., Hartwigsen, G., & Knösche, T. R. (2020). A novel approach to localize cortical TMS effects. Neuroimage, 209, 116486. '
                    'or Gomez, L. J., Dannhauer, M., & Peterchev, A. V. (2020). Fast computational optimization of TMS coil placement for individualized electric field targeting. NeuroImage 2021; 228: 117696.
                    'Cite for tDCS: Saturnino, G. B., Siebner, H. R., Thielscher, A., & Madsen, K. H. (2019). Accessibility of cortical regions to focal TES: Dependence on spatial position, safety, and practical constraints. NeuroImage, 203, 116183.') # TODO: with logic to start this, "if not args.optimization_target"
parser.add_argument('--participant_label', help='The label(s) of the participant(s) that should be analyzed. The label '
                   'corresponds to sub-<participant_label> from the BIDS spec '
                   '(so it does not include "sub-"). If this parameter is not '
                   'provided all subjects should be analyzed. Multiple '
                   'participants can be specified with a space separated list.',
                   nargs="+")
parser.add_argument('--stimulation_type', choices=['TMS','tDCS'], required=True,
                    help='The stimulation type of the transcranial stimulation. Choices are "TMS" and "tDCS".')
parser.add_argument('-v', '--version', action='version',
                    version='SimNIBS version {}'.format(__version__))

args = parser.parse_args()

# unit/error checking submitted arguments
if args.stimulation_type == 'TMS':
    assert args.coil_name, print('If TMS is specified a coil_name must be specified following the argument "--coil_name". Choices are "Magstim_70mm_fig8" and "MagVenture_MC_B70".')
if not args.coil_name:
    args.coil_name = ''
if not args.optimization_target:
    args.optimization_target = False
if not args.low_mem:
    args.low_mem = False

subjects_to_analyze = []
# only for a subset of subjects
if args.participant_label:
    subjects_to_analyze = args.participant_label
# for all subjects
else:
    subject_dirs = glob(os.path.join(args.bids_dir, "sub-*"))
    subjects_to_analyze = [subject_dir.split("-")[-1] for subject_dir in subject_dirs]

# running participant level
if args.analysis_level == "participant":
    # pass to SimNIBS_simulation class for head reconstruction and simulation
    for subject_label in subjects_to_analyze:
        session_dirs = glob(os.path.join(args.bids_dir,'sub-'+str(subject_label),'ses-*'))
        if len(session_dirs) > 0:
            sessions_to_analyze = [session_dir.split("-")[-1] for session_dir in session_dirs]
            for session_label in sessions_to_analyze:
                SimNIBS_simulation(BIDS_folder=args.bids_folder,
                    coil_name=args.coil_name,
                    derivatives_folder=args.output_dir,
                    low_mem=args.low_mem
                    optimization=args.optimization_target,
                    session=session_label,
                    stimulation_type=args.stimulation_type,
                    subject=subject_label)    
        else:   
             SimNIBS_simulation(BIDS_folder=args.bids_folder,
                coil_name=args.coil_name,
                derivatives_folder=args.output_dir,
                low_mem=args.low_mem,
                optimization=args.optimization_target,
                stimulation_type=args.stimulation_type,
                subject=subject_label)                   
# running group level
elif args.analysis_level == "group":
    pass