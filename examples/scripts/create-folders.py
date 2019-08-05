import os
import sys


prefix = sys.argv[1]


if not os.path.exists('meta'): os.mkdir('meta')
if not os.path.exists('running_times'): os.mkdir('running_times')
if not os.path.exists('running_times/refinement'): os.mkdir('running_times/refinement')
if not os.path.exists('running_times/skeletons'): os.mkdir('running_times/skeletons')
if not os.path.exists('width-errors'): os.mkdir('width-errors')
if not os.path.exists('baselines'): os.mkdir('baselines')
if not os.path.exists('volumetric_somae'): os.mkdir('volumetric_somae')


output_directories = ['baselines/teasers', 'baselines/topological-thinnings', 'connectomes', 'distances', 'segmentations', 'skeletons', 'somae', 'surfaces', 'synapses', 'volumetric_somae/segmentations', 'widths']


for output_directory in output_directories:
    if not os.path.exists(output_directory): os.mkdir(output_directory)
    if not os.path.exists('{}/{}'.format(output_directory, prefix)): os.mkdir('{}/{}'.format(output_directory, prefix))


