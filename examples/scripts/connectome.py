from synapseaware.isthmus import topological_thinning
from synapseaware.teaser import teaser
from synapseaware.connectome import wiring




prefix = 'Fib25'
label = 1



topological_thinning.TopologicalThinning(prefix, label)
teaser.TEASER(prefix, label)
wiring.GenerateSkeleton(prefix, label)
wiring.RefineSkeleton(prefix, label)
