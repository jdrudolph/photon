import phos.data.read
import phos.data.network

def get_experiment(experiment_type):
    if experiment_type == 'silac':
        from phos.data.silac import experiment
        return experiment
    elif experiment_type == 'labelfree':
        from phos.data.labelfree import experiment
        return experiment
    else:
        raise NotImplementedError('experiment type {} unknown'.format(experiment_type))
