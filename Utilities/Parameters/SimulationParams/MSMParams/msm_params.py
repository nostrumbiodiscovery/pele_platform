
class MSMParams(object):


    def __init__(self, args):
        self.lagtime = args.lagtime
        self.msm_clust = args.msm_clust
        self.sasa_max = None
        self.sasa_min = None
        self.clusters = args.clust = args.clust if not args.test else 2
