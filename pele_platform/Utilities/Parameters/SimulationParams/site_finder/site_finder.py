class SiteFinderParams(object):
    def generate_site_finder_params(self, args):
        self.site_finder = args.site_finder
        self.n_components = (
            args.n_components
            if args.n_components
            else self.simulation_params.get("n_components", 10)
        )
