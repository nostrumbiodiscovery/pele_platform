class GenerativeModel:

    def __init__(self, core, resname, iterations, only_grow=True, only_rank=False):
        self.core = core
        self.iterations = iterations
        self.only_grow = only_grow
        self.only_rank = only_rank
        self.resname = resname

    def run(self):
        from growai import grow
        return grow.main(self.core, self.resname, self.only_grow, self.only_rank, self.iterations)