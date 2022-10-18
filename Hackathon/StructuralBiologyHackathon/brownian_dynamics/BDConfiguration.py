class BDConfiguration(object):
    def __init__(self, particles, borders):
        self.particles = particles
        self.borders = borders  # In pixels

    def get_particles(self):
        return self.particles

    def get_borders(self):
        return self.borders
