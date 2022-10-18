class Particle(object):
    """
    Base class for Particle Objects in 2D space
    """

    diffusion_coef = 0

    def __init__(self, x, y):
        """

        :param x: initial x coordinate value
        :param y: initial y coordinate value
        """
        self.x_cord = x
        self.y_cord = y
        self.x_tmp = x
        self.y_tmp = y
        self.is_phos = False

    def get_location(self):
        """

        :return: x,y coordinates of particle
        """
        return self.x_cord, self.y_cord

    def set_location(self):
        """

        :param x: new x_cord location
        :param y: new y_cord location
        :return:
        """
        self.x_cord = self.x_tmp
        self.y_cord = self.y_tmp

    def set_tmp_location(self,x,y):
        self.x_tmp = x
        self.y_tmp = y

    def get_phos(self):
        """

        :return: True - if particle is phosphorated, False - if not phosphorated
        """
        return self.is_phos

    def set_phos(self, phos_status):
        """
        Change phosphorylation status of particle
        :param phos_status: boolean value indicating if particles is
        phosphorylated
        :return:
        """
        self.is_phos = phos_status

    def get_particle_type(self):
        return self.__class__.__name__


    def get_diffusion_coef(self):
        return self.diffusion_coef



class MHC(Particle):
    diffusion_coef = 1e-10

    def get_particle_type(self):
        return self.__class__.__name__



class TCR(Particle):
    diffusion_coef = 1e-10

    def get_particle_type(self):
        return self.__class__.__name__


class CD45(Particle):
    diffusion_coef = 3.8e-10

    def get_particle_type(self):
        return self.__class__.__name__

