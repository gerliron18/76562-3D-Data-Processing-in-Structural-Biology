from Forces import force_calculation

ALL_FORCE_FUNCTIONS = ["location_interaction"]


class MC_coefficient(object):
    """
    Class representing the MC coefficients
    """

    def __init__(self, name, min, max, step):
        """
        :param name: Name of the coefficient
        :param min: The minimum value for the coefficient
        :param max: The maximum value for the coefficient
        :param step: The step size for the MC
        """
        self.name = name
        self.min = min
        self.max = max
        self.step = step

    def get_name(self):
        return self.name


def get_force_function(func_name):
    """
    For the requested func_name calculates the required coefficients, and returns a tuple of the function and a list of
    the parameters.
    :param func_name: Name of the function, should be in the ALL_FORCE_FUNCTIONS list
    :return: Tuple of the function and a list of the parameters.
    """
    if func_name not in ALL_FORCE_FUNCTIONS:
        raise NotImplementedError
    if func_name == "location_interaction":
        a1 = MC_coefficient("interaction", -0.1, 0.1, 0.01)
        a2 = MC_coefficient("location_a", -0.1, 0.1 ,0.01)
        b = MC_coefficient("location_b", -10, 10, 0.1)
        return force_calculation.location_interaction, [a1, a2, b]
