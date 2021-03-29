import yaml


def read_yaml(f):
    with open(f) as temp:
        output = yaml.load(temp.read())
    return output


def logparabola(E, theta):
    N0, Gamma = theta
    dNdE = 10 ** N0 * (E / 0.1) ** Gamma
    return dNdE
