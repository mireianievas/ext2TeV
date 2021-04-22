import yaml


def read_yaml(f):
    with open(f) as temp:
        output = yaml.load(temp.read())
    return output
