"""Helper functions."""


def convert_dic_to_list(data):
    """Returns a list containing the values of the given list as strings."""
    return [str(i) for i in data]


def add_arrays(data):
    """Returns a list containing the sum of the given lists."""
    return [sum(float(i) for i in sublist) for sublist in zip(*data)]


def get_max_result(data, factor):
    """Returns the maximum value of the given list multiplied by the given factor."""
    return round(max(data) * factor, 3)
