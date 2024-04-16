"""Helper functions."""


def convert_list_to_float(data):
    """Returns a list containing the values of the given list as floats."""
    return [float(i) for i in data]


def convert_list_to_string(data):
    """Returns a list containing the values of the given list as strings."""
    return [str(i) for i in data]


def add_arrays(data):
    """Returns a list containing the sum of the given lists."""
    return [sum(float(i) for i in sublist) for sublist in zip(*data)]


def get_max_result(data, factor):
    """Returns the maximum value of the given list multiplied by the given factor."""
    return round(max(add_arrays(data)) * factor, 3)


def devide_arrays_in_percentage(list1, list2):
    """Returns the percentage of each element from list1 to list2."""
    return [round(float(i) / float(j) * 100, 2) for i, j in zip(list1, list2)]
