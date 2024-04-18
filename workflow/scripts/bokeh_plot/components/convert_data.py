"""Functions to convert the data for the plots."""

from components.helpers import add_arrays, convert_list_to_float, devide_arrays_in_percentage


def convert_time_data(data, key, time_format):
    """Returns a dictionary containing the given data in the format for the time plot."""
    result = {}
    result["SUBKEY"] = [f"{key} = {value}" for value in data[0]]
    result["value"] = data[0]
    result["all_times"] = add_arrays(data[1:])
    for i, element in enumerate(data[1:]):
        result[time_format[i + 1]] = convert_list_to_float(element)
        result[f"{time_format[i + 1]}_percentage"] = devide_arrays_in_percentage(
            result[time_format[i + 1]], result["all_times"]
        )
    return result


def convert_size_data(data, key, size_format):
    """Returns a dictionary containing the given data in the format for the size plot."""
    result = {}
    result[size_format[0]] = [f"{key} = {value}" for value in data[0]]
    result["value"] = data[0]
    result["sizes"] = [sum(float(i) for i in sublist) for sublist in zip(*data[1:5])]
    for i, element in enumerate(data[1:5]):
        result[size_format[i + 1]] = convert_list_to_float(element)
        result[f"{size_format[i+1]}_percentage"] = devide_arrays_in_percentage(
            result[size_format[i + 1]], result["sizes"]
        )
    for i, element in enumerate(data[5:]):
        result[f"{size_format[i+1]}_avg_load_factor"] = element
    return result
