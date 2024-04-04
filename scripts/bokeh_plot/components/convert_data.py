from components.helpers import add_arrays, convert_list_to_float, devide_arrays_in_percentage

def convert_time_data(data, key, TIME_FORMAT):
    """Returns a dictionary containing the given data in the format for the time plot."""
    export = {}
    export["SUBKEY"] = [f"{key} = {value}" for value in data[0]]
    export["value"] = data[0]
    export["all_times"] = add_arrays(data[1:])
    for i, element in enumerate(data[1:]):
        export[TIME_FORMAT[i + 1]] = convert_list_to_float(element)
        export[f"{TIME_FORMAT[i + 1]}_percentage"] = devide_arrays_in_percentage(
            export[TIME_FORMAT[i + 1]], export["all_times"]
        )
    return export


def convert_size_data(data, key, SIZE_FORMAT):
    """Returns a dictionary containing the given data in the format for the size plot."""
    export = {}
    export[SIZE_FORMAT[0]] = [f"{key} = {value}" for value in data[0]]
    export["value"] = data[0]
    export["sizes"] = [sum(float(i) for i in sublist) for sublist in zip(*data[1:5])]
    for i, element in enumerate(data[1:5]):
        export[SIZE_FORMAT[i + 1]] = convert_list_to_float(element)
        export[f"{SIZE_FORMAT[i+1]}_percentage"] = devide_arrays_in_percentage(
            export[SIZE_FORMAT[i + 1]], export["sizes"]
        )
    for i, element in enumerate(data[5:]):
        export[f"{SIZE_FORMAT[i+1]}_avg_load_factor"] = element
    return export
