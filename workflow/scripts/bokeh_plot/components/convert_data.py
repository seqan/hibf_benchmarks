"""Functions to convert the data for the plots."""


def prepare_time_data(time_data_list, keypair, time_configs):
    """Returns a dictionary containing the given data in the format for the time plot."""
    key, keyname = keypair
    data = time_data_list[time_data_list.loc[:, "KEY"] == key].copy()
    key_value = data.iloc[:, 1].str.split("=").str[1].astype(float).tolist()
    data["TOTAL_TIME"] = data.loc[:, time_configs["FORMAT"]].sum(axis=1)
    for col in data.columns:
        if "seconds" in col:
            data.insert(
                data.columns.get_loc(col) + 1,
                col.replace("in_seconds", "percentage"),
                (data[col] / data["TOTAL_TIME"] * 100).round(3),
            )
    data.loc[:, "KEY"] = keyname
    data.insert(2, "SUBKEY_VALUE", key_value)
    for name in time_configs["NAMES"]:
        data.insert(data.columns.get_loc(name) + 1, time_configs["NAMES"][name], data[name])
    data = data.sort_values(by=["SUBKEY_VALUE"])
    # print("Available data for time plot:")
    # print(data.columns)
    export = data.to_dict(orient="list")
    return export


def prepare_size_data(size_data, keypair, size_configs):
    """Returns a dictionary containing the given data in the format for the size plot."""
    key, keyname = keypair
    key_data = size_data[size_data.loc[:, "KEY"] == key].copy().sort_values(by=["SUBKEY", "LEVEL"])
    key_data = key_data.reset_index(drop=True)
    for i in range(4):
        key_data.loc[-(i + 1)] = i
    key_data.insert(5, "GB_SIZE", round((key_data["BIT_SIZE"].astype(float)) / (1000 * 1000 * 1000 * 8), 3))
    data = key_data.pivot_table(
        index="SUBKEY",
        columns="LEVEL",
        values=[col for col in key_data.columns if col not in ["KEY", "SUBKEY", "LEVEL", "SUBKEY_VALUE"]],
        aggfunc="first",
    )
    data = data.drop([0, 1, 2, 3])
    data.reset_index(inplace=True)
    data.fillna(0, inplace=True)
    data.columns = [f"LEVEL_{level}_{value}" if value != "SUBKEY" else value for (value, level) in data.columns]
    data["GB_TOTAL_SIZE"] = data.loc[:, size_configs["FORMAT"]].sum(axis=1)
    for col in data.columns:
        if "GB_SIZE" in col:
            data[col.replace("GB_SIZE", "GB_SIZE_percentage")] = (data[col] / data["GB_TOTAL_SIZE"] * 100).round(3)
    subkey_value_list = data["SUBKEY"].str.split("=").str[1].astype(float).tolist()
    data.insert(0, "SUBKEY_VALUE", subkey_value_list)
    data.insert(0, "KEY", keyname)
    for name in size_configs["NAMES"]:
        data.insert(data.columns.get_loc(name) + 1, size_configs["NAMES"][name], data[name])
    data = data.sort_values(by=["SUBKEY_VALUE"])
    # print("Available data for size plot:")
    # print(data.columns)
    export = data.to_dict(orient="list")
    return export
