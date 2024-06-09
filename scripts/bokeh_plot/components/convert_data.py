def prepare_time_data(time_data_list, keypair, timepairs):
    """Returns a dictionary containing the given data in the format for the time plot."""
    key, keyname = keypair
    data = time_data_list[time_data_list.loc[:, "KEY"] == key]
    key_value = data.iloc[:, 1].str.split("=").str[1].astype(float).tolist()
    for col in data.columns:
        if "seconds" in col:
            data.insert(data.columns.get_loc(col)+1, col.replace("in_seconds", "percentage"), (data[col] / data["wall_clock_time_in_seconds"] * 100).round(3))
    data.loc[:, "KEY"] = keyname
    data.insert(2, "SUBKEY_VALUE", key_value)
    for name in timepairs:
        data.insert(data.columns.get_loc(name)+1, timepairs[name], data[name])
    data = data.sort_values(by=["SUBKEY_VALUE"])
    # print("Possible data for time plot:")
    # print(data.columns)
    export = data.to_dict(orient="list")
    return export


def prepare_size_data(size_data, keypair, sizepairs):
    """Returns a dictionary containing the given data in the format for the size plot."""
    key, keyname = keypair
    key_data = size_data[size_data.loc[:, "KEY"] == key].sort_values(by=["SUBKEY", "LEVEL"])
    key_data = key_data.reset_index(drop=True)
    for i in range(4):
        key_data.loc[-(i+1)] = i
    key_data.insert(5, "GB_SIZE", round((key_data["BIT_SIZE"].astype(float)) / (1000*1000*1000*8), 3))
    data = key_data.pivot_table(
        index="SUBKEY",
        columns="LEVEL",
        values=[col for col in key_data.columns if col not in ["KEY", "SUBKEY", "LEVEL", "SUBKEY_VALUE"]],
        aggfunc="first"
    )
    data = data.drop([0, 1, 2, 3])
    data.reset_index(inplace=True)
    data.fillna(0, inplace=True)
    data.columns = [f"LEVEL_{level}_{value}" if value != "SUBKEY" else value for (value, level) in data.columns]  
    data["GB_TOTAL_SIZE"] = [data["LEVEL_0_GB_SIZE"][i] + data["LEVEL_1_GB_SIZE"][i] + data["LEVEL_2_GB_SIZE"][i] + data["LEVEL_3_GB_SIZE"][i] for i in range(len(data))]
    data["BIT_TOTAL_SIZE"] = [data["LEVEL_0_BIT_SIZE"][i] + data["LEVEL_1_BIT_SIZE"][i] + data["LEVEL_2_BIT_SIZE"][i] + data["LEVEL_3_BIT_SIZE"][i] for i in range(len(data))]
    for col in data.columns:
        if "GB_SIZE" in col:
            data[col.replace("GB_SIZE", "GB_SIZE_percentage")] = (data[col] / data["GB_TOTAL_SIZE"] * 100).round(3)
    subkey_value_list = data["SUBKEY"].str.split("=").str[1].astype(float).tolist()
    data.insert(0, "SUBKEY_VALUE", subkey_value_list)
    data.insert(0, "KEY", keyname)
    for name in sizepairs:
        data.insert(data.columns.get_loc(name)+1, sizepairs[name], data[name])
    data = data.sort_values(by=["SUBKEY_VALUE"])
    # print("Possible data for size plot:")
    # print(data.columns)
    export = data.to_dict(orient="list")
    return export