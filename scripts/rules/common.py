"""Generate parameters for the snakemake rules."""


def get_params(constrained):
    """Generate parameters for the snakemake rules."""
    # pylint: disable=undefined-variable
    result = []
    result += [f"alpha={str(param).replace('.', '_')}" for param in config["PARAMS"]["ALPHA"]]
    result += [
        f"{param1}={param2}"
        for param1 in config["PARAMS"]["MODE"] or [config["DEFAULT_PARAMS"]["MODE"]]
        for param2 in config["PARAMS"]["T_MAX"] or [config["DEFAULT_PARAMS"]["T_MAX"]]
    ]
    if constrained:
        result += [f"hash={param}" for param in config["PARAMS"]["NUM_HASHES"] if param < 6]
    else:
        result += [f"hash={param}" for param in config["PARAMS"]["NUM_HASHES"]]
    result += [f"kmer={param}" for param in config["PARAMS"]["KMER_SIZE"]]
    result += [f"relaxed-fpr={str(param).replace('.', '_')}" for param in config["PARAMS"]["RELAXED_FPR"]]
    return result
