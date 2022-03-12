 
def ask_simulation_prefixes(simulation_prefixes):
    print("File prefixes found:")
    for i, prefix in enumerate(simulation_prefixes):
        print(f"{i + 1}:\t" + prefix)

    print("Select all prefixes? (empty input)")
    print("Or input comma separated numbers and " +
          "dash separated intervals to select prefixes.")

    return input(">")


def select_simulation_prefixes(simulation_prefixes, input_string):
    range_descriptors = [input_string]
    if "," in input_string:
        range_descriptors = [
            v.strip()
            for v in input_string.split(",")
        ]

    OutputPrefixes = []

    if range_descriptors in [[], ["all"]]:
        return simulation_prefixes

    for v in range_descriptors:
        try:
            if "-" in v:
                limits = v.split("-")
                assert len(limits) == 2
                F, T = [int(k) for k in limits]
                V = list(range(F, T + 1))
            else:
                V = [int(v)]

        except (ValueError, AssertionError):
            raise Exception("Invalid input.")

        for prefix_idx in V:
            OutputPrefixes.append(simulation_prefixes[prefix_idx - 1])

    for prefix in OutputPrefixes:
        print('\t' + prefix)

    return OutputPrefixes
