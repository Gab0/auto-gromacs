from typing import List


def ask_simulation_prefixes(simulation_prefixes) -> str:
    print("File prefixes found:")
    for i, prefix in enumerate(simulation_prefixes):
        print(f"{i + 1}:\t" + prefix)

    print("Select all prefixes? (empty input)")
    print("Or input comma separated numbers and " +
          "dash separated intervals to select prefixes.")

    return input(">")


def process_range_descriptors(input_string: str, total_size: int) -> List[int]:
    sector_splitter = ";"
    range_splitter = "-"

    range_descriptors = [input_string]
    if sector_splitter in input_string:
        range_descriptors = [
            v.strip()
            for v in input_string.split(sector_splitter)
        ]

    if range_descriptors in [[], ["all"]]:
        return list(range(total_size))

    selected_indexes = []
    for v in range_descriptors:
        try:
            if range_splitter in v:
                limits = v.split(range_splitter)
                assert len(limits) == 2
                F, T = [int(k) for k in limits]
                V = list(range(F, T + 1))
            else:
                V = [int(v)]
            selected_indexes += [idx - 1 for idx in V]
        except (ValueError, AssertionError) as error:
            raise Exception("Invalid input.") from error

    return selected_indexes


def select_simulation_prefixes(simulation_prefixes, input_string) -> List[str]:

    selected_indexes = process_range_descriptors(input_string, len(simulation_prefixes))

    output_prefixes = [simulation_prefixes[idx] for idx in selected_indexes]
    for prefix in output_prefixes:
        print('\t' + prefix)

    return output_prefixes
