
from typing import Dict
import re
import os


def read_settings_from_mdp(fpath, overrides={}):
    pattern = r"([\w\-\_]+)[\t ]*=[\t ]*([\w\d\.\-\_\s]+)[\t ]*;*(.*)"
    with open(fpath) as mdp:
        for line in mdp.read().splitlines():
            cat = re.findall(pattern, line)
            if cat:
                name, parameter, comment = cat[0]
                parameter = parameter.strip()
                for o_name, o_parameter in overrides.items():
                    if name == o_name:
                        parameter = o_parameter
                yield name, parameter, comment


def build_settings_summary(self, arguments):
    settings_file = "settings.csv"

    Summary = {
        "force field": arguments.FF,
        "water model": arguments.solvent
    }

    mdp_prefixes = [
        "ions",
        "nvt",
        "npt",
        "md"
    ]

    header = ["Stage", "Parameter", "Value"]
    mdp_parameters = []

    for k in sorted(Summary.keys()):
        mdp_parameters.append(["*", k, Summary[k]])

    for mdp_prefix in mdp_prefixes:
        data = read_settings_from_mdp(
            self.to_wd(mdp_prefix + ".mdp"))
        for param, value in data:
            k = [mdp_prefix.upper(), param, value]
            mdp_parameters.append(k)

    mdp_parameters = self.compact_parameter_list(mdp_parameters)
    with open(self.to_wd(settings_file), 'w') as f:
        for parameter in [header] + mdp_parameters:
            f.write(",".join(parameter) + "\n")


def compact_parameter_list(parameters):
    compat = {}

    def make_key(parameter, value):
        return f"{parameter}:{value}"

    for stage, parameter, value in parameters:
        K = make_key(parameter, value)
        if K not in compat.keys():
            compat[K] = []
        compat[K].append(stage)

    output_parameters = []
    for stage, parameter, value in parameters:
        K = make_key(parameter, value)
        if K in compat.keys():
            stages = compat[K]
            message = "+".join(stages)
            output_parameters.append([message, parameter, value])
            del compat[K]

    return output_parameters


def load_mdp(self, arguments, mdpname):

    mdp_identifier = mdpname.split(".")[0]
    overrides = load_mdp_overrides_from_options(
        arguments, mdp_identifier
    )

    if os.path.isfile(mdpname):
        print(">Using user-defined %s" % mdpname)
        Source = mdpname
    else:
        print(">Writing built-in %s" % mdpname)
        Source = os.path.join(self.module_dir, "mdp", mdpname)

    Target = os.path.join(self.working_dir, mdpname)
    settings = read_settings_from_mdp(Source, overrides)
    write_mdp(Target, settings)
    # shutil.copy2(Source, Target)

def write_mdp(fpath, settings):
    with open(fpath, 'w') as f:
        for parameter, value, comment in settings:
            f.write(f"{parameter}\t\t=\t{value}\n")


def load_mdp_overrides_from_options(options, prefix) -> Dict[str, str]:
    overrides = {}
    m = "Override" + prefix
    for p in dir(options):
        if p.startswith(m):
            variable_name = p.replace(m, "")
            value = getattr(options, p)
            if value is not None:
                overrides[variable_name] = value

    return overrides


def add_option_override(parser, prefix, parameter):
    parser.add_argument(
        f"--{prefix.lower()}-{parameter}",
        dest=f"Override{prefix.upper()}{parameter}",
        type=int,
        help=f"Override the '{parameter}' parameter"
        f"inside the .mdp files for the {prefix} step."
    )
