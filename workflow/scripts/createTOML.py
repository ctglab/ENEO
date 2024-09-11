#!/usr/bin/python3

# Helper function to generate the TOML file required for the 
# annotation with vcfanno.

import argparse
from pip._vendor import tomli
import toml
import yaml

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-y", "--yaml", help="yaml file with the configuration", required=True)
    parser.add_argument("-t", "--template", help="TOML template file", required=True)
    parser.add_argument("-o", "--output", help="TOML output file", required=True)
    return parser.parse_args()


def main(yaml_file, template_file, output_file):
    # read the yaml file
    with open(yaml_file, 'r') as f:
        conf_main = yaml.safe_load(f)
    # load the toml template
    with open(template_file, 'rb') as f:
        conf_template = tomli.load(f)
    # update the template with the configuration
    for field in conf_template['annotation']:
        if "gnomAD" in field['names'][0]:
            field['file'] = conf_main['resources']['gnomad']
        elif "rs_ids" in field['names'][0]:
            field['file'] = conf_main['resources']['dbsnps']
        elif "REDI" in field['names'][0]:
            field['file'] = conf_main['resources']['REDI']
        elif "indel" in field['names'][0]:
            field['file'] = conf_main['resources']['indel']
        elif "cosmic" in field['names'][0]:
            field['file'] = conf_main['resources']['cosmic']
        elif "Unmet" in field['names'][0]:
            field['file'] = conf_main['resources']['unmet_bed']
        else:
            raise ValueError("Unknown field: {}".format(field['names'][0]))
    # write the output file
    with open(output_file, "w") as f:
        toml.dump(conf_template, f)


if __name__ == "__main__":
    args = parse_args()
    main(args.yaml, args.template, args.output)
