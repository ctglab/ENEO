#!/usr/bin/env python3

"""
Generate a TOML file for vcfanno by merging a YAML configuration
with a TOML template. All resource paths are taken from the YAML.
"""

import argparse
import toml
import yaml


def parse_args():
    parser = argparse.ArgumentParser(description="Generate TOML annotation file for vcfanno.")
    parser.add_argument("-y", "--yaml", required=True, help="YAML configuration file")
    parser.add_argument("-t", "--template", required=True, help="TOML template file")
    parser.add_argument("-o", "--output", required=True, help="Output TOML file")
    return parser.parse_args()


def load_yaml(path):
    with open(path, "r") as f:
        return yaml.safe_load(f)


def load_toml(path):
    with open(path, "r") as f:
        return toml.load(f)


def update_template(template, config):
    """Replace the input template paths using rules based on field name."""
    resources = config.get("resources", {})
    mapping_rules = [
        ("gnomAD", "gnomad"),
        ("rs_ids", "dbsnps"),
        ("REDI", "REDI"),
        ("indel", "indel"),
        ("Unmet", "unmet_bed"),
    ]
    for field in template.get("annotation", []):
        name = field["names"][0]
        matched = False
        for key_substring, resource_key in mapping_rules:
            if key_substring in name:
                field["file"] = resources[resource_key]
                matched = True
                break
        if not matched:
            raise ValueError(f"Unknown field in template: {name}")

    return template


def main():
    args = parse_args()
    config = load_yaml(args.yaml)
    template = load_toml(args.template)
    updated = update_template(template, config)
    with open(args.output, "w") as f:
        toml.dump(updated, f)


if __name__ == "__main__":
    main()
