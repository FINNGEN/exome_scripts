#!/usr/bin/env python3

import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Check overlap and mapping between two sample lists using an ID mapping file."
    )
    parser.add_argument(
        "list1",
        type=str,
        help="First sample list file (one sample ID per line)"
    )
    parser.add_argument(
        "list2",
        type=str,
        help="Second sample list file (one sample ID per line)"
    )
    parser.add_argument(
        "--mapping",
        type=str,
        default="/mnt/disks/data/samples/exclusions/finngen_R12_duplicate_list.txt",
        help="Tab-delimited mapping file (default: %(default)s)"
    )
    args = parser.parse_args()

    # Load sample lists
    set1 = set(line.strip() for line in open(args.list1) if line.strip())
    set2 = set(line.strip() for line in open(args.list2) if line.strip())

    # Load mapping
    mapping = [set(line.strip().split()) for line in open(args.mapping) if line.strip()]
    id2grp = {i: g for g in mapping for i in g}

    # Find matches
    direct = set1 & set2
    indirect = set(
        s for s in set1 - direct
        if any(x in set2 for x in id2grp.get(s, {s}) - {s})
    )
    no_match = set1 - direct - indirect

    print(f"Direct matches: {len(direct)}")
    print(f"Indirect matches via mapping: {len(indirect)}")
    print(f"No match: {len(no_match)}")

if __name__ == "__main__":
    main()
