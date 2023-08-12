from collections import defaultdict
from copy import deepcopy
from argparse import ArgumentParser


def main_parser(
    ori_to_mut: defaultdict,
    input_file: str,
    output_file: str,
    is_single_muts: bool = False,
    ) -> None:
 
    with open(input_file) as fh:
        all_lines = sorted(fh.readlines(), key=lambda x: x.split(",")[1])
        for line in all_lines:
            if not line.startswith(",name"):
                _, id, seq, _, score = line.strip().split(",")
                id = id.split(".pdb")
                if id[1] == "":
                    ori_to_mut[id[0]].append(["ori", seq, score])
                else:
                    if is_single_muts:
                        pos = "".join(id[1][2:-1])
                    answer = [seq, score, pos] if is_single_muts else [seq, score]
                    ori_to_mut[id[0]].append(answer)

    # dataset validation
    validated_dataset = deepcopy(ori_to_mut)
    for ori_id, protein_info in ori_to_mut.items():
        # no mutation proteins validation
        if len(protein_info) == 1:
            del validated_dataset[ori_id]
            print("no mutation proteins", ori_id)
            continue
        
        # several native proteins with different stability score
        # flag ori must be only in the first position
        if "ori" in protein_info[1] and protein_info[0][2] != protein_info[1][2]:
            del validated_dataset[ori_id]
            print("several native proteins with different stability score", ori_id)
            continue
        # no native protein validation
        if "ori" not in protein_info[0]:
            del validated_dataset[ori_id]
            print("no native protein", ori_id)
            continue
        
        # check that native protein and all mutation proteins have stability score
        if protein_info[0][2] == "":
            del validated_dataset[ori_id]
            print("native protein without stability score", ori_id)
            continue
        else:
            for mut_protein_info in protein_info[1:]:
                if mut_protein_info[1] == "":
                    validated_dataset[ori_id].remove(mut_protein_info)
                    print("mutation proteins without stability score", ori_id, mut_protein_info[0])
                    continue

    # lets find multiple mutation positions
    if not is_single_muts:
        for ori_id, protein_info in validated_dataset.items():
            native_seq = protein_info[0][1]
            for mut_protein_info in protein_info[1:]:
                mut_seq = mut_protein_info[0]
                mutation_positions = []
                if len(native_seq) != len(mut_seq):
                    print("native seq length and mutatoin seq length are not equal", ori_id)
                for idx in range(len(native_seq)):
                    if native_seq[idx] != mut_seq[idx]:
                        mutation_positions.append(idx)
                index_mut_item = validated_dataset[ori_id].index(mut_protein_info)
                validated_dataset[ori_id][index_mut_item].append(mutation_positions)

    # write result
    with open(output_file, "a+") as fw:
        for key, items in validated_dataset.items():
            for item in items:
                _, ori_seq, ori_score = items[0]
                if "ori" not in item:
                    mut_seq, mut_score, mut_positions = item
                    if isinstance(mut_positions, list):
                        mut_positions = "_".join([str(x) for x in mut_positions])
                    fw.write(f"{ori_seq},{mut_seq},{float(ori_score)-float(mut_score)},{mut_positions}\n")


if __name__ == "__main__":
    ori_to_mut = defaultdict(list)

    parser = ArgumentParser(description='Input files parser')

    parser.add_argument('--input', '-i', type=str, help='Path to input file', default="./single_muts_train.csv")
    parser.add_argument('--output', '-o', type=str, help='Path to output file', default="./train_dataset_single.csv")
    parser.add_argument('--is_single_muts', '-s', type=bool, help='Is single mutations', default=True)
    
    args = parser.parse_args()

    is_single_muts = args.is_single_muts
    input_file = args.input
    output_file = args.output
    # input_file = "./single_muts_train.csv" if is_single_muts else "./multiple_muts_train.csv"
    # output_file = "./train_dataset_single.csv" if is_single_muts else "./train_dataset_multiple.csv"

    if is_single_muts:
        # manually add native proteins which are not in file
        # find there stability scire and sequence in other file
        ori_to_mut["HEEH_rd2_0779"].append(("ori", "TLDEARELVERAKKEGTGVDVNGQRFEDWREAERWVREQEKNK", 1.08))
        ori_to_mut["HHH_rd2_0134"].append(("ori", "SKDEAQREAERAIRSGNKEEARRILEEAGYSPEQAERIIRKLG", 1.43))
        ori_to_mut["HHH_rd3_0138"].append(("ori", "ERRKIEEIAKKLYQSGNPEAARRFLRKAGISEEEIERILQKAG", 1.92))

    main_parser(ori_to_mut, input_file, output_file, is_single_muts)
