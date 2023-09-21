import pandas as pd
import re
import RGB_Rule
import numpy as np


# function to read the protein fasta into a string for regular expression
def read_fasta_protein(file_name):
    with open(file_name) as file_open:
        file_content = file_open.read()
        file = file_content.split("\n")
        protein_string = "".join(file[1:])
        return protein_string


# function to read the CSV into list
# still working
def read_csv_peptide(file_name, i):
    columns = [0, 1, 2]
    file_content = pd.read_excel(file_name, header=None, skiprows=1, names=columns)
    file = [element for element in file_content[i] if str(element) != "nan"]
    return file


# function to create a dictioanry to have peptide list with no PTMs( for matching)
# The value would be the PTM amino acids and the position
def split_in_two(pep_list):
    pep_dict = {re.compile('[^a-zA-Z]').sub("", i): [re.findall(r"\w\[[\d]{3}\]", i),
                                                     [pos.start()for pos in re.finditer(r"\w\[[\d]{3}\]", i)]] for i in
                pep_list}
    for values in pep_dict.values():
        if len(values[1]) > 0:
            for i in range(0, len(values[1])):
                values[1][i] = values[1][i]-5*i
    return pep_dict


# function to match the peptide onto the protein seq
def match_on_peptide(pep_dict, protein_seq):
    position_list = [[pos.start(), pos.end(), pep_dict[i][1], pep_dict[i][0]] for i in pep_dict.keys() for pos in
                     re.finditer(i, protein_seq)]
    return position_list


def overlap_result_ptm(pep_position_list: list, protein_fasta: str) -> list:
    matrix_of_ptm = ["_" for i in range(0, len(protein_fasta))]
    for i in pep_position_list:
        for j in range(i[0], i[1]):
            if len(i[2]) > 0:
                for ptm in range(0, len(i[2])):
                    matrix_of_ptm[i[0] + i[2][ptm]] = i[3][ptm]
            else:
                continue
    return matrix_of_ptm


def overlap_result_all(pep1_position_list: list, pep2_position_list: list, pep3_position_list: list,
                       protein_fasta: str) -> list:
    matrix_of_0 = [0 for i in range(0, len(protein_fasta))]
    for i in pep1_position_list:
        if len(i[2]) > 0:
            for ptm in i[2]:
                matrix_of_0[i[0] + ptm] += 0.1
        for j in range(i[0], i[1]):
            matrix_of_0[j] += 1
    for i in pep2_position_list:
        if len(i[2]) > 0:
            for ptm in i[2]:
                matrix_of_0[i[0] + ptm] += 0.01
        for j in range(i[0], i[1]):
            matrix_of_0[j] += 10
    for i in pep3_position_list:
        if len(i[2]) > 0:
            for ptm in i[2]:
                matrix_of_0[i[0] + ptm] += 0.001
        for j in range(i[0], i[1]):
            matrix_of_0[j] += 100
    return matrix_of_0


def pos_list_to_matrix(pos_list, protein_list):
    matrix_of_space = ["_" for i in range(0, len(protein_list))]
    for i in pos_list:
        for j in range(i[0], i[1]):
            matrix_of_space[j] = "—"
        if len(i[2]) > 0:
            for ptm in range(0, len(i[2])):
                matrix_of_space[i[0] + i[2][ptm]] = i[3][ptm][0]
    return matrix_of_space


def ptm_set(pep1_dict, pep2_dict, pep3_dict):
    all_set = set()
    for i in pep1_dict.values():
        if i[0]:
            all_set.add(i[0][0])
    for i in pep2_dict.values():
        if i[0]:
            all_set.add(i[0][0])
    for i in pep3_dict.values():
        if i[0]:
            all_set.add(i[0][0])
    all_set_list = list(sorted(all_set))
    return all_set_list


def ptm_colors(ptm_set_in):
    ptm_color = {i: tuple(np.random.randint(0, 255, size=3)) for i in ptm_set_in}
    return ptm_color


def colors_align(pep_align, ptm_color):
    color_list = [ptm_color[i] if len(i) > 1 else (0, 0, 0) for i in pep_align]
    return color_list

def replace_non_zero_digits_with_one(input_str):
    input_str = str(input_str)
    result_str = ""
    for char in input_str:
        if char.isdigit() and char != '0':
            result_str += '1'
        else:
            result_str += char
    return result_str


def apply_color(matrix):
    color = [RGB_Rule.color[round(float(i), 0)] for i in matrix]
    return color


def apply_ptm(matrix):
    ptm = [RGB_Rule.ptm[round(float(i) % 1, 3)] for i in matrix]
    return ptm

def domains(domain_seq, protein_seq):
    position_list_domain = [[i.start(), i.end()] for i in re.finditer(domain_seq, protein_seq)]
    return position_list_domain

def domain_colors(domain_position, domain, protein_seq, domain_color):
    matrix_of_white = [0 for i in range(0, len(protein_seq))]
    for i in domain_position:
        for j in range(i[0], i[1]):
            matrix_of_white[j] = 1
    color_matrix = [domain_color[domain] if i == 1 else "white" for i in matrix_of_white]
    return color_matrix

def merge_domain(domain_list, pro_seq):
    matrix_of_white = ["white" for i in range(0, len(pro_seq))]
    for i in range(0, len(domain_list)):
        for color in range(0, len(domain_list[i])):
            if domain_list[i][color] != "white":
                matrix_of_white[color] = domain_list[i][color]
    return matrix_of_white


if __name__ == "__main__":
    file_name = "C:/Users/Orangeeee/Downloads/sequence.fasta"
    csv_name = "C:/Users/Orangeeee/Downloads/SCV-input-PSMs-081823.xlsx"
    protein_seq = read_fasta_protein(file_name)
    eda = "KGLAFTDVDVDSIKIAWESPQGQVSRYRVTYSSPEDGIRELFPAPDGEDDTAELQGLRPGSEYTVSVVALHDDMESQPLIGIQSTAIP"
    edb = "TDLSFVDITDSSIGLRWTPLNSSTIIGYRITVVAAGEGIPIFEDFVDSSVGYYTVTGLEPGIDYDISVITLINGGESAPTTLTQQTAVP"

    pep1_list = read_csv_peptide(csv_name, 0)
    pep2_list = read_csv_peptide(csv_name, 1)
    pep3_list = read_csv_peptide(csv_name, 2)

    pep1_dict = split_in_two(pep1_list)
    pep2_dict = split_in_two(pep2_list)
    pep3_dict = split_in_two(pep3_list)

    pep1_pos_list = match_on_peptide(pep1_dict, protein_seq)
    pep2_pos_list = match_on_peptide(pep2_dict, protein_seq)
    pep3_pos_list = match_on_peptide(pep3_dict, protein_seq)

    pep_ptm_1 = overlap_result_ptm(pep1_pos_list, protein_seq)
    pep_ptm_2 = overlap_result_ptm(pep2_pos_list, protein_seq)
    pep_ptm_3 = overlap_result_ptm(pep3_pos_list, protein_seq)

    ptm_sets = ptm_set(pep1_dict, pep2_dict, pep3_dict)
    ptm_color = ptm_colors(ptm_sets)

    colors_pep1 = colors_align(pep_ptm_1, ptm_color)
    colors_pep2 = colors_align(pep_ptm_2, ptm_color)
    colors_pep3 = colors_align(pep_ptm_3, ptm_color)

    eda_align = domains(eda, protein_seq)
    edb_align = domains(edb, protein_seq)

    eda_color = domain_colors(eda_align, "eda", protein_seq, RGB_Rule.domain_color)
    edb_color = domain_colors(edb_align, "edb", protein_seq, RGB_Rule.domain_color)

    merged_domain = [eda_color, edb_color]
    domain_col = merge_domain(merged_domain, protein_seq)

    ptm_sets = ptm_set(pep1_dict, pep2_dict, pep3_dict)
    overlap_matrix = overlap_result_all(pep1_pos_list, pep2_pos_list, pep3_pos_list, protein_seq)

    final_matrix = [replace_non_zero_digits_with_one(i) for i in overlap_matrix]
