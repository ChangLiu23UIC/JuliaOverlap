from flask import Flask, render_template
from Overlap_estimate import *

file_name = "C:/Users/Orangeeee/Downloads/sequence.fasta"
csv_name = "C:/Users/Orangeeee/Downloads/SCV-input-PSMs-081823.xlsx"
eda = "KGLAFTDVDVDSIKIAWESPQGQVSRYRVTYSSPEDGIRELFPAPDGEDDTAELQGLRPGSEYTVSVVALHDDMESQPLIGIQSTAIP"
edb = "TDLSFVDITDSSIGLRWTPLNSSTIIGYRITVVAAGEGIPIFEDFVDSSVGYYTVTGLEPGIDYDISVITLINGGESAPTTLTQQTAVP"
protein_seq = read_fasta_protein(file_name)

pep1_list = read_csv_peptide(csv_name, 0)
pep2_list = read_csv_peptide(csv_name, 1)
pep3_list = read_csv_peptide(csv_name, 2)

pep1_dict = split_in_two(pep1_list)
pep2_dict = split_in_two(pep2_list)
pep3_dict = split_in_two(pep3_list)

pep1_pos_list = match_on_peptide(pep1_dict, protein_seq)
pep2_pos_list = match_on_peptide(pep2_dict, protein_seq)
pep3_pos_list = match_on_peptide(pep3_dict, protein_seq)

pep_align1 = pos_list_to_matrix(pep1_pos_list, protein_seq)
pep_align2 = pos_list_to_matrix(pep2_pos_list, protein_seq)
pep_align3 = pos_list_to_matrix(pep3_pos_list, protein_seq)

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

overlap_matrix = overlap_result_all(pep1_pos_list, pep2_pos_list, pep3_pos_list, protein_seq)

final_matrix = [replace_non_zero_digits_with_one(i) for i in overlap_matrix]

color_list = apply_color(final_matrix)
ptm_list = apply_ptm(final_matrix)

app = Flask(__name__)

@app.route('/')
def index():
    return render_template('index.html', protein_string = protein_seq, char_color = color_list, size = len(protein_seq), ptm_list = ptm_list,
                           pep1 = pep_align1, pep2 = pep_align2, pep3 = pep_align3, color1 = colors_pep1, color2 = colors_pep2, color3 = colors_pep3, ptm_color = ptm_color,
                           merged = domain_col)


if __name__ == "__main__":
    app.run(debug=True)
