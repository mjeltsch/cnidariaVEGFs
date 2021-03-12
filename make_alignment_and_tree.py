#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
import os, sqlite3, ntpath, types, shutil, glob
from Bio import SeqIO
from Bio import Entrez
from pathlib import Path
from phylolib import execute_subprocess
from shutil import copyfile
from ete3 import Tree, TreeStyle, AttrFace, TextFace, NodeStyle, SequenceFace, ImgFace, SVGFace, faces, add_face_to_node, PhyloTree

# A Psi blast  limited to the taxon Viridae was run against the starting sequence AAD03735.1 until
# no new sequences were found above the 0.005 threshold. yielding 84 sequences.
# Fasta names were ajusted to include the host species (resulting in file seqdump_modified_names.fasta).
#
# PCMA:
# wget http://prodata.swmed.edu/download/pub/PCMA/pcma.tar.gz
# tar -xvzf pcma.tar.gz
# cd pcma
# make
# sudo cp pcma /bin
#
# Fix dialign-t and clustalw:
# sudo ln -s /usr/bin/dialign-tx /bin/dialign-t\
# sudo ln -s /usr/bin/clustalw /bin/clustalw2
#
# Something like the following might needed to fix multiprocessor phyml if you get an error like this:
# "error while loading shared libraries: libmpi.so.12: cannot open shared object file: No such file or directory":
# cd /usr/lib/x86_64-linux-gnu/
# sudo ln -s  openmpi/lib/libmpi.so.20.10.1 libmpi.so.12
#

def do_alignment(infile, outfile):
    execute_subprocess(
        'Generating multiple sequence alignment with the following command:',
        't_coffee {0} -outfile {1} -output=fasta_aln -mode mcoffee'.format(infile, outfile))

def trim_difficult_streches(file):
    execute_subprocess(
        '',
        # -b4 = minimum length of a block
        # -b5 = allowed gap positions (none, half, all = n,h,a)
        'Gblocks {0} -t=p -b4=10'.format(file))

def remove_extra_text_after_fasta_description(infile, outfile):
    execute_subprocess(
        'Removing everything from fasta description after the first blank space:',
        'sed \'s/\s.*$//\' {0} > {1}'.format(infile, outfile))

def encode_fasta_descriptions(fasta_file_name, list_file_name, encoded_aligned_fasta_file):
    # Make code name list
    execute_subprocess(
        'Converting fasta descriptions part 1 (creating code list) with t_coffee using the following command:',
        't_coffee -other_pg seq_reformat -in {0} -output code_name > {1}'.format(fasta_file_name, list_file_name))
    # Encode fasta file
    execute_subprocess(
        'Converting fasta descriptions part 2 (replacing fasta descriptions with codes) with t_cofeee using the following command:',
        't_coffee -other_pg seq_reformat -code {0} -in {1} > {2}'.format(list_file_name, fasta_file_name, encoded_aligned_fasta_file))

def convert_into_phylip(infile, outfile):
    execute_subprocess(
        'Convert into phylip using the following command:',
        't_coffee -other_pg seq_reformat -in {0} -output phylip_aln > {1}'.format(infile, outfile))

def make_tree(file):
    # This accepts only a phylip file
    # Detect whether parallel bootstrapping should be performed
    # b (bootstrap): -1 is aLRT statistics, for final analysis use: -b 100 or 1000 (on a fast computer)
    mpirun_path = shutil.which('mpirun')
    phymlmpi_path = shutil.which('phyml-mpi')
    #
    # -x LD_LIBRARY_PATH /usr/lib/x86_64-linux-gnu/libmpi.so.20.10.1
    #
    if mpirun_path != '' and phymlmpi_path != '':
        phylo_command = 'mpirun -n 4 phyml-mpi -i {0} -d aa -b -1'.format(file)
    else:
        phylo_command = 'phyml -i {0} -d aa -b -1'.format(file)

    execute_subprocess(
        'Make tree with the following command:',
        phylo_command)

    # phyml adds or doesn't add the .txt extension to the output file (depending on the version) and we need to check for this!
    phyml_output_file = '{0}_phyml_tree'.format(file)
    if os.path.isfile(phyml_output_file):
        os.rename(phyml_output_file, '{0}.txt'.format(phyml_output_file))
    return '{0}.txt'.format(phyml_output_file)

def decode_fasta_descriptions(list_file, infile, outfile):
    execute_subprocess(
        'Decoding tree file file into human-readable format using the following command:',
        't_coffee -other_pg seq_reformat -decode {0} -in {1} > {2}'.format(list_file, infile, outfile))

def vegfe_tree_layout(node):
    nameFace1 = faces.TextFace(" ", fsize=12)
    nameFace2 = faces.AttrFace("name", fsize=12, fgcolor="Black")

    if node.is_leaf():
        # Add an static face that handles the node name
        faces.add_face_to_node(nameFace1, node, column=0)
        faces.add_face_to_node(nameFace2, node, column=1)

        # Sets the style of leaf nodes
        node.img_style["size"] = 10
        node.img_style["shape"] = "circle"

        #If node is an internal node
    else:
        # Sets the style of internal nodes
        node.img_style["size"] = 8
        node.img_style["shape"] = "circle"
        node.img_style["fgcolor"] = "#000000"

    # If leaf is "Hsa" (homo sapiens), highlight it using a
    # different background.
    if node.is_leaf() and node.name.startswith("Hsa"):
        node.img_style["bgcolor"] = "#9db0cf"

def run():
    do_calculations = True
    sequences = 'cnidaria_VEGFs.fasta'
    aligned_sequences = 'cnidaria_VEGFs_aligned.fasta'
    list_file_name = "cnidaria_VEGFs.lst"
    encoded_aligned_sequences = 'cnidaria_VEGFs_aligned_encoded.fasta'
    phylip_file = 'cnidaria_VEGFs.phylip'
    phyml_output_file = 'phyml_result'
    final_tree = 'cnidaria_VEGFs.tree'
    svg_file = 'cnidaria_VEGFs_tree.svg'

    if do_calculations == True:
        do_alignment(sequences, aligned_sequences)
        encode_fasta_descriptions(aligned_sequences, list_file_name, encoded_aligned_sequences)
        #remove_extra_text_after_fasta_description(encoded_aligned_fasta_file, encoded_aligned_trimmed_fasta_file)
        convert_into_phylip(encoded_aligned_sequences, phylip_file)
        phyml_output_file = make_tree(phylip_file)
        decode_fasta_descriptions(list_file_name, phyml_output_file, final_tree)

    t = Tree(final_tree)

    # Coloring regions
    # White = NodeStyle()
    # White["bgcolor"] = "White"
    # #n1 = t.get_common_ancestor("SSV_Lagothrix_lagotricha_P01128.1", "SSV_Lagothrix_lagotricha_YP_001165471.3")
    # #n1.set_style(White)
    # MCV
    # Green = NodeStyle()
    # Green["bgcolor"] = "#dce6dc"
    # n2 = t.get_common_ancestor("MCV_Squalius_cephalus_QCQ67828.1", "MCV_Lates_calcarifer_YP_009163772.1")
    # n2.set_style(Green)
    # # PPV
    # Gray = NodeStyle()
    # Gray["bgcolor"] = "Gainsboro"
    # # PPV main clade
    # n11 = t.get_common_ancestor("RDPPV_Cervus_elaphus_YP_009112871.1", "PCPV_Rangifer_tarandus_ADC53897.1")
    # n11.set_style(Gray)
    # # PPV separate clade
    # n12 = t.get_common_ancestor("PCPV_Bos_taurus_ACJ06526.1", "PCPV_Bos_taurus_ACJ06525.1")
    # n12.set_style(Gray)
    # # ORFV
    # Beige = NodeStyle()
    # Beige["bgcolor"] = "#efebe4"
    # # Orf region 1
    # n3 = t.get_common_ancestor("ORFV_Capra_hircus_AWN09372.1", "ORFV_Ovis_aries_P52585.1")
    # n3.set_style(Beige)
    # # Orf region 2
    # n13 = t.get_common_ancestor("ORFV_Rupicapra_rupicapra_QJD15044.1", "ORFV_Capricornis_crispus_BAI50009.1")
    # n13.set_style(Beige)
    # # # Mammalian VEGFB
    # # # THERE IS A BUG HERE!!! We have to miscolor this branch on purpose in white to get the overall background white!
    # Yellow = NodeStyle()
    # Yellow["bgcolor"] = "#ffff99"
    # n5 = t.get_common_ancestor("VEGFB_Homo_sapiens_AAL79000.1", "VEGFB_Cervus_elaphus_OWK17446.1")
    # n5.set_style(Yellow)
    # # n5.swap_children()
    # # # ???
    # # n9 = t.get_common_ancestor("VEGFB_Capra_hircus_XP_017898275.1", "VEGFB_Homo_sapiens_AAL79000.1")
    # # n9.set_style(White)
    # # Fish VEGFB
    # Yellow_desaturated = NodeStyle()
    # Yellow_desaturated["bgcolor"] = "#ffffe6"
    # n4 = t.get_common_ancestor("VEGFB_Lates_calcarifer_XP_018515937.1", "VEGFBB_Danio_rerio_tr|E9QDL3|E9QDL3_DANRE")
    # n4.set_style(Yellow_desaturated)
    # # Mammalian VEGFA
    # Orange = NodeStyle()
    # Orange["bgcolor"] = "#ffa280"
    # n6 = t.get_common_ancestor("VEGFA_Ovis_aries_AAC23608.1", "VEGFA_Zalophus_californianus_XP_027459870.1")
    # n6.set_style(Orange)
    # # Fish VEGFA
    # Orange_desaturated = NodeStyle()
    # Orange_desaturated["bgcolor"] = "#ffece6"
    # n9 = t.get_common_ancestor("VEGFA_Lates_calcarifer_XP_018536247.1", "VEGFA_Cyprinus_carpio_XP_018949927.1")
    # n9.set_style(Orange_desaturated)
    # # Mammalian PGF
    # Violet = NodeStyle()
    # Violet["bgcolor"] = "#d98dbd"
    # n7 = t.get_common_ancestor("PGF_Zalophus_californianus_XP_027427310.1", "PGF_Homo_sapiens_NP_001193941.1")
    # n7.set_style(Violet)
    # # Fish PGF
    # Violet_desaturated = NodeStyle()
    # Violet_desaturated["bgcolor"] = "#d9c3d1"
    # n8 = t.get_common_ancestor("PGF_Lates_calcarifer_XP_018545725.1", "PGF_Cyprinus_carpio_XP_018975916.1")
    # n8.set_style(Violet_desaturated)
    # # BPSV
    # cyangrey = NodeStyle()
    # cyangrey["bgcolor"] = "#dae6ee"
    # n10 = t.get_common_ancestor("BPSV_Bos_taurus_AAS16477.1", "SePPV_Halichoerus_grypus_YP_009389414.1")
    # n10.set_style(cyangrey)

    t.set_outgroup(t & 'TGF-beta1')
    # Swapping children
    # n20 = t.get_common_ancestor("MCV_Lates_calcarifer_YP_009163772.1", "SePPV_Halichoerus_grypus_YP_009389414.1")
    # n20.swap_children()
    # n21 = t.get_common_ancestor("VEGFBA_Danio_rerio_tr|E7EZR4|E7EZR4_DANRE", "VEGFA_Homo_sapiens_NP_001020539.2")
    # n21.swap_children()
    # n22 = t.get_common_ancestor("VEGFB_Ovis_aries_XP_027815351.1", "VEGFA_Homo_sapiens_NP_001020539.2")
    # n22.swap_children()
    # n23 = t.get_common_ancestor("PGF_Cyprinus_carpio_XP_018975916.1", "VEGFA_Homo_sapiens_NP_001020539.2")
    # n23.swap_children()
    #n24 = t.get_common_ancestor("SePPV_Halichoerus_grypus_YP_009389414.1", "VEGFA_Homo_sapiens_NP_001020539.2")
    #n24.swap_children()
    #n25 = t.get_common_ancestor("ORFV_Rupicapra_rupicapra_QJD15044.1", "VEGFA_Ovis_aries_AAC23608.1")
    #n25.swap_children()
    #n26 = t.get_common_ancestor("ORFV_Capra_hircus_ASL69179.1", "SePPV_Halichoerus_grypus_YP_009389414.1")
    #n26.swap_children()

    ts = TreeStyle()
    ts.show_leaf_name = False
    # Zoom in x-axis direction
    ts.scale = 270
    # This makes all branches the same length!!!!!!!
    #ts.force_topology = True

    ts.show_branch_support = True
    ts.show_branch_length = False
    ts.draw_guiding_lines = True
    #ts.branch_vertical_margin = 10 # 10 pixels between adjacent branches
    ts.layout_fn = vegfe_tree_layout
    print(t)
    t.render(svg_file, tree_style = ts, units = "mm", h = 340)

    print('Drawing tree completed.')

if __name__ == '__main__':
    run()
