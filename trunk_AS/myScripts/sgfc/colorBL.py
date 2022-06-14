# This pymol script will color FtsB and FtsL appropriately
from pymol import cmd

# Functions


def count_mols_in_sel(sel="sele"):
    # Returns the number of distinct molecules in a given selection

    sel_copy = "__selcopy"
    cmd.select(sel_copy, sel)
    num_objs = 0
    atoms_in_sel = cmd.count_atoms(sel_copy)

    while atoms_in_sel > 0:
        num_objs += 1
        cmd.select(sel_copy, "%s and not (bm. first %s)" % (sel_copy, sel_copy))
        atoms_in_sel = cmd.count_atoms(sel_copy)

    print "There are %d distinct molecules in the selection '%s'." % (num_objs, sel) 
    return num_objs

# Main
def color_BL():
    obj_list = cmd.get_names('objects')
    ftsb_color = 'slate'
    ftsb_hotspot_color = 'cyan'
    ftsb_hotspot = [15, 12, 16, 5, 19] # Coevolutionary hotspots


    ftsl_color = 'paleyellow'
    ftsl_hotspot_color = 'orange'
    ftsl_hotspot = [46, 45, 49, 39, 56, 52, 53] # Coevolutionary hotspots

    chains = cmd.get_chains(obj_list[0])
    numChains = len(chains)

    # chainIds = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K',
    #             'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V',
    #             'W', 'X', 'Y', 'Z']

    # Make the selections for ftsb and ftsl
    ftsb_chains = ""
    ftsl_chains = ""
    for i in range(0, numChains/2):
        ftsb_chains += chains[i] + "+"

    ftsb_chains = ftsb_chains[:-1]
    print ftsb_chains

    for i in range(numChains/2, numChains):
        ftsl_chains += chains[i] + "+"

    ftsl_chains = ftsl_chains[:-1]
    print ftsl_chains

    cmd.select("ftsb", "chain " + ftsb_chains)
    cmd.select("ftsl", "chain " + ftsl_chains)

    # Color the chains appropriately
    cmd.color(ftsb_color, "ftsb and elem c")
    cmd.color(ftsl_color, "ftsl and elem c")

    # Select and color the hotspot residues
    ftsb_hotspot_sele = "ftsb and resi "
    for i in ftsb_hotspot:
        ftsb_hotspot_sele += str(i) + "+"
    ftsb_hotspot_sele = ftsb_hotspot_sele[:-1]

    ftsl_hotspot_sele = "ftsl and resi "
    for i in ftsl_hotspot:
        ftsl_hotspot_sele += str(i) + "+"
    ftsl_hotspot_sele = ftsl_hotspot_sele[:-1]

    cmd.select("ftsb_hotspot", ftsb_hotspot_sele)
    cmd.select("ftsl_hotspot", ftsl_hotspot_sele)
    cmd.color(ftsb_hotspot_color, "ftsb_hotspot and elem c")
    cmd.color(ftsl_hotspot_color, "ftsl_hotspot and elem c")

    # Select and color the coiled-coil heptads
    cmd.select("ftsb_a", "ftsb and resi 22+29+36+43+50+57+64+71+78")
    cmd.select("ftsb_b", "ftsb and resi 23+30+37+44+51+58+65+72+79")
    cmd.select("ftsb_c", "ftsb and resi 24+31+38+45+52+59+66+73+80")
    cmd.select("ftsb_d", "ftsb and resi 25+32+39+46+53+60+67+74+81")
    cmd.select("ftsb_e", "ftsb and resi 26+33+40+47+54+61+68+75+82")
    cmd.select("ftsb_f", "ftsb and resi 27+34+41+48+55+62+69+76+83")
    cmd.select("ftsb_g", "ftsb and resi 21+28+35+42+49+56+63+70+77")

    cmd.select("ftsl_a", "ftsl and resi 53+60+67+74+81+88+95+102+109")
    cmd.select("ftsl_b", "ftsl and resi 54+61+68+75+82+89+96+103+110")
    cmd.select("ftsl_c", "ftsl and resi 55+62+69+76+83+90+97+104+111")
    cmd.select("ftsl_d", "ftsl and resi 56+63+70+77+84+91+98+105+112")
    cmd.select("ftsl_e", "ftsl and resi 57+64+71+78+85+92+99+106+113")
    cmd.select("ftsl_f", "ftsl and resi 58+65+72+79+86+93+100+107+114")
    cmd.select("ftsl_g", "ftsl and resi 52+59+66+73+80+87+94+101+108")

    cmd.color("palegreen", "(ftsb_a or ftsb_d or ftsl_a or ftsl_d) and elem c")

cmd.extend("colorBL", color_BL)
