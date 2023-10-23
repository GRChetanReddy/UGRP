import time, os, pandas
import peptideutils as pu

import Bio.PDB
from PeptideBuilder import Geometry
import PeptideBuilder

def readin(fpath):
	f = open(fpath)
	content = f.read()
	f.close()
	return content

def calcAP(peptide):
    try:
        APs = pandas.read_csv("/users/rkb19187/Desktop/AP_50ns.txt", sep=" ", index_col=0, header=None)
    except: #It might be being used right now
        time.sleep(100)
        APs = pandas.read_csv("/users/rkb19187/Desktop/AP_50ns.txt", sep=" ", index_col=0, header=None)
    if pu.translate1to3(peptide) in APs.index:
        return APs.at[pu.translate1to3(peptide), 1]
    
    if len(peptide) == 3:
        root = "/users/rkb19187/Desktop/APMD/tripeptides/"
    elif len(peptide) == 4:
        root = "/users/rkb19187/Desktop/APMD/tetrapeptides/"
    elif len(peptide) == 5:
        root = "/users/rkb19187/Desktop/APMD/pentapeptides/"
    elif len(peptide) == 6:
        root = "/users/rkb19187/Desktop/APMD/hexapeptides/"
    elif len(peptide) == 7:
        root = "/users/rkb19187/Desktop/APMD/heptapeptides/"
    elif len(peptide) == 8:
        root = "/users/rkb19187/Desktop/APMD/octapeptides/"
    if not os.path.exists(root+peptide+"/"+"SASA.xvg"):
        print("Peptide MD experiment not done/complete (SASA.xvg not found)".format(peptide))
        return False
    
    if not os.path.exists(root+peptide+"/"+peptide+"_eq.gro"):
        print("No {}_eq.gro file, calcAP failed".format(peptide))
        return False
    
    if not os.path.exists(root+peptide+"/SASA.xvg"):
        if not os.path.exists(root+peptide+"/"+peptide+"_startend.gro"):
            cwd = os.getcwd()
            os.chdir("/users/rkb19187/Desktop/APMD")
            os.system("python ConcatGRO.py "+peptide)
            os.chdir(cwd)
        cwd = os.getcwd()
        os.chdir(root+peptide)
        os.system("echo Protein | gmx_mpi sasa -f {}_startend.gro -s {}_eq.tpr -o SASA.xvg".format(peptide, peptide))
        os.chdir(cwd)
        
    csv_SASA = readin(root+peptide+"/SASA.xvg")
    initial = -1
    AP = -1
    for line in csv_SASA.split("\n"):
        if len(line) > 0:
            if line[0] != "#" and line[0] != "@":
                line = line.split(" ")
                for i in range(10000):
                    try:
                        line.remove("")
                    except:
                        break
                if len(line) < 2:
                    print("calc_AP: len(line) < 2 for", peptide, "(", line, ")")
                    AP = -1
                if initial == -1:
                    initial = float(line[1])
                if line[0] == "50000.000":
                    AP = initial / float(line[1])
                    break
    if AP == -1:
        print("Didn't get to t = 50 ns when reading SASA.xvg")
        os.remove(root+peptide+"/SASA.xvg")
        return False
    return AP

def mkpdb(sequence):
    for index in range(len(sequence)):
        letter = sequence[index]
        geo = Geometry.geometry(letter)
        #geo.phi=-180
        #geo.psi_im1=180
        if index == 0:
            structure = PeptideBuilder.initialize_res(geo)
        else:
            structure = PeptideBuilder.add_residue(structure, geo)
    out = Bio.PDB.PDBIO()
    out.set_structure(structure)
    out.save(sequence+"_aa.pdb")

def vmd_atomistic(vmd, sequence):
    cmd = vmd+" -dispdev text -e {}.pgn".format(sequence)
    os.system(cmd)

def count_water(fname):
    water_residues, water_count = set(('W', 'WF')), 0
    with open(fname) as gro_file:
    	for line in gro_file:
    		try:
    			resname = line[5:10].rstrip(' ')
    		except ValueError:
    			pass
    		else:
    			if resname in water_residues:
    				water_count += 1
    return water_count

def replace_in_file(fin, fout, original, new):
    if os.name == 'nt':
        cmd = "powershell -Command \"(gc {}) -replace '{}', '{}' | Out-File -encoding ASCII {}\" " .format(fin, original, new, fout)
    else:
        if fin == fout:
            cmd = "sed -i \"s/{}/{}/g\" {}".format(original, new, fin)
        else:
            cmd = "sed \"s/{}/{}/g\" {} > {}".format(original, new, fin, fout)
    print("replace_in_file:", cmd)
    os.system(cmd)

def runMD(peptide, cpus=4):
    APs = pandas.read_csv("/users/rkb19187/Desktop/AP_50ns.txt", sep=" ", index_col=0, header=None)
    if pu.translate1to3(peptide) in APs.index:
        print("runMD: peptide in APs.index")
        return (pu.translate1to3(peptide), APs.at[pu.translate1to3(peptide), 1])
    
    root = "/users/rkb19187/Desktop/APMD/"
    odir = os.getcwd()
    os.chdir(root)

    if "-" in peptide:
        peptide = pu.translate3to1(peptide)
    arguements = ""
    for letter in list(peptide):
        arguements = arguements + letter + "  "

    cmd = "sh apmd.sh "+arguements
    print(cmd)
    os.system(cmd)
    
    os.chdir(odir)
    
    return peptide