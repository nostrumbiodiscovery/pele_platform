import argparse, os, sys
from string import *
import glob
from math import sqrt

def randomize_starting_position(clean_ligand_pdb, input_ligand, ligname, rec_file, rec_com, lig_com, env):
    ### Randomize ligand starting position for outside-inside
    import pymol
    import numpy
    contact = -1
    n = 0
    n0 = 0
    if not lig_com :
        rec_com = True
    pymol.cmd.load(clean_ligand_pdb, 'ligand') 
    pymol.cmd.save(input_ligand,'ligand')
    #clean_ligand_pdb, origin = clean_lig_pdb(input_ligand) ## Do this step again in case Pymol changes PDB names (maybe not needed?)
    #pymol.cmd.delete('ligand')
    #pymol.cmd.load(, 'ligand') #Reload ligand
    COM_lig=pymol.cmd.centerofmass('ligand')
    if lig_com :
        print("Ligand COM is %s" %COM_lig)
        print("Sampling 25A spherical box centered around ligand")
        D = 25.0
        sphere_cent = COM_lig

    pymol.cmd.load(rec_file, 'receptor') 
    if rec_com :
        COM=pymol.cmd.centerofmass('receptor')
        pymol.cmd.create('rec_zero','receptor') 
        pymol.cmd.translate((-1*COM[0],-1*COM[1],-1*COM[2]),'rec_zero') # Move rec to the origin
        rec_min_max = pymol.cmd.get_extent('rec_zero')
        maxdist1 = numpy.sqrt(rec_min_max[0][0]**2 + rec_min_max[0][1]**2 + rec_min_max[0][2]**2)
        maxdist2 = numpy.sqrt(rec_min_max[1][0]**2 + rec_min_max[1][0]**2 + rec_min_max[1][0]**2)
        maxdist = numpy.maximum(maxdist1,maxdist2)
        print("Receptor COM is %s" %COM)
        print("Receptor max sphere radius is %s" %maxdist)
        D=numpy.ceil(6.0+maxdist) #radius of the sphere from the origin
        print("Sampling %sA spherical box centered around receptor COM" %D)
        sphere_cent = COM

    output = []
    while (n < 10) :
        n0 += 1
        phi = numpy.random.uniform(0,2*numpy.pi)
        costheta = numpy.random.uniform(-1,1)
        u = numpy.random.uniform(0,1)

        theta = numpy.arccos( costheta )
        r = D * numpy.cbrt( u )
        x = r * numpy.sin( theta) * numpy.cos( phi )
        y = r * numpy.sin( theta) * numpy.sin( phi )
        z = r * numpy.cos( theta )
                
        translation = (x,y,z)
        pymol.cmd.delete('ligand')
        pymol.cmd.load(clean_ligand_pdb, 'ligand') #Reload ligand
        pymol.cmd.translate((-1*COM_lig[0],-1*COM_lig[1],-1*COM_lig[2]),'ligand') # Move ligand to the origin
        pymol.cmd.translate(translation,'ligand')

        # Make sure new ligand position is at least 5A but not more than 8A from any protein atom
        # Also make sure it is within the sampling sphere 
        pymol.cmd.translate((sphere_cent[0],sphere_cent[1],sphere_cent[2]),'ligand') # Move ligand back to the original frame
        COM_lig_new=pymol.cmd.centerofmass('ligand')
        dist = numpy.sqrt( ((COM_lig_new[0] - sphere_cent[0])**2) + ((COM_lig_new[1] - sphere_cent[1])**2)  + ((COM_lig_new[2] - sphere_cent[2])**2) )
        if dist < D :
            rotation = (numpy.random.uniform(0,360),numpy.random.uniform(0,360),numpy.random.uniform(0,360)) #add a random rotation
            pymol.cmd.rotate('x',rotation[0],'ligand') #add a random rotation
            pymol.cmd.rotate('y',rotation[1],'ligand') #add a random rotation
            pymol.cmd.rotate('z',rotation[2],'ligand') #add a random rotation
            contact_near=pymol.cmd.select('contact','(receptor and (ligand around 5))')
            contact_far=pymol.cmd.select('contact','(receptor and (ligand around 8))')
            if int(contact_near) == 0 and int(contact_far) > 0:
                if env:
                    rand_lig_pdb = os.path.join(env.pele_dir, '%s_rand_position0%d.pdb' %(ligname,n))
                else:
                    rand_lig_pdb = os.path.join(".", '%s_rand_position0%d.pdb' %(ligname,n))
                pymol.cmd.save(rand_lig_pdb,'ligand')
                output.append(rand_lig_pdb)
                n += 1

    return output

def join(receptor, ligands, env, output="input{}.pdb"):

    """
    Join receptor&ligand pdb in one conserving old formatting
    and not repiting atomnumbers
    """
    
    with open(receptor, "r") as f:
        receptor_content = f.readlines()
        initial_atomnum = max([ int(line[6:11]) for line in receptor_content if line.startswith("ATOM") or line.startswith("HETATM")])

    outputs = []
    for i, ligand in enumerate(ligands):
        with open(ligand, "r") as fin:
            #exclude connects
            ligand_content = [line for line in fin if line.startswith("ATOM") or line.startswith("HETATM")]
            
            current_atomnum = initial_atomnum +1
            for j, line in enumerate(ligand_content):
                ligand_content[j] = line[:6] + "{:>5}".format(current_atomnum) + line[11:]
                current_atomnum += 1
                
        content_join_file = receptor_content + ligand_content + ["END"]
        if env:
            output_path = os.path.join(env.pele_dir, output.format(i))
        else:
            output_path = os.path.join(".", output.format(i))
        with open(output_path, "w") as fout:
            fout.write("".join(content_join_file))
        outputs.append( output_path )

    return outputs

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ligand", type=str, required=True, help="Ligand pdb file")
    parser.add_argument("--receptor", type=str, required=True, help="Receptor pdb file")
    parser.add_argument("--resname", type=str, required=True, help="Ligand resname")
    args = parser.parse_args()
    return os.path.abspath(args.ligand), os.path.abspath(args.receptor), args.resname

if __name__ == "__main__":
   ligand, receptor, resname  = parse_args()
   output = randomize_starting_position(ligand, "ligand_input.pdb", resname,
           receptor, None, None, None)
   join(receptor, output, None)

