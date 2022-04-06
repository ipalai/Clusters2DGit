# This code creates LAMMPS input files for fluid clusters as described in the following paper:
#
# I. Palaia and A. Saric,
# Controlling cluster size in 2D phase-separating binary mixtures with specific interactions
# https://www.biorxiv.org/content/10.1101/2022.02.10.479877v1
#
import argparse
from make_initconfig import make_particles
import os
import numpy as np

def write_in_script(eps_pp, eps_ll, real, folderPattern, filePattern):

    filename = "Input_%s_%i.in" % (filePattern, real)

    f = open(filename, "w")


    f.write("units                  lj\n")
    f.write("dimension              2\n")
    f.write("atom_style             full  \n")
    f.write("boundary               p p p   \n")
    f.write("read_data              %s \n"%configName)

    f.write("neighbor               0.5 bin\n")
    f.write("neigh_modify           every 1 delay 1\n")


    f.write("pair_style             hybrid/overlay cosine/squared %lf lj/cut %lf\n" %(global_cutoff, global_cutoff))
    f.write("pair_modify            shift yes\n")
    f.write("pair_coeff             * * cosine/squared %lf %lf %lf wca \n" % (0, linker_sigma, linker_sigma))  ## repulsive only                   

    centreslist = [1, 2+q_A]
    for id in centreslist:
        for id2 in centreslist:
            if id2 >= id:
                f.write(f"pair_coeff             {id} {id2} lj/cut %lf %lf %lf\n" % (eps_pp, protein_sigma*2.0**(-1.0/6.0), protein_sigma ))  ## LJ purely repulsive (with lj/cut, sigma is not the same as cosine/squared wca, there's a factor 2**-1/6 difference)
                f.write(f"pair_coeff             {id} {id2} cosine/squared %lf %lf %lf\n" % (eps_pp_attr, protein_sigma, protein_sigma + protein_isotropic_attr_range))  ## cosine squared purely attractive

    ligands_A = range(2,q_A+2)
    ligands_B = range(q_A+3,q_A+q_B+3) 
    for iA in ligands_A:
        for iB in ligands_B:
            f.write(f"pair_coeff             {iA} {iB} cosine/squared {eps_ll} {linker_sigma} {ligand_range + linker_sigma} \n")               ## LAT pY132 - PLC nSH2


    f.write(f"velocity               all create {Temperature} 1 \n")

    f.write("log                    %s/Log_%s_%i.dat\n"%(folderPattern,filePattern,real))
    
    f.write(f"fix                    fLANG all langevin {Temperature} {Temperature} {LangevinDamping} %i\n"%seed)
    # f.write("fix                    fNVE all nve\n")
    f.write("neigh_modify           exclude molecule/intra all\n")
    f.write("fix                    rigidNVE all rigid/nve molecule\n")


    f.write("thermo                 %i\n"%dump_step)
    f.write("thermo_style           custom step temp press etotal epair\n")
    f.write("thermo_modify          flush yes\n")

    f.write("timestep               0.01\n")

    f.write("fix                    enforce_2d all enforce2d\n")
    f.write("dump                   1 all custom %i %s/Movie_%s_%i.xyz id type mol x y z \n"%(dump_step, folderPattern, filePattern, real))

    f.write("run                    %i\n"%run_steps)
    f.close()







if __name__ == "__main__":

    dump_step = 5e5
    eq_steps = 1e5
    run_steps = 5e7

    ## Simulation parameters

    eps_pp = 10.0 
    eps_ll = 10.0 
    
    num_A = 200
    num_B = 200
    density_A_target=0.03 #0.0075
    density_total=-1 # =-1 if simulation runs at constant A density (i.e. if nA and densA are given, while densB is defined by nB), >0 if simulation runs at constant total density

    q_A=4
    q_B=4

    protein_sigma = 2.0
    protein_isotropic_attr_range=1.0
    eps_pp_attr=1.0

    linker_sigma = 0.1
    dfactor=0.5*0.95
    d_pl=dfactor*protein_sigma
    ligand_range = 0.2
    thetadeg = 25

    Temperature = 1
    LangevinDamping = 1.0



    parser = argparse.ArgumentParser(description="Script for fluid clusters.")
    parser.add_argument('--eps_pp', '-eps_pp', dest='eps_pp', action='store', type=float, default=eps_pp, help='Energy depth of repulsion between proteins')
    parser.add_argument('--eps_ll', '-eps_ll', dest='eps_ll', action='store', type=float, default=eps_ll, help='Energy depth of attraction/repulsion between cross linkers')

    parser.add_argument('--num_A', '-num_A', '-nA', dest='num_A', action='store', type=int, default=num_A, help='num_A = number of A')
    parser.add_argument('--num_B', '-num_B', '-nB', dest='num_B', action='store', type=int, default=num_B, help='num_B = number of B')
    parser.add_argument('--density_A_target', '-density_A_target', '-dens', dest='density_A_target', action='store', type=float, default=density_A_target, help='density_A_target = surface density of A in units of sigma_p')
    parser.add_argument('--density_total', '-density_total', '-denstot', dest='density_total', action='store', type=float, default=density_total, help='density_total = surface density of A+B in units of sigma_p (overrides --density_A_target)')

    parser.add_argument('--q_A', '-q_A', '-qA', dest='q_A', action='store', type=int, default=q_A, help='q_A = valence of A')
    parser.add_argument('--q_B', '-q_B', '-qB', dest='q_B', action='store', type=int, default=q_B, help='q_B = valence of B')
    
    parser.add_argument('--protein_sigma', '-protein_sigma', dest='protein_sigma', action='store', type=float, default=protein_sigma, help='protein_sigma = sigma_p')
    parser.add_argument('--linker_sigma', '-linker_sigma', dest='linker_sigma', action='store', type=float, default=linker_sigma, help='linker_sigma = sigma_l')
    parser.add_argument('--ligand_range', '-ligand_range', dest='ligand_range', action='store', type=float, default=ligand_range, help='ligand_range = r_l')
    parser.add_argument('--dfactor', '-dfactor', dest='dfactor', action='store', type=float, default=dfactor, help='dfactor = position-of-crosslinker / sigma_p. (dfactor = 0.5 corresponds to the surface of the central particle)')

    parser.add_argument('--protein_isotropic_attr_range', '-ra', dest='protein_isotropic_attr_range', action='store', type=float, default=protein_isotropic_attr_range, help='Attraction range of the cosine-squared isotropic attraction between molecules')
    parser.add_argument('--eps_pp_attr', '-ea', dest='eps_pp_attr', action='store', type=float, default=eps_pp_attr, help='Strength of the cosine-squared isotropic attraction between molecules')

    parser.add_argument('-real','--real', dest='real', action='store', type=int, default=0, help='Number of realisation (statistics). Sets the seed.')
    parser.add_argument('-runsteps','--runsteps', dest='run_steps', action='store', type=int, default=run_steps, help='Number of time steps')
    parser.add_argument('-dumpstep','--dumpstep', dest='dump_step', action='store', type=int, default=dump_step, help='Number of skipped time steps in dump file')

    args = parser.parse_args()

    eps_pp = args.eps_pp
    eps_ll = args.eps_ll

    num_A = args.num_A
    num_B = args.num_B
    density_A_target = args.density_A_target
    density_total = args.density_total
    if density_total>0:
        density_A_target=density_total*num_A/(num_A+num_B)

    q_A = args.q_A
    q_B = args.q_B

    protein_sigma = args.protein_sigma
    linker_sigma = args.linker_sigma
    ligand_range = args.ligand_range
    dfactor = args.dfactor

    protein_isotropic_attr_range = args.protein_isotropic_attr_range
    eps_pp_attr = args.eps_pp_attr

    real = args.real
    run_steps = args.run_steps
    dump_step = args.dump_step
    seed = real+100

    global_cutoff = protein_sigma+protein_isotropic_attr_range


    ## Initial configuration

    num_tot = num_A + num_B
    side = np.ceil(np.sqrt(num_tot))
    lattice_factor = np.sqrt(1.0*num_A/(side*side)/density_A_target)       ## Ensures that A surface density is 0.0075 (protein_sigma)^-2
    density_A = num_A/(side*lattice_factor)**2                         ## In units of (protein_sigma)^-2

    system = make_particles(protein_sigma, linker_sigma, dfactor, num_A, q_A, num_B, q_B, side, lattice_factor, real)
    system.make_A()
    system.make_B()


    filePattern = f"2d_ljcossq_ep{eps_pp}_sp{protein_sigma}_ea{eps_pp_attr}_ra{protein_isotropic_attr_range}_el{eps_ll}_sl{linker_sigma}_rl{ligand_range}_d{dfactor}_nA{num_A}_qA{q_A}_nB{num_B}_qB{q_B}_T{Temperature}_Tdamp{LangevinDamping}_dA{density_A:.4f}"
    folderPattern = f"Results_el{eps_ll}_qA{q_A}_qB{q_B}"

    if not os.path.exists(folderPattern):
        os.makedirs(folderPattern)

    configName = "Config_%s_%i.dat" %(filePattern,real)

    header = ["LAMMPS Description \n \n",
              "\t " + str(system.numAll) + " atoms \n \t " + str(system.numBonds) +
              " bonds \n \t " + str(system.numAngles) + " angles \n \t 0 dihedrals \n \t 0 impropers \n",
              "\n \t "+str(system.numTypes)+" atom types \n \t 2 bond types \n \t 0 angle types \n \t 0 dihedral types \n \t 0 improper types \n",
              "\n \t " + str(-system.Lx*0.5) + " " + str(system.Lx*0.5) + " xlo xhi\n \t", str(-system.Ly*0.5) + " " + str(system.Ly*0.5) + " ylo yhi \n \t",
              str(-system.Lz*0.5) + " " + str(system.Lz*0.5) + " zlo zhi\n"]

    header.append("\nMasses \n \n")
    for i in range(len(system.type_mass_list)):
        header.append("\t %i %.4f \n"%(system.type_mass_list[i][0],system.type_mass_list[i][1]))

    f = open(configName, "w")

    for item in header:
        f.write("%s " % item)

    for item in system.coords:
        f.write("%s " % item)


    f.close()

    ## Write input script
    write_in_script(eps_pp, eps_ll, real, folderPattern, filePattern)
