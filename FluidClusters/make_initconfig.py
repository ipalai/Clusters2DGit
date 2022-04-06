import numpy as np


class make_particles(object):


    def __init__(self, sigma_proteins, sigma_linkers, dfactor, num_A, q_A, num_B, q_B, side, lattice_factor, real):

        self.num_A = num_A
        self.num_B = num_B
        self.num_tot = self.num_A + self.num_B

        self.q_A = q_A
        self.q_B = q_B

        self.sigma_proteins = sigma_proteins
        self.sigma_linkers = sigma_linkers
        self.dfactor = dfactor

        self.lattice_factor = lattice_factor
        self.side = side
     
        self.lattice_constant = lattice_factor*self.sigma_proteins
        self.Lx = self.side*self.lattice_constant
        self.Ly = self.side*self.lattice_constant
        self.Lz = 0.1


        self.coords = ["\nAtoms \n \n"]
        self.bonds = ["\nBonds \n \n"]
        self.angles= ["\nAngles \n \n"]

        # m is number of atoms and k is number of 3 atom-molecules
        self.numAll = 0
        self.k = 0

        self.numBonds = 0
        self.bondId = 0

        self.numAngles = 0
        self.angleId = 0

        self.numTypes = 0
        self.type_mass_list =[]


        def init_lattice_2d(lat_con):
            num_x = int(self.side)
            num_y = int(self.side)

            lattice_out = np.zeros(shape=(num_x*num_y, 3))
            counter = 0
            for i in range(num_x):
                for j in range(num_y):

                    lattice_out[counter, 0] = lat_con * i - lat_con*num_x * 0.5
                    lattice_out[counter, 1] = lat_con * j - lat_con*num_y * 0.5
                    lattice_out[counter, 2] = 0.0
                    counter += 1

            return lattice_out



        self.lattice_sites = init_lattice_2d(self.lattice_constant)

        assert self.lattice_factor>1.1, "ERROR: lattice_factor<1.1, too many particles."

        assert self.num_tot <= self.lattice_sites.shape[0], "ERROR: not enough lattice sites, increase box size."

        lattice_inds = np.arange(self.lattice_sites.shape[0])

        np.random.seed(100+real)
        self.start_inds_A = np.random.choice(lattice_inds, size = self.num_A, replace = False)

        self.lattice_inds_remain = np.array([item for item in lattice_inds if item not in set(list(self.start_inds_A))])

        self.start_inds_B = np.random.choice(self.lattice_inds_remain, size = self.num_B, replace = False)

        self.lattice_inds_remain = np.array([item for item in lattice_inds if item not in set(list(self.start_inds_A) +  list(self.start_inds_B))])

        test_list = list(self.start_inds_A) + list(self.start_inds_B)
        assert len(test_list) == self.num_tot


    def make_A(self):

        particles_list = range(1,self.q_A+2)
        particle_types = len(particles_list)
        self.numTypes += len(np.unique(particles_list))

        for i in range(particle_types):
             to_add = [particles_list[i], 1]
             if to_add not in self.type_mass_list:
                 self.type_mass_list.append(to_add)

        for i in range(len(self.start_inds_A)):  ## go through chains
            indBuf = self.start_inds_A[i]

            self.k += 1
            x = self.lattice_sites[indBuf, 0]
            y = self.lattice_sites[indBuf, 1]
            z = self.lattice_sites[indBuf, 2]

            for n in range(self.q_A+1):
                self.numAll += 1
                if (n == 0):

                    self.coords.append(
                        "\t " + str(self.numAll) + " " + str(self.k) + " "+str(particles_list[n])+" 0 " + str(x) + " " + str(y) + " " + str(
                            z) + " 0 0 0 \n")

                radial_dist = self.dfactor*self.sigma_proteins

                position_on_circle = radial_dist*np.array([[np.cos(i*2*np.pi/self.q_A), np.sin(i*2*np.pi/self.q_A)] for i in range(0,self.q_A)])
                if n > 0:

                    self.coords.append(
                        "\t " + str(self.numAll) + " " + str(self.k) + " " +  str(particles_list[n]) + " 0 " + str(x - position_on_circle[n-1,0]) + " " + str(y - position_on_circle[n-1,1]) + " " + str(
                            z) + " 0 0 0 \n")



    def make_B(self):

        particles_list = range(self.q_A+2,self.q_A+self.q_B+3)  ## first one is central body
        particle_types = len(particles_list)
        self.numTypes += len(np.unique(particles_list))
   
        for i in range(particle_types):
             self.type_mass_list.append([particles_list[i], 1])

        for i in range(len(self.start_inds_B)):  ## go through chains
            indBuf = self.start_inds_B[i]

            self.k += 1
            x = self.lattice_sites[indBuf, 0]
            y = self.lattice_sites[indBuf, 1]
            z = self.lattice_sites[indBuf, 2]

            radial_dist = self.dfactor*self.sigma_proteins
            position_on_circle = radial_dist*np.ones(shape = (self.q_B,2 ))
            ligand_angles = 2.0*np.pi*np.array([i/self.q_B for i in range(self.q_B)])
            for i in range(len(ligand_angles)):
                position_on_circle[i, 0] *= np.sin(ligand_angles[i])
                position_on_circle[i, 1] *= np.cos(ligand_angles[i])


            for n in range(particle_types):
                self.numAll += 1
                if (n == 0):
                    self.coords.append(
                        "\t " + str(self.numAll) + " " + str(self.k) + " " + str(
                            particles_list[n]) + " 0 " + str(x) + " " + str(y) + " " + str(
                            z) + " 0 0 0 \n")

                if n > 0:
                    self.coords.append(
                        "\t " + str(self.numAll) + " " + str(self.k) + " " + str(
                            particles_list[n]) + " 0 " + str(
                            x - position_on_circle[n - 1, 0]) + " " + str(
                            y - position_on_circle[n - 1, 1]) + " " + str(
                            z) + " 0 0 0 \n")
