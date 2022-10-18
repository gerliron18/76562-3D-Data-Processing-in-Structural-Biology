import random
import statistics
import sys

from Bio.PDB import PDBList, PDBParser, Superimposer, MMCIFIO

MIN_NUM_OF_ARGS = 5
NUM_OF_ARGS_WITH_BONUS = 6

PDB_ID1_INDEX = 1
CHAIN1_INDEX = 2
PDB_ID2_INDEX = 3
CHAIN2_INDEX = 4
BONUS_INDEX = 5

BONUS_9_1_NUM_OF_RUNS = 1000
BONUS_9_1_PERC_OF_OUTLIERS = 0.1


def get_ca_atoms(chain):
    """
    Filter out all atoms which aren't CA from all the residues
    :param chain: The chain to work on
    :return: A list with all the CA atoms
    """
    ca_list = []
    for residue in chain.get_residues():
        if "CA" in residue:
            ca_list.append(residue["CA"])

    return ca_list


def get_structure(pdb_id):
    """
    Download and get the structure of a pdb file
    :param pdb_id: The id of the pdb file
    :return: The pdb obj
    """
    pdb_list = PDBList()
    pdb_parser = PDBParser(PERMISSIVE=1, QUIET=True)
    pdb_f = pdb_list.retrieve_pdb_file(pdb_id, file_format="pdb")
    pdb = pdb_parser.get_structure(pdb_id, pdb_f)
    return pdb


def align_chains(pdb1_id, pdb2_id, chain1, chain2):
    """
    A Python script that takes two PDB identifiers and chain identifiers
    as input. The script will retrieve the appropriate files from the PDB, and
    align the chains from each file by their set of C-alpha (“CA”) atoms (filtering
    out residues that lack CA atoms). The output will include the RMSD
    between the C-alpha atoms (CA) of the two structures, and the aligned
    PDB files in mmCIF format. It is assumed the two chains have the same
    number of CA atoms, in the same order
    :param pdb1_id: The id of pdb file 1
    :param pdb2_id: The id of pdb file 2
    :param chain1: The chain of pdb 1
    :param chain2: The chain of pdb 2
    :return: None
    """
    # Init the objects
    sup = Superimposer()
    io = MMCIFIO()

    # Download and read the pdb files
    pdb1 = get_structure(pdb1_id)
    pdb2 = get_structure(pdb2_id)

    # Extract the chains and atoms from the file
    chain1 = pdb1[0][chain1]
    chain2 = pdb2[0][chain2]
    ca1 = get_ca_atoms(chain1)
    ca2 = get_ca_atoms(chain2)

    # Align the atoms, calculate the RMSD
    sup.set_atoms(ca1, ca2)
    sup.apply(pdb2[0].get_atoms())

    print("The RMSD value for alignment of %s chain %s and %s chain %s is %s" %
          (pdb1_id, chain1.id, pdb2_id, chain2.id, sup.rms))

    # Save the align files
    io.set_structure(pdb1)
    io.save("%s.cif" % pdb1_id)
    io.set_structure(pdb2)
    io.save("%s.cif" % pdb2_id)


# BONUS 9.1 IMPLEMENTED
def bonus_9_1(pdb1_id, pdb2_id, chain1, chain2):
    """
    Add an option to use your alternative measure for RMSD from 8e, or a
    different measure
    We decided to use the second option we gave for RMSD in 8e which mean taking the mean of several runs
    of RMSD without 10% of the atoms to try and ignore outliers
    :param pdb1_id: The id of pdb file 1
    :param pdb2_id: The id of pdb file 2
    :param chain1: The chain of pdb 1
    :param chain2: The chain of pdb 2
    :return: None
    """
    # Init the objects
    sup = Superimposer()

    # Download and read the pdb files
    pdb1 = get_structure(pdb1_id)
    pdb2 = get_structure(pdb2_id)

    # Extract the chains and atoms from the file
    chain1 = pdb1[0][chain1]
    chain2 = pdb2[0][chain2]
    ca1 = get_ca_atoms(chain1)
    ca2 = get_ca_atoms(chain2)

    size_to_take = int(len(ca1) * (1 - BONUS_9_1_PERC_OF_OUTLIERS))
    indexes_list = [i for i in range(len(ca1))]
    rmsd_list = []

    # Align the atoms, calculate the RMSD for BONUS_9_1_NUM_OF_RUNS times
    for i in range(BONUS_9_1_NUM_OF_RUNS):
        temp_index = random.sample(indexes_list, size_to_take)
        temp_ca1 = [ca1[index] for index in temp_index]
        temp_ca2 = [ca2[index] for index in temp_index]
        sup.set_atoms(temp_ca1, temp_ca2)
        rmsd_list.append(sup.rms)

    print("The RMSD using our system is: %s" % (statistics.mean(rmsd_list)))


def main():
    # Check that we have enough arguments
    if len(sys.argv) < MIN_NUM_OF_ARGS:
        print("Not enough arguments to run the script\n")
        print("my_align.py <pdbid1> <chain1> <pdbid2> <chain2>")
        sys.exit(-1)

    # Extract the argument s
    chain1 = sys.argv[CHAIN1_INDEX]
    chain2 = sys.argv[CHAIN2_INDEX]
    pdb1 = sys.argv[PDB_ID1_INDEX]
    pdb2 = sys.argv[PDB_ID2_INDEX]

    align_chains(pdb1, pdb2, chain1, chain2)

    if len(sys.argv) == NUM_OF_ARGS_WITH_BONUS and sys.argv[BONUS_INDEX] == "True":
        bonus_9_1(pdb1, pdb2, chain1, chain2)


if __name__ == '__main__':
    main()
