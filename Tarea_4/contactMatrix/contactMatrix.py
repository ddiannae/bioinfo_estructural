import Bio.PDB
import numpy
import pylab
import seaborn as sns

pdb_code = "1AYI"
pdb_filename = "P1.pdb"

def calc_residue_dist(residue_one, residue_two) :
    """Returns the C-alpha distance between two residues"""
    min_dist =  10000.0
    for atom1 in residue_one.get_atoms():
        for atom2 in residue_two.get_atoms():
            diff_vector = atom1.coord - atom2.coord
            adist = numpy.sqrt(numpy.sum(diff_vector * diff_vector))
            if adist < min_dist:
                min_dist = adist
    return min_dist

def calc_dist_matrix(chain_one, chain_two, dim) :
    """Returns a matrix of C-alpha distances between two chains"""

    answer = numpy.zeros((dim, dim), numpy.float)
    for row, residue_one in enumerate(chain_one) :
        if residue_one.get_resname() != "HOH":
            for col, residue_two in enumerate(chain_two) :
                if residue_two.get_resname() != "HOH":
                    answer[row, col] = calc_residue_dist(residue_one, residue_two)
    return answer

structure = Bio.PDB.PDBParser().get_structure(pdb_code, pdb_filename)
model = structure[0]
dist_matrix = calc_dist_matrix(model["A"], model["A"], 86)

sns.distplot(dist_matrix, hist=True, kde=True)

pylab.matshow(numpy.transpose(dist_matrix))
pylab.colorbar()
pylab.show()

contact_map = dist_matrix < 15.0
pylab.matshow(numpy.transpose(contact_map))
pylab.colorbar()
pylab.show()

A = numpy.asarray(dist_matrix).reshape(-1)
sns.distplot(A, hist=True, kde=True,
             color = 'darkblue',
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 4})
