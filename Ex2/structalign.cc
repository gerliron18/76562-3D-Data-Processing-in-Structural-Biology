/**
 * Dror Bar (drordod) 203523352
 * Liron Gershuny (gerliron18) 308350503
 */


#include "Vector3.h"
#include "Atom.h"
#include "Triangle.h"
#include "RigidTrans3.h"
#include "Matrix3.h"
#include "Molecule.h"
#include "PDB.h"
#include "Match.h"
#include "GeomHash.h"
#include <chrono>
#include <iostream>

int main(int argc, char *argv[])
{
	// measure the run time
	auto start = std::chrono::system_clock::now();

	if (argc != 4)
	{
		std::cerr << "Usage: " << argv[0] << " epsilon pdb1 pdb2" << std::endl;
		exit(1);
	}

	//********Parameters********************
	float m_fDistThr = atof(argv[1]); // distance threshold on atoms in correspondence

	// read the two files into Molecule
	std::ifstream fileModel(argv[3]);
	std::ifstream fileTarget(argv[2]);

	if (!fileModel)
	{
		std::cout << "File " << argv[2] << "does not exist." << std::endl;
		return 0;
	}
	if (!fileTarget)
	{
		std::cout << "File " << argv[3] << "does not exist." << std::endl;
		return 0;
	}

	Molecule<Atom> molModel, molTarget;
	molModel.readPDBfile(fileModel, PDB::CAlphaSelector());
	molTarget.readPDBfile(fileTarget, PDB::CAlphaSelector());

	fileModel.close();
	fileTarget.close();

	//	Check if we have Ca or not, if not than probably RNA
	if (molModel.size() == 0)
	{
		std::ifstream fileModel(argv[3]);
		std::ifstream fileTarget(argv[2]);
		molModel.readPDBfile(fileModel, PDB::PSelector());
		molTarget.readPDBfile(fileTarget, PDB::PSelector());
		fileModel.close();
		fileTarget.close();
	}

	// next we insert the target molecule into hash
	// this will help us to find atoms that are close faster
	GeomHash<Vector3, int> gHash(3,
								 m_fDistThr); // 3 is a dimension and m_fDistThr is the size of the hash cube
	for (unsigned int i = 0; i < molTarget.size(); i++)
	{
		gHash.insert(molTarget[i].position(), i); // coordinate is the key to the hash, we store
		// atom index
	}

	// now we try random rotations and choose the best alignment from random rotations
	unsigned int iMaxSize = 0;
	float maxRMSD = 0;
	RigidTrans3 rtransBest, rtransTemp;
	Triangle target_tr, model_tr;


	for (unsigned int i = 0; i < molTarget.size() - 2; i++)
	{
		target_tr = Triangle(molTarget[i], molTarget[i + 1], molTarget[i + 2]);
		for (unsigned int j = 0; j < molModel.size() - 2; j++)
		{
			Match match;

			model_tr = Triangle(molModel[j], molModel[j + 1], molModel[j + 2]);
			rtransTemp = target_tr | model_tr;

			// apply rotation on each atom in the model molecule and
			// add the pairs of atoms (one from target and one from model)
			// that are close enough to the match list
			for (unsigned int i = 0; i < molModel.size(); i++)
			{
				Vector3 mol_atom = rtransTemp * molModel[i].position(); // transformation

				// find close target molecule atoms using the hash
				HashResult<int> result;
				gHash.query(mol_atom, m_fDistThr, result); // key is mol atom coordinate

				// check if the atoms in the result are inside the distance threshold
				// the hash is a cube shape, there can be atoms further that the threshold
				for (auto x = result.begin(); x != result.end(); x++)
				{
					float dist = mol_atom.dist(molTarget[*x].position());
					if (dist <= m_fDistThr)
					{
						float score = (1 / (1 + dist));
						match.add(*x, i, score, score);
					}
				}
				result.clear();
			}

			//calculates transformation that is a little better than "rotation"
			match.calculateBestFit(molTarget, molModel);

			if (iMaxSize < match.size())
			{
				iMaxSize = match.size();
				maxRMSD = match.rmsd();
				rtransBest = match.rigidTrans();
			}
		}
	}

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;

	// End of the loop print the match size, rmsd and the transformation
	std::cout << iMaxSize << " " << maxRMSD << " " << rtransBest;

	// Save the pdb2 after transformation in a new file
	std::ifstream fromPath(argv[3]);
	Molecule<Atom> fromModel;
	fromModel.readAllPDBfile(fromPath, PDB::Selector());
	fromPath.close();

	fromModel.rigidTrans(rtransBest);
	std::ofstream target_f;
	target_f.open("transformed.pdb");
	target_f << fromModel;
	target_f.close();
}
