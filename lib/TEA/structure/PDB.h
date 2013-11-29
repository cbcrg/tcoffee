/*
 * PDB.h
 *
 *  Created on: Feb 16, 2012
 *      Author: Carsten Kemena
 *
 *   Copyright 2011 Carsten Kemena
 *
 *
 * This file is part of BioTools++.
 *
 * BioTools++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BioTools++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BioTools++.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef PDB_H_
#define PDB_H_

#include <vector>
#include <string>


namespace BioTools {
namespace Structure {

typedef struct
{
	int atom_id;
	std::string name;
	double x;
	double y;
	double z;
} Atom;


typedef struct
{
	int amino_id;
	std::string name;
	std::vector<Atom> atoms;

} AminoAcid;

typedef struct
{
	std::vector<AminoAcid> amino_acids;
} Model;

typedef struct
{
	char chain_id;
	bool use_model;
	std::vector<Model> models;
} Chain;





class PDB {

private:


public:
	PDB();
	virtual ~PDB();

	std::vector<Chain> chains;
};

} /* namespace Structure */
} /* namespace BioTools */
#endif /* PDB_H_ */
