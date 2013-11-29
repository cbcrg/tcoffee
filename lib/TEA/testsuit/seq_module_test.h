/*
 * aln_test.cpp
 *
 *  Created on: Oct 13, 2011
 *      Author: Carsten Kemena
 *
 *   Copyright 2011 Carsten Kemena
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


// C header
#include <cstdlib>

// C++ header
#include <vector>
#include <string>
#include <utility>

// CxxTest header
#include <cxxtest/TestSuite.h>

#include "../src/lib/Sequence/Alignment.h"
#include "../src/lib/utils/ScoringMatrix.h"
#include "../src/lib/Sequence/aln_analysis.h"

using namespace std;
using namespace BioTools;
using namespace BioTools::Seq;
using namespace BioTools::Utils;


class AlignmentTest : public CxxTest::TestSuite
{
public:

	void test_FASTA()
	{
		printf("\nRunning alignment IO tests\n");
	//	TS_TRACE("FASTA format handling");
		Alignment aln;
		aln.read("../../testsuit/data/aln.fa");
		TS_ASSERT_EQUALS(aln.seq(1)->sequence(), "MAHARVLLLALAVLATAAVAVASSSSFADSNPIRPVTDRAASTLESAVLGALGRTRHALRFARFAVRYGKSYESAAEVRRRFRIFSESLEEVRSTN----RKGLPYRLGINRFSDMSWEEFQATRL-GAAQTCSATLAGNHLMRDA--AALPETKDWREDG-IVSPVKNQAHCGSCWTFSTTGALEAAYTQATGKNISLSEQQLVDCAGGFNNF--------GCNGGLPSQAFEYIKYNGGIDTEESYPYKGVNGV-CHYKAENAAVQVLDSVNITLNAEDELKNAVGLVRPVSVAFQVIDGFRQYKSGVYTSDHCGTTPDDVNHAVLAVGYGVENGV-----PYWLIKNSWGADWGDNGYFKMEMGKNMCAIATCASYPVVAA");

		Alignment aln2;
		vector<string> names;
		names.push_back("ALEU_HORVU");
		names.push_back("CATH_HUMAN");
		aln2.read("../../testsuit/data/aln.fa",names);
		TS_ASSERT_EQUALS(aln2.n_seqs(), 2);
		TS_ASSERT_EQUALS(aln2.seq(1)->name(), "CATH_HUMAN");

		string t = "------MWATLPLLCAGAWLLGV--------PVCGAAELSVNSLEK------------FHFKSWMSKHRKTY-STEEYHHRLQTFASNWRKINAHNNGNHTFKMALNQFSDMSFAEIKHKYLWSEPQNCSAT--KSNYLRGTGPYPPSVDWRKKGNFVSPVKNQGACGSCWTFSTTGALESAIAIATGKMLSLAEQQLVDCAQDFNNYGCQGGLPSQAFEYILYNKGIMGEDTYPYQGKDGYCKFQPGKAIGFVKDVANITIYDEEAMVEAVALYNPVSFAFEVTQDFMMYRTGIYSSTSCHKTPDKVNHAVLAVGYGEKNGIPYWIVKNSWGPQWGMNGYFLIERGKNMCGLAACASYPIPLV";
		TS_ASSERT_EQUALS(aln2.seq(1)->sequence(), t);

		aln2.write("../../testsuit/data/out.fa", "fasta");
		aln.clear();
		aln.read("../../testsuit/data/out.fa");
		TS_ASSERT_EQUALS(aln.seq(1)->sequence(), t);
	}

	void test_MSF()
	{
	//	TS_TRACE("MSF format handling");
		Alignment aln;
		aln.read("../../testsuit/data/aln.msf");
		TS_ASSERT_EQUALS(aln.seq(0)->sequence(), "MAGAL-CALLLLQLLGRGEGKNEELRLYHYLFDTYDPGRRPVQEPEDTVTISLKVTLTNLISLNEKEETLTTSVWIGIDWQDYRLNYSKGDFGGVETLRV");

		Alignment aln2;
		vector<string> names;
		names.push_back("ACHE_BOVIN");
		names.push_back("ACHE_MOUSE");
		aln2.read("../../testsuit/data/aln.msf",names);
		TS_ASSERT_EQUALS(aln2.n_seqs(), 2);
		TS_ASSERT_EQUALS(aln2.seq(1)->name(), "ACHE_MOUSE");

		string t = "MAGALCALLLLQLLGRGEGKNEELRLYHYLFDTYDPGRRPVQEPEDTVTISLKVTLTNLISLNEKEETLTTSVWIGIDWQDYRLNYSKGDFGGVETLRV";
		TS_ASSERT_EQUALS(aln2.seq(0)->sequence(), t);

	//	aln2.write("./testsuit/data/out.msf", "msf");

	}

	void test_CLUSTALW()
	{
		Alignment aln;
		aln.read("../../testsuit/data/aln.aln");
		TS_ASSERT_EQUALS(aln.seq(0)->sequence(), "-----MKVILLFVLAVFTVFVSS---------------RGIPPEEQ------------SQFLEFQDKFNKKY-SHEEYLERFEIFKSNLGKIEELNLIAINHKADTKFGVNKFADLSSDEE-WQFYIGGVF-DIPCN--PNSLDHGILIVGYSAKNTIFRKNMPYWIVKNSWGADWGEQGYIYLRRGKNTCGVSNFVSTSII--");

		Alignment aln2;
		vector<string> names;
		names.push_back("CYS1_DICDI");
		names.push_back("CATH_HUMAN");
		aln2.read("../../testsuit/data/aln.aln",names);
		TS_ASSERT_EQUALS(aln2.n_seqs(), 2);
		TS_ASSERT_EQUALS(aln2.seq(0)->name(), "CYS1_DICDI");

		string t = "MKVILLFVLAVFTVFVSS-------RGIPPEEQSQFLEFQDKFNKKYSHEEYLERFEIFKSNLGKIEELNLIAINHKADTKFGVNKFADLSSDEE-WQFYIGGVF-DIPCN--PNSLDHGILIVGYSAKNTIFRKNMPYWIVKNSWGADWGEQGYIYLRRGKNTCGVSNFVSTSII--";
		TS_ASSERT_EQUALS(aln2.seq(0)->sequence(), t);

		aln2.write("../../testsuit/data/out.aln", "clustalw");
		aln.clear();
		aln.read("../../testsuit/data/out.aln");
		TS_ASSERT_EQUALS(aln.seq(0)->sequence(), t);
	}

	void test_STOCKHOLM()
	{
		Alignment aln;
		aln.read("../../testsuit/data/aln.stk");
		TS_ASSERT_EQUALS(aln.seq(8)->sequence(), "DIIENCKY--CGSFDIE---KVKDIYTCGDCTQTYTT");
		TS_ASSERT_EQUALS(aln.n_seqs(), 14);

		Alignment aln2;
		vector<string> names;
		names.push_back("Q8QNH7_ESV1/101-133");
		names.push_back("VLTF3_VACCC/1-32");
		names.push_back("YR429_MIMIV/213-247");
		aln2.read("../../testsuit/data/aln.stk",names);
		TS_ASSERT_EQUALS(aln2.n_seqs(), 3);
		TS_ASSERT_EQUALS(aln2.seq(1)->name(), "VLTF3_VACCC/1-32");
		string t = "MNLRLCSG--CRHNGIV-SEQGYEYCIFCESVFQK";
		TS_ASSERT_EQUALS(aln2.seq(1)->sequence(), t);
		aln.clear();
		aln.read("../../testsuit/data/aln2.stk");
		TS_ASSERT_EQUALS(aln.seq(0)->sequence(), "UUAAUCGAGCUCAACACUCUUCGUAUAUCCUC-UCAAUAUGG-GAUGAGGGUCUCUAC-AGGUA-CCGUAAA-UACCUAGCUACGAAAAGAAUGCAGUUAAUGU");
	}

	void test_PHYLIP()
	{
		Alignment aln;
		aln.read("../../testsuit/data/aln.phy_i");
		TS_ASSERT_EQUALS(aln.n_seqs(), 3);
		TS_ASSERT_EQUALS(aln.seq(1)->sequence(), "MAHARVLLLALAVLATAAVAVASSSSFADSNPIRPVTDRAASTLESAVLGALGRTRHALRFARFAVRYGKSYESAAEVRRRFRIFSESLEEVRSTN----RKGLPYRLGINRFSDMSWEEFQATRL-GAAQTCSATLAGNHLMRDA--AALPETKDWREDG-IVSPVKNQAHCGSCWTFSTTGALEAAYTQATGKNISLSEQQLVDCAGGFNNF--------GCNGGLPSQAFEYIKYNGGIDTEESYPYKGVNGV-CHYKAENAAVQVLDSVNITLNAEDELKNAVGLVRPVSVAFQVIDGFRQYKSGVYTSDHCGTTPDDVNHAVLAVGYGVENGV-----PYWLIKNSWGADWGDNGYFKMEMGKNMCAIATCASYPVVAA");

		aln.clear();
		vector<string> names;
		names.push_back("CYS1_DICDI");
		names.push_back("ALEU_HORVU");
		aln.read("../../testsuit/data/aln.phy_i",names);

		TS_ASSERT_EQUALS(aln.n_seqs(), 2);
		TS_ASSERT_EQUALS(aln.seq(1)->sequence(), "MAHARVLLLALAVLATAAVAVASSSSFADSNPIRPVTDRAASTLESAVLGALGRTRHALRFARFAVRYGKSYESAAEVRRRFRIFSESLEEVRSTN----RKGLPYRLGINRFSDMSWEEFQATRL-GAAQTCSATLAGNHLMRDA--AALPETKDWREDGIVSPVKNQAHCGSCWTFSTTGALEAAYTQATGKNISLSEQQLVDCAGGFNNF--------GCNGGLPSQAFEYIKYNGGIDTEESYPYKGVNGV-CHYKAENAAVQVLDSVNITLNAEDELKNAVGLVRPVSVAFQVIDGFRQYKSGVYTSDHCGTTPDDVNHAVLAVGYGVENGV-----PYWLIKNSWGADWGDNGYFKMEMGKNMCAIATCASYPVVAA");

		aln.clear();
		aln.read("../../testsuit/data/aln.phy_s");
		TS_ASSERT_EQUALS(aln.n_seqs(), 4);
		TS_ASSERT_EQUALS(aln.seq(2)->sequence(), "MKKLKKHPDFPKKPLTPYFRFFMEKRAKYAKLHP--EMSNLDLTKILSKKYKELPEKKKMKYIQDFQREKQEFERNLARFREDHPDLIQNAKK---");


		aln.clear();
		names.clear();
		names.push_back("2lef_A");
		names.push_back("1k99_A");
		aln.read("../../testsuit/data/aln.phy_s",names);
		TS_ASSERT_EQUALS(aln.n_seqs(), 2);
		TS_ASSERT_EQUALS(aln.seq(0)->sequence(), "MKKLKKHPDFPKKPLTPYFRFFMEKRAKYAKLHPEMSNLDLTKILSKKYKELPEKKKMKYIQDFQREKQEFERNLARFREDHPDLIQNAKK---");
	}

	void test_CODATA()
	{
		Alignment aln;
		aln.read("../../testsuit/data/aln.cd");
		TS_ASSERT_EQUALS(aln.seq(1)->sequence(), "TSPASIRPPAGPSSR---------RPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGW----------RASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE");
		TS_ASSERT_EQUALS(aln.n_seqs(), 4);

		Alignment aln2;
		vector<string> names;
		names.push_back("IXI_234");
		names.push_back("IXI_237");
		aln2.read("../../testsuit/data/aln.cd",names);
		TS_ASSERT_EQUALS(aln2.n_seqs(), 2);
		TS_ASSERT_EQUALS(aln2.seq(1)->name(), "IXI_237");

		string t = "TSPASLRPPAGPSSRPAMVSSRR-RPSPPGPRRPT----CSAAPRRPQATGGYKTCSGTCTTSTSTRHRGRSGYSARTTTAACLRASRKSMRAACSR--GSRPNRFAPTLMSSCLTSTTGPPAYAGDRSHE";
		TS_ASSERT_EQUALS(aln2.seq(1)->sequence(), t);
	}

	void test_AMPS()
	{
		Alignment aln;
		aln.read("../../testsuit/data/aln.ap");
		TS_ASSERT_EQUALS(aln.seq(2)->sequence(), "NAP--LV");
		TS_ASSERT_EQUALS(aln.n_seqs(), 6);

		Alignment aln2;
		vector<string> names;
		names.push_back("A0231");
		names.push_back("GLobin_C");
		aln2.read("../../testsuit/data/aln.ap",names);
		TS_ASSERT_EQUALS(aln2.n_seqs(), 2);
		TS_ASSERT_EQUALS(aln2.seq(1)->name(), "GLobin_C");

		string t = "QAPPWLL";
		TS_ASSERT_EQUALS(aln2.seq(1)->sequence(), t);
	}
};


class AlignmentModificationTest : public CxxTest::TestSuite
{
public:

	void test_col_trim()
	{
		printf("\nRunning alignment modification tests\n");
		Alignment aln;
		aln.read("../../testsuit/data/col_trim.fa");
		aln.column_trim(0.5);
		TS_ASSERT_EQUALS(aln[0].sequence(), "GAAAAAAAAC");
		aln.clear();
		aln.read("../../testsuit/data/col_trim.fa");
		aln.column_trim(0.75);
		TS_ASSERT_EQUALS(aln[2].sequence(), "AGAAAAAAA-AC");
		TS_ASSERT_EQUALS(aln[3].sequence(), "-GAAAAAAA-AC");
	}

	void test_merge()
	{
		Alignment aln1, aln2;
		aln1.read("../../testsuit/data/merge.fa");
		aln2.read("../../testsuit/data/col_trim.fa");
		vector<pair<unsigned int, unsigned int> > gap1, gap2;
		gap1.push_back(pair<unsigned int, unsigned int>(2,2));
		gap1.push_back(pair<unsigned int, unsigned int>(5,3));
		gap2.push_back(pair<unsigned int, unsigned int>(5,1));
		aln1.merge(gap1, aln2, gap2);
		TS_ASSERT_EQUALS(aln1[0].sequence(), "GG--GGG---GGGGG");
		TS_ASSERT_EQUALS(aln1[2].sequence(), "-GAAA-AAAAAAAAC");
	}

	void test_keep_seq()
	{
		Alignment aln;
		aln.read("../../testsuit/data/col_trim.fa");
		vector<size_t> ids;
		ids.push_back(1);
		ids.push_back(2);
		ids.push_back(3);
		aln.keep_seqs(ids);
		TS_ASSERT_EQUALS(aln.n_seqs(), 3);
		TS_ASSERT_EQUALS(aln[0].sequence(), "GGAAAAAAAAAC");
		TS_ASSERT_EQUALS(aln[1].sequence(), "AGAAAAAAA-AC");
	}


	void test_gapinsert()
	{
		Sequence seq("TEST", "", "ACGT");
		vector<pair<unsigned int, unsigned int> > vec;
		vec.push_back(pair<unsigned int, unsigned int>(0,3));
		vec.push_back(pair<unsigned int, unsigned int>(2,1));
		vec.push_back(pair<unsigned int, unsigned int>(4,2));
		seq.insert_gaps(vec);
		TS_ASSERT_EQUALS(seq.sequence(), "---AC-GT--");
	}

};




class AlignmentAnalysisTest : public CxxTest::TestSuite
{
public:
	void test_scoring()
	{
		printf("\nRunning alignment analysis tests\n");
		Scoring_Matrix matrix;
		matrix.read("../../data/BLOSUM62");
		Alignment aln;
		aln.read("../../testsuit/data/score_test.fa");
		TS_ASSERT_EQUALS(sum_of_pairs_score(aln, matrix, -2, -1),80.0);
	}

	void test_identity()
	{
		Alignment aln;
		aln.read("../../testsuit/data/score_test.fa");
		pair<size_t, size_t> id = identity(aln);
		TS_ASSERT_EQUALS(id.first, 15);
		TS_ASSERT_EQUALS(id.second, 17);
	}

	void test_consensus()
	{
		Alignment aln;
		aln.read("../../testsuit/data/consensus.fa");
		string result=simple_consensus(aln, 0.51, 'N')->sequence();
		TS_ASSERT_EQUALS(result, "AAAAACCCNN");
	}

};


