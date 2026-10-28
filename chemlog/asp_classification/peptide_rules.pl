% ASP (clingo) rules for the ChemLog peptide classification.
%
% This program mirrors the algorithmic classifiers in chemlog/alg_classification (and the MSOL / FOL
% definitions they implement) and is expected to produce the same results:
%   - number of amino acid residues connected via amide bonds (peptide_size/2)
%   - charge category (charge_category/2)
%   - proteinogenic amino acid residues (proteinogenic_residue/5)
%   - emericellamides and 2,5-diketopiperazines
%   - the resulting ChEBI classes (chebi_class/2)
%
% Input facts are produced by chemlog.preprocessing.mol_to_asp (see there for the full list). Each
% molecule M is a fact molecule(M), its atoms are linked via has_atom(M, A); all atom-level
% predicates (c/1, n/1, bSINGLE/2, has_bond_to/2, num_hs/2, formal_charge/2, cip_code_S/1, ...)
% range over atoms, global properties (net_charge_neutral/1, ...) over molecules. Several molecules
% can be classified in one program, since no rule connects atoms of different molecules.
%
% The peptide size requires a choice (which carbon fragment an amino group belongs to) and is
% determined by optimisation: only the optimal answer set has the correct peptide_size/2 (and the
% classes depending on it).

% input predicates that do not occur in every molecule
#defined c/1. #defined n/1. #defined o/1. #defined s/1. #defined se/1. #defined h/1.
#defined aromatic/1. #defined cip_code_S/1. #defined cip_code_R/1.
#defined has_bond_to/2. #defined bSINGLE/2. #defined bDOUBLE/2. #defined bAROMATIC/2.
#defined net_charge_positive/1. #defined net_charge_negative/1. #defined net_charge_neutral/1.

% ---------------------------------------------------------------------------------------------------
% General helpers
% ---------------------------------------------------------------------------------------------------

% total degree (SMARTS X): bonds to heavy atoms (explicit H atoms are already counted in num_hs/2) + Hs
xdeg(A,H+K) :- num_hs(A,H), K = #count{ B : has_bond_to(A,B), not h(B) }.
% unspecified SMARTS bond: single or aromatic
sb(A,B) :- bSINGLE(A,B).
sb(A,B) :- bAROMATIC(A,B).

% ---------------------------------------------------------------------------------------------------
% Functional groups (chemlog.alg_classification.peptide_size_classifier)
% ---------------------------------------------------------------------------------------------------

% X-O where O has a single bond to X and either one H atom or a negative charge (imidic form of an amide)
single_o(X,O) :- bSINGLE(X,O), o(O), num_hs(O,1).
single_o(X,O) :- bSINGLE(X,O), o(O), formal_charge(O,-1).
double_o(X,O) :- bDOUBLE(X,O), o(O).

% the nitrogen of an amide bond may only be connected to C and H atoms, using single bonds (or a double
% bond in the imidic form N=C-O)
amide_n_violation(N) :- n(N), has_bond_to(N,X), not c(X), not h(X).
amide_n_violation(N) :- n(N), has_bond_to(N,X), not bSINGLE(N,X), not imidic_bond(N,X).
imidic_bond(N,X) :- n(N), bDOUBLE(N,X), single_o(X,_).

% amide bond C(=O)N or imidic form C(O)=N
amide_bond(C,O,N) :- n(N), c(C), bSINGLE(C,N), not amide_n_violation(N), double_o(C,O).
amide_bond(C,O,N) :- n(N), c(C), bDOUBLE(C,N), not amide_n_violation(N), single_o(C,O).
amide_c(C) :- amide_bond(C,_,_).

% amino group: N with only C / H neighbours, single bonds (or a double bond to the C of an amide bond)
amino_violation(N) :- n(N), has_bond_to(N,X), not c(X), not h(X).
amino_violation(N) :- n(N), has_bond_to(N,X), not bSINGLE(N,X), not amide_c(X).
amino(N) :- n(N), not amino_violation(N).

% carboxylic acid derivative C(=O)Y with Y not C or H, represented by its carbon
carboxy(C) :- c(C), double_o(C,_), bSINGLE(C,Y), not c(Y), not h(Y).

% ---------------------------------------------------------------------------------------------------
% Peptide size (chemlog.alg_classification.peptide_size_classifier.get_n_amino_acid_residues)
% ---------------------------------------------------------------------------------------------------

% carbon fragments ("chunks"): connected components of carbon atoms that have at least one C-C bond.
% Each chunk is represented by its smallest atom.
cc_bond(A,B) :- has_bond_to(A,B), c(A), c(B).
c_reach(A,A) :- cc_bond(A,_).
c_reach(A,B) :- c_reach(A,X), cc_bond(X,B).
not_chunk_rep(A) :- c_reach(A,B), B < A.
in_chunk(A,R) :- c_reach(A,R), not not_chunk_rep(R).
chunk_has_carboxy(R) :- in_chunk(C,R), carboxy(C).

% an amino group can belong to each chunk that contains a carboxy group and one of its neighbours
% (except for the neighbour it forms an amide bond with)
amino_chunk_candidate(N,R) :- amino(N), has_bond_to(N,X), not amide_bond(X,_,N), in_chunk(X,R),
    chunk_has_carboxy(R).
% choose one assignment for each amino group (amino groups without candidates remain unassigned)
{ amino_in_chunk(N,R) : amino_chunk_candidate(N,R) } = 1 :- amino_chunk_candidate(N,_).

% amino acid residue: chunk with a carboxy group and at least one amino group assigned to it
aar(R) :- amino_in_chunk(_,R).
% two residues are linked if the C of an amide bond belongs to one and its N (amino group) to the other
aar_link(RC,RN) :- amide_bond(C,_,N), amino_in_chunk(N,RN), in_chunk(C,RC), aar(RC).
aar_link(RN,RC) :- aar_link(RC,RN).

% a molecule only counts as peptide candidate if it has an amide bond and at least two residues
has_amide_bond(M) :- molecule(M), has_atom(M,C), amide_c(C).
peptide_candidate(M) :- has_amide_bond(M), 2 <= #count{ R : aar(R), has_atom(M,R) }.

% the longest chain is the connected component (of the aar_link graph) around a chosen root residue
{ chain_root(R) : aar(R), has_atom(M,R) } = 1 :- peptide_candidate(M).
in_chain(R) :- chain_root(R).
in_chain(R2) :- in_chain(R1), aar_link(R1,R2).
#maximize { 1,R : in_chain(R) }.

peptide_size(M,K) :- molecule(M), K = #count{ R : in_chain(R), has_atom(M,R) }.

% ---------------------------------------------------------------------------------------------------
% Charge category (chemlog.alg_classification.charge_classifier)
% ---------------------------------------------------------------------------------------------------

charged(A) :- formal_charge(A,Q), Q != 0.
% fragments (connected components) are only relevant if they contain charged atoms -> only compute
% reachability from charged atoms; each fragment is represented by its smallest charged atom
frag_reach(Q,Q) :- charged(Q).
frag_reach(Q,B) :- frag_reach(Q,A), has_bond_to(A,B).
not_frag_rep(Q) :- frag_reach(Q,Q2), charged(Q2), Q2 < Q.
frag_rep(Q) :- charged(Q), not not_frag_rep(Q).
frag_charge(Q,S) :- frag_rep(Q), S = #sum{ C,A : frag_reach(Q,A), formal_charge(A,C) }.

ion_neighbour(A) :- has_bond_to(A,B), charged(B).
isolated_anion(A) :- formal_charge(A,C), C < 0, not ion_neighbour(A).
isolated_cation(A) :- formal_charge(A,C), C > 0, not ion_neighbour(A).

anionic_fragment(M) :- frag_rep(Q), has_atom(M,Q), frag_charge(Q,S), S < 0.
cationic_fragment(M) :- frag_rep(Q), has_atom(M,Q), frag_charge(Q,S), S > 0.
% neutral fragment with a positive and a negative charge on atoms that are not adjacent to other charges
zwitterionic_fragment(M) :- frag_rep(Q), has_atom(M,Q), frag_charge(Q,0),
    frag_reach(Q,A), isolated_anion(A), frag_reach(Q,C), isolated_cation(C).

salt(M) :- anionic_fragment(M), cationic_fragment(M).
zwitterion(M) :- zwitterionic_fragment(M), net_charge_neutral(M), not salt(M).
anion(M) :- anionic_fragment(M), not salt(M), not zwitterion(M).
cation(M) :- cationic_fragment(M), not salt(M), not zwitterion(M), not anion(M).
neutral(M) :- molecule(M), not salt(M), not zwitterion(M), not anion(M), not cation(M).

charge_category(M,"SALT") :- salt(M).
charge_category(M,"ZWITTERION") :- zwitterion(M).
charge_category(M,"ANION") :- anion(M).
charge_category(M,"CATION") :- cation(M).
charge_category(M,"NEUTRAL") :- neutral(M).

% ---------------------------------------------------------------------------------------------------
% Proteinogenic amino acids (chemlog.alg_classification.proteinogenics_classifier)
% ---------------------------------------------------------------------------------------------------

% alpha amino acid backbone: amino group N, alpha carbon (only bonded to carbons and N) and carboxy C
alpha_hetero_neighbour(N,Al) :- amino(N), has_bond_to(N,Al), c(Al), has_bond_to(Al,X), X != N, not c(X).
backbone(N,Al,Ca) :- amino(N), has_bond_to(N,Al), c(Al), not alpha_hetero_neighbour(N,Al),
    has_bond_to(Al,Ca), carboxy(Ca).

% side chains start at a carbon neighbour of the alpha carbon (not part of the backbone)
side_chain_start_candidate(N,Al,Ca,S) :- backbone(N,Al,Ca), has_bond_to(Al,S), c(S), S != Ca,
    not has_bond_to(Ca,S).
has_side_chain_start(N,Al,Ca) :- side_chain_start_candidate(N,Al,Ca,_).
side_chain_start(N,Al,Ca,S) :- side_chain_start_candidate(N,Al,Ca,S), num_hs(Al,1),
    #count{ S2 : side_chain_start_candidate(N,Al,Ca,S2) } = 1.

% side chains, translated from the SMARTS patterns in proteinogenics_classifier.side_chains
% (sm_* are the SMARTS atom primitives, defined at the end of this section)
% W: [CH2X4][CX3]1=[CX3H][NX3H+0][CX3]2=[CX3]1[CX3H]=[CX3H][CX3H]=[CX3H]2
side_chain(N,Al,Ca,"W") :- side_chain_start(N,Al,Ca,S0), sm_c_h2_x4(S0), sm_c_x3(S1), sm_c_h1_x3(S2),
    sm_n_h1_x3_q0(S3), sm_c_x3(S4), sm_c_x3(S5), sm_c_h1_x3(S6), sm_c_h1_x3(S7), sm_c_h1_x3(S8),
    sm_c_h1_x3(S9), sb(S0,S1), bDOUBLE(S1,S2), sb(S2,S3), sb(S3,S4), bDOUBLE(S4,S5), sb(S5,S6),
    bDOUBLE(S6,S7), sb(S7,S8), bDOUBLE(S8,S9), sb(S5,S1), sb(S9,S4), S1!=S4, S1!=S6, S1!=S7, S1!=S8, S1!=S9,
    S2!=S4, S2!=S5, S2!=S6, S2!=S7, S2!=S8, S2!=S9, S4!=S6, S4!=S7, S4!=S8, S5!=S7, S5!=S8, S5!=S9, S6!=S8,
    S6!=S9, S7!=S9, S1!=Al, S2!=Al, S4!=Al, S5!=Al, S6!=Al, S7!=Al, S8!=Al, S9!=Al.

% Y: [CH2X4][CX3]1=[CX3H][CX3H]=[CX3]([OX2H+0])[CX3H]=[CX3H]1
side_chain(N,Al,Ca,"Y") :- side_chain_start(N,Al,Ca,S0), sm_c_h2_x4(S0), sm_c_x3(S1), sm_c_h1_x3(S2),
    sm_c_h1_x3(S3), sm_c_x3(S4), sm_o_h1_x2_q0(S5), sm_c_h1_x3(S6), sm_c_h1_x3(S7), sb(S0,S1), bDOUBLE(S1,S2),
    sb(S2,S3), bDOUBLE(S3,S4), sb(S4,S5), sb(S4,S6), bDOUBLE(S6,S7), sb(S7,S1), S1!=S3, S1!=S4, S1!=S6,
    S2!=S4, S2!=S6, S2!=S7, S3!=S6, S3!=S7, S4!=S7, S1!=Al, S2!=Al, S3!=Al, S4!=Al, S6!=Al, S7!=Al.

% F: [CH2X4][CX3]1=[CX3H][CX3H]=[CX3H][CX3H]=[CX3H]1
side_chain(N,Al,Ca,"F") :- side_chain_start(N,Al,Ca,S0), sm_c_h2_x4(S0), sm_c_x3(S1), sm_c_h1_x3(S2),
    sm_c_h1_x3(S3), sm_c_h1_x3(S4), sm_c_h1_x3(S5), sm_c_h1_x3(S6), sb(S0,S1), bDOUBLE(S1,S2), sb(S2,S3),
    bDOUBLE(S3,S4), sb(S4,S5), bDOUBLE(S5,S6), sb(S6,S1), S1!=S3, S1!=S4, S1!=S5, S2!=S4, S2!=S5, S2!=S6,
    S3!=S5, S3!=S6, S4!=S6, S1!=Al, S2!=Al, S3!=Al, S4!=Al, S5!=Al, S6!=Al.

% H: [CH2X4][CX3]1=[CX3H][NHX3][CX3H]=[NX2]1
side_chain(N,Al,Ca,"H") :- side_chain_start(N,Al,Ca,S0), sm_c_h2_x4(S0), sm_c_x3(S1), sm_c_h1_x3(S2),
    sm_n_h1_x3(S3), sm_c_h1_x3(S4), sm_n_x2(S5), sb(S0,S1), bDOUBLE(S1,S2), sb(S2,S3), sb(S3,S4),
    bDOUBLE(S4,S5), sb(S5,S1), S1!=S4, S2!=S4, S1!=Al, S2!=Al, S4!=Al.

% H_iso: [CH2X4][CX3]1=[CX3H][NX2]=[CX3H][NHX3]1
side_chain(N,Al,Ca,"H") :- side_chain_start(N,Al,Ca,S0), sm_c_h2_x4(S0), sm_c_x3(S1), sm_c_h1_x3(S2),
    sm_n_x2(S3), sm_c_h1_x3(S4), sm_n_h1_x3(S5), sb(S0,S1), bDOUBLE(S1,S2), sb(S2,S3), bDOUBLE(S3,S4),
    sb(S4,S5), sb(S5,S1), S1!=S4, S2!=S4, S1!=Al, S2!=Al, S4!=Al.

% R: [CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])N
side_chain(N,Al,Ca,"R") :- side_chain_start(N,Al,Ca,S0), sm_c_h2_x4(S0), sm_c_h2_x4(S1), sm_c_h2_x4(S2),
    sm_n_h1_x3(S3), sm_c_h0_x3(S4), sm_n_h2_x3_q1__n_h1_x2_q0(S5), sm_n(S6), sb(S0,S1), sb(S1,S2), sb(S2,S3),
    sb(S3,S4), bDOUBLE(S4,S5), sb(S4,S6), amino(S3), amino(S6), S0!=S2, S3!=S6, S5!=S6, S1!=Al, S2!=Al,
    S4!=Al.

% R_iso: [CH2X4][CH2X4][CH2X4][NX2]=[CH0X3](N)N
side_chain(N,Al,Ca,"R") :- side_chain_start(N,Al,Ca,S0), sm_c_h2_x4(S0), sm_c_h2_x4(S1), sm_c_h2_x4(S2),
    sm_n_x2(S3), sm_c_h0_x3(S4), sm_n(S5), sm_n(S6), sb(S0,S1), sb(S1,S2), sb(S2,S3), bDOUBLE(S3,S4),
    sb(S4,S5), sb(S4,S6), amino(S5), amino(S6), S0!=S2, S3!=S5, S3!=S6, S5!=S6, S1!=Al, S2!=Al, S4!=Al.

% O: [CH2X4][CH2X4][CH2X4][CH2X4]N[CX3][CH1X4]1[NX2]=[CH1X3][CH2X4][CH1X4]([CH3X4])1
side_chain(N,Al,Ca,"O") :- side_chain_start(N,Al,Ca,S0), sm_c_h2_x4(S0), sm_c_h2_x4(S1), sm_c_h2_x4(S2),
    sm_c_h2_x4(S3), sm_n(S4), sm_c_x3(S5), sm_c_h1_x4(S6), sm_n_x2(S7), sm_c_h1_x3(S8), sm_c_h2_x4(S9),
    sm_c_h1_x4(S10), sm_c_h3_x4(S11), sb(S0,S1), sb(S1,S2), sb(S2,S3), sb(S3,S4), sb(S4,S5), sb(S5,S6),
    sb(S6,S7), bDOUBLE(S7,S8), sb(S8,S9), sb(S9,S10), sb(S10,S11), sb(S10,S6), amino(S4), carboxy(S5), S0!=S2,
    S0!=S3, S0!=S9, S1!=S3, S1!=S9, S2!=S9, S3!=S9, S4!=S7, S5!=S8, S1!=Al, S2!=Al, S3!=Al, S5!=Al, S6!=Al,
    S8!=Al, S9!=Al, S10!=Al, S11!=Al.

% K: [CH2X4][CH2X4][CH2X4][CH2X4]N
side_chain(N,Al,Ca,"K") :- side_chain_start(N,Al,Ca,S0), sm_c_h2_x4(S0), sm_c_h2_x4(S1), sm_c_h2_x4(S2),
    sm_c_h2_x4(S3), sm_n(S4), sb(S0,S1), sb(S1,S2), sb(S2,S3), sb(S3,S4), amino(S4), S0!=S2, S0!=S3, S1!=S3,
    S1!=Al, S2!=Al, S3!=Al.

% M: [CH2X4][CH2X4][SX2][CH3X4]
side_chain(N,Al,Ca,"M") :- side_chain_start(N,Al,Ca,S0), sm_c_h2_x4(S0), sm_c_h2_x4(S1), sm_s_x2(S2),
    sm_c_h3_x4(S3), sb(S0,S1), sb(S1,S2), sb(S2,S3), S1!=Al, S3!=Al.

% Q: [CH2X4][CH2X4][CX3](=[OX1+0])[NH2X3]
side_chain(N,Al,Ca,"Q") :- side_chain_start(N,Al,Ca,S0), sm_c_h2_x4(S0), sm_c_h2_x4(S1), sm_c_x3(S2),
    sm_o_x1_q0(S3), sm_n_h2_x3(S4), sb(S0,S1), sb(S1,S2), bDOUBLE(S2,S3), sb(S2,S4), amino(S4), S1!=Al,
    S2!=Al.

% E: [CH2X4][CH2X4][CX3]
side_chain(N,Al,Ca,"E") :- side_chain_start(N,Al,Ca,S0), sm_c_h2_x4(S0), sm_c_h2_x4(S1), sm_c_x3(S2),
    sb(S0,S1), sb(S1,S2), carboxy(S2), S1!=Al, S2!=Al.

% P: [CH2X4][CH2X4][CH2X4]N
side_chain(N,Al,Ca,"P") :- side_chain_start(N,Al,Ca,S0), sm_c_h2_x4(S0), sm_c_h2_x4(S1), sm_c_h2_x4(S2),
    sm_n(S3), sb(S0,S1), sb(S1,S2), sb(S2,S3), amino(S3), S0!=S2, S1!=Al, S2!=Al,
    S3=N.  % proline: the ring closes at the backbone N

% N: [CH2X4][CX3](=[OX1+0])[NH2X3]
side_chain(N,Al,Ca,"N") :- side_chain_start(N,Al,Ca,S0), sm_c_h2_x4(S0), sm_c_x3(S1), sm_o_x1_q0(S2),
    sm_n_h2_x3(S3), sb(S0,S1), bDOUBLE(S1,S2), sb(S1,S3), amino(S3), S1!=Al.

% D: [CH2X4][CX3]
side_chain(N,Al,Ca,"D") :- side_chain_start(N,Al,Ca,S0), sm_c_h2_x4(S0), sm_c_x3(S1), sb(S0,S1), carboxy(S1),
    S1!=Al.

% I: [CHX4]([CH3X4])[CH2X4][CH3X4]
side_chain(N,Al,Ca,"I") :- side_chain_start(N,Al,Ca,S0), sm_c_h1_x4(S0), sm_c_h3_x4(S1), sm_c_h2_x4(S2),
    sm_c_h3_x4(S3), sb(S0,S1), sb(S0,S2), sb(S2,S3), S1!=S3, S1!=Al, S2!=Al, S3!=Al.

% L: [CH2X4][CHX4]([CH3X4])[CH3X4]
side_chain(N,Al,Ca,"L") :- side_chain_start(N,Al,Ca,S0), sm_c_h2_x4(S0), sm_c_h1_x4(S1), sm_c_h3_x4(S2),
    sm_c_h3_x4(S3), sb(S0,S1), sb(S1,S2), sb(S1,S3), S2!=S3, S1!=Al, S2!=Al, S3!=Al.

% C: [CH2X4][SHX2+0,SX1-]
side_chain(N,Al,Ca,"C") :- side_chain_start(N,Al,Ca,S0), sm_c_h2_x4(S0), sm_s_h1_x2_q0__s_x1_qm1(S1),
    sb(S0,S1).

% U: [CH2X4][SeHX2+0,SeX1-]
side_chain(N,Al,Ca,"U") :- side_chain_start(N,Al,Ca,S0), sm_c_h2_x4(S0), sm_se_h1_x2_q0__se_x1_qm1(S1),
    sb(S0,S1).

% T: [CHX4]([CH3X4])[OH1X2]
side_chain(N,Al,Ca,"T") :- side_chain_start(N,Al,Ca,S0), sm_c_h1_x4(S0), sm_c_h3_x4(S1), sm_o_h1_x2(S2),
    sb(S0,S1), sb(S0,S2), S1!=Al.

% V: [CHX4]([CH3X4])[CH3X4]
side_chain(N,Al,Ca,"V") :- side_chain_start(N,Al,Ca,S0), sm_c_h1_x4(S0), sm_c_h3_x4(S1), sm_c_h3_x4(S2),
    sb(S0,S1), sb(S0,S2), S1!=S2, S1!=Al, S2!=Al.

% S: [CH2X4][OH1X2]
side_chain(N,Al,Ca,"S") :- side_chain_start(N,Al,Ca,S0), sm_c_h2_x4(S0), sm_o_h1_x2(S1), sb(S0,S1).

% A: [CH3X4]
side_chain(N,Al,Ca,"A") :- side_chain_start(N,Al,Ca,S0), sm_c_h3_x4(S0).

% the first side chain (in this order) that matches is used
side_chain_rank("W",1). side_chain_rank("Y",2). side_chain_rank("F",3). side_chain_rank("H",4).
side_chain_rank("R",5). side_chain_rank("O",6). side_chain_rank("K",7). side_chain_rank("M",8).
side_chain_rank("Q",9). side_chain_rank("E",10). side_chain_rank("P",11). side_chain_rank("N",12).
side_chain_rank("D",13). side_chain_rank("I",14). side_chain_rank("L",15). side_chain_rank("C",16).
side_chain_rank("U",17). side_chain_rank("T",18). side_chain_rank("V",19). side_chain_rank("S",20).
side_chain_rank("A",21).
side_chain_preceded(N,Al,Ca,I) :- side_chain(N,Al,Ca,L), side_chain_rank(L,I),
    side_chain(N,Al,Ca,L2), side_chain_rank(L2,J), J < I.
first_side_chain(N,Al,Ca,L) :- side_chain(N,Al,Ca,L), side_chain_rank(L,I), not side_chain_preceded(N,Al,Ca,I).

% L-amino acids: CIP code S, except for cysteine and selenocysteine (R)
expected_cip("C","R"). expected_cip("U","R").
expected_cip(L,"S") :- side_chain_rank(L,_), L != "C", L != "U".
cip(A,"S") :- cip_code_S(A).
cip(A,"R") :- cip_code_R(A).

% N-formyl group (SMARTS N[CH1X3]=[OX1])
formylated(N) :- sm_n(N), sb(N,F), sm_c_h1_x3(F), bDOUBLE(F,O), sm_o_x1(O).

proteinogenic(N,Al,Ca,"G") :- backbone(N,Al,Ca), num_hs(Al,2), not has_side_chain_start(N,Al,Ca).
proteinogenic(N,Al,Ca,L) :- first_side_chain(N,Al,Ca,L), L != "M", expected_cip(L,X), cip(Al,X).
proteinogenic(N,Al,Ca,"M") :- first_side_chain(N,Al,Ca,"M"), not formylated(N), cip(Al,"S").
proteinogenic(N,Al,Ca,"fMet") :- first_side_chain(N,Al,Ca,"M"), formylated(N), cip(Al,"S").
proteinogenic_residue(M,L,N,Al,Ca) :- proteinogenic(N,Al,Ca,L), has_atom(M,N).

% SMARTS atom primitives: element (upper case = aliphatic), H = total number of Hs, X = total degree
sm_o_x1(A) :- o(A), not aromatic(A), xdeg(A,1).
sm_c(A) :- c(A), not aromatic(A).
sm_c_h0_x3(A) :- c(A), not aromatic(A), num_hs(A,0), xdeg(A,3).
sm_c_h1__c_h2(A) :- c(A), not aromatic(A), num_hs(A,1).
sm_c_h1__c_h2(A) :- c(A), not aromatic(A), num_hs(A,2).
sm_c_h1_x3(A) :- c(A), not aromatic(A), num_hs(A,1), xdeg(A,3).
sm_c_h1_x4(A) :- c(A), not aromatic(A), num_hs(A,1), xdeg(A,4).
sm_c_h2_x4(A) :- c(A), not aromatic(A), num_hs(A,2), xdeg(A,4).
sm_c_h3_x4(A) :- c(A), not aromatic(A), num_hs(A,3), xdeg(A,4).
sm_c_x3(A) :- c(A), not aromatic(A), xdeg(A,3).
sm_n(A) :- n(A), not aromatic(A).
sm_n_h1_x3(A) :- n(A), not aromatic(A), num_hs(A,1), xdeg(A,3).
sm_n_h1_x3_q0(A) :- n(A), not aromatic(A), num_hs(A,1), xdeg(A,3), formal_charge(A,0).
sm_n_h2_x3(A) :- n(A), not aromatic(A), num_hs(A,2), xdeg(A,3).
sm_n_h2_x3_q1__n_h1_x2_q0(A) :- n(A), not aromatic(A), num_hs(A,2), xdeg(A,3), formal_charge(A,1).
sm_n_h2_x3_q1__n_h1_x2_q0(A) :- n(A), not aromatic(A), num_hs(A,1), xdeg(A,2), formal_charge(A,0).
sm_n_x2(A) :- n(A), not aromatic(A), xdeg(A,2).
sm_o_h0_x1(A) :- o(A), not aromatic(A), num_hs(A,0), xdeg(A,1).
sm_o_h0_x2(A) :- o(A), not aromatic(A), num_hs(A,0), xdeg(A,2).
sm_o_h1_x2(A) :- o(A), not aromatic(A), num_hs(A,1), xdeg(A,2).
sm_o_h1_x2_q0(A) :- o(A), not aromatic(A), num_hs(A,1), xdeg(A,2), formal_charge(A,0).
sm_o_x1_q0(A) :- o(A), not aromatic(A), xdeg(A,1), formal_charge(A,0).
sm_s_h1_x2_q0__s_x1_qm1(A) :- s(A), not aromatic(A), num_hs(A,1), xdeg(A,2), formal_charge(A,0).
sm_s_h1_x2_q0__s_x1_qm1(A) :- s(A), not aromatic(A), xdeg(A,1), formal_charge(A,-1).
sm_s_x2(A) :- s(A), not aromatic(A), xdeg(A,2).
sm_se_h1_x2_q0__se_x1_qm1(A) :- se(A), num_hs(A,1), xdeg(A,2), formal_charge(A,0).
sm_se_h1_x2_q0__se_x1_qm1(A) :- se(A), xdeg(A,1), formal_charge(A,-1).

% ---------------------------------------------------------------------------------------------------
% Substructures (chemlog.alg_classification.substructure_classifier), translated from SMARTS
% ---------------------------------------------------------------------------------------------------

% 2,5-diketopiperazine: C1[NX3H][CX3H0](=[OX1H0])C[NX3H][CX3H0]1=[OX1H0]
diketopiperazine_match(M) :- has_atom(M,X0), sm_c(X0), sm_n_h1_x3(X1), sm_c_h0_x3(X2), sm_o_h0_x1(X3),
    sm_c(X4), sm_n_h1_x3(X5), sm_c_h0_x3(X6), sm_o_h0_x1(X7), sb(X0,X1), sb(X1,X2), bDOUBLE(X2,X3), sb(X2,X4),
    sb(X4,X5), sb(X5,X6), bDOUBLE(X6,X7), sb(X6,X0), X0!=X2, X0!=X4, X1!=X5, X2!=X6, X3!=X7, X4!=X6.

% emericellamide (stereo information in the SMARTS is ignored, as in RDKit's default substructure match)
emericellamide_match(M) :- has_atom(M,X0), sm_c_h1__c_h2(X0), sm_o_h0_x2(X1), sm_c_h0_x3(X2), sm_o_h0_x1(X3),
    sm_c_h1_x4(X4), sm_c_h3_x4(X5), sm_n_h1_x3(X6), sm_c_h0_x3(X7), sm_o_h0_x1(X8), sm_c_h1_x4(X9),
    sm_c_h3_x4(X10), sm_n_h1_x3(X11), sm_c_h0_x3(X12), sm_o_h0_x1(X13), sm_c_h1_x4(X14), sm_c_h2_x4(X15),
    sm_c_h1_x4(X16), sm_c_h3_x4(X17), sm_c_h3_x4(X18), sm_n_h1_x3(X19), sm_c_h0_x3(X20), sm_o_h0_x1(X21),
    sm_c_h1_x4(X22), sm_n_h1_x3(X23), sm_c_h0_x3(X24), sm_o_h0_x1(X25), sm_c_h2_x4(X26), sm_n_h1_x3(X27),
    sm_c_h0_x3(X28), sm_o_h0_x1(X29), sm_c_h1__c_h2(X30), sm_c_h1_x4(X31), sm_c_h3_x4(X32), sm_c_h3_x4(X33),
    sb(X0,X1), sb(X1,X2), bDOUBLE(X2,X3), sb(X2,X4), sb(X4,X5), sb(X4,X6), sb(X6,X7), bDOUBLE(X7,X8),
    sb(X7,X9), sb(X9,X10), sb(X9,X11), sb(X11,X12), bDOUBLE(X12,X13), sb(X12,X14), sb(X14,X15), sb(X15,X16),
    sb(X16,X17), sb(X16,X18), sb(X14,X19), sb(X19,X20), bDOUBLE(X20,X21), sb(X20,X22), sb(X22,X23),
    sb(X23,X24), bDOUBLE(X24,X25), sb(X24,X26), sb(X26,X27), sb(X27,X28), bDOUBLE(X28,X29), sb(X28,X30),
    sb(X22,X31), sb(X31,X32), sb(X31,X33), sb(X30,X0), X0!=X4, X0!=X9, X0!=X14, X0!=X15, X0!=X16, X0!=X22,
    X0!=X26, X0!=X31, X2!=X7, X2!=X12, X2!=X20, X2!=X24, X2!=X28, X3!=X8, X3!=X13, X3!=X21, X3!=X25, X3!=X29,
    X4!=X9, X4!=X14, X4!=X16, X4!=X22, X4!=X30, X4!=X31, X5!=X10, X5!=X17, X5!=X18, X5!=X32, X5!=X33, X6!=X11,
    X6!=X19, X6!=X23, X6!=X27, X7!=X12, X7!=X20, X7!=X24, X7!=X28, X8!=X13, X8!=X21, X8!=X25, X8!=X29,
    X9!=X14, X9!=X16, X9!=X22, X9!=X30, X9!=X31, X10!=X17, X10!=X18, X10!=X32, X10!=X33, X11!=X19, X11!=X23,
    X11!=X27, X12!=X20, X12!=X24, X12!=X28, X13!=X21, X13!=X25, X13!=X29, X14!=X16, X14!=X22, X14!=X30,
    X14!=X31, X15!=X26, X15!=X30, X16!=X22, X16!=X30, X16!=X31, X17!=X18, X17!=X32, X17!=X33, X18!=X32,
    X18!=X33, X19!=X23, X19!=X27, X20!=X24, X20!=X28, X21!=X25, X21!=X29, X22!=X30, X23!=X27, X24!=X28,
    X25!=X29, X26!=X30, X30!=X31, X32!=X33.

% substructures are only checked for molecules of the matching size
emericellamide(M) :- peptide_size(M,5), emericellamide_match(M).
diketopiperazine(M) :- peptide_size(M,2), diketopiperazine_match(M).

% ---------------------------------------------------------------------------------------------------
% ChEBI classes (chemlog.cli.resolve_chebi_classes)
% ---------------------------------------------------------------------------------------------------

peptide_like(M) :- peptide_size(M,K), K >= 2.

peptide_anion(M) :- peptide_like(M), anion(M).
peptide_cation(M) :- peptide_like(M), cation(M).
peptide_zwitterion(M) :- peptide_like(M), zwitterion(M).
dipeptide_zwitterion(M) :- peptide_zwitterion(M), peptide_size(M,2).
tripeptide_zwitterion(M) :- peptide_zwitterion(M), peptide_size(M,3).
% neutral peptides
peptide(M) :- peptide_like(M), neutral(M).
dipeptide(M) :- peptide(M), peptide_size(M,2).
tripeptide(M) :- peptide(M), peptide_size(M,3).
tetrapeptide(M) :- peptide(M), peptide_size(M,4).
pentapeptide(M) :- peptide(M), peptide_size(M,5).
oligopeptide(M) :- peptide(M), peptide_size(M,K), K < 10.
polypeptide(M) :- peptide(M), peptide_size(M,K), K >= 10.

chebi_class(M,"24866") :- salt(M).
chebi_class(M,"22563") :- anion(M).
chebi_class(M,"36961") :- cation(M).
chebi_class(M,"27369") :- zwitterion(M).
chebi_class(M,"60334") :- peptide_anion(M).
chebi_class(M,"60194") :- peptide_cation(M).
chebi_class(M,"60466") :- peptide_zwitterion(M).
chebi_class(M,"90799") :- dipeptide_zwitterion(M).
chebi_class(M,"155837") :- tripeptide_zwitterion(M).
chebi_class(M,"16670") :- peptide(M).
chebi_class(M,"46761") :- dipeptide(M).
chebi_class(M,"47923") :- tripeptide(M).
chebi_class(M,"48030") :- tetrapeptide(M).
chebi_class(M,"48545") :- pentapeptide(M).
chebi_class(M,"25676") :- oligopeptide(M).
chebi_class(M,"15841") :- polypeptide(M).
chebi_class(M,"64372") :- emericellamide(M).
chebi_class(M,"65061") :- diketopiperazine(M).

#show peptide_size/2.
#show charge_category/2.
#show proteinogenic_residue/5.
#show emericellamide_match/1.
#show diketopiperazine_match/1.
#show chebi_class/2.
