import re

from rdkit import Chem

from chemlog.preprocessing.mol_to_fol import mol_to_fol_atoms

# predicates of the FOL structure that are not needed as ASP facts (equality / domain predicates)
_SKIPPED_PREDICATES = {"EQ", "atom", "global"}
_ASP_IDENTIFIER = re.compile(r"^_*[a-z][A-Za-z0-9_']*$")
_HS_PREDICATE = re.compile(r"^has_(\d+)_hs$")
_CHARGE_PREDICATE = re.compile(r"^charge(_m)?(\d+)$")


def mol_id_to_term(mol_id) -> str:
    """Turn a molecule identifier into a clingo term (integer, constant or string)."""
    if isinstance(mol_id, int) or (isinstance(mol_id, str) and mol_id.isdigit()):
        return str(int(mol_id))
    mol_id = str(mol_id)
    if _ASP_IDENTIFIER.match(mol_id):
        return mol_id
    return '"' + mol_id.replace("\\", "\\\\").replace('"', '\\"') + '"'


def atom_term(mol_term: str, atom_idx: int) -> str:
    return f"a({mol_term},{atom_idx})"


def mol_to_asp_facts(mol: Chem.Mol, mol_id) -> list[str]:
    """
    Translate a molecule into ASP facts, using the FOL structure from `mol_to_fol_atoms`.

    Atoms are represented as terms a(M, Idx). The facts are:
    - molecule(M) and has_atom(M, A) for each atom
    - each unary / binary predicate of the FOL structure on atoms (e.g. c(A), charge0(A), has_1_hs(A),
      cip_code_S(A), bSINGLE(A, B), has_bond_to(A, B))
    - each global predicate of the FOL structure on the molecule (e.g. net_charge_neutral(M))
    - num_hs(A, N) and formal_charge(A, Q), the numeric versions of has_N_hs and chargeQ / charge_mQ
    - aromatic(A) for atoms with an aromatic flag (not part of the FOL structure, but needed to translate
      SMARTS patterns, where e.g. C only matches aliphatic carbon atoms)
    """
    m = mol_id_to_term(mol_id)
    universe, extensions = mol_to_fol_atoms(mol)
    n_atoms = universe - 1  # the last element of the universe represents the molecule (global properties)
    facts = [f"molecule({m})."]
    facts += [f"has_atom({m},{atom_term(m, i)})." for i in range(n_atoms)]
    for predicate, extension in extensions.items():
        if predicate in _SKIPPED_PREDICATES or not _ASP_IDENTIFIER.match(predicate):
            continue
        if extension.ndim == 1:
            for i in extension.nonzero()[0]:
                facts.append(f"{predicate}({m if i == n_atoms else atom_term(m, i)}).")
                hs_match = _HS_PREDICATE.match(predicate)
                if hs_match and i < n_atoms:
                    facts.append(f"num_hs({atom_term(m, i)},{hs_match.group(1)}).")
                charge_match = _CHARGE_PREDICATE.match(predicate)
                if charge_match and i < n_atoms:
                    charge = int(charge_match.group(2)) * (-1 if charge_match.group(1) else 1)
                    facts.append(f"formal_charge({atom_term(m, i)},{charge}).")
        elif extension.ndim == 2:
            for i, j in zip(*extension.nonzero()):
                facts.append(f"{predicate}({atom_term(m, i)},{atom_term(m, j)}).")
    for atom in mol.GetAtoms():
        if atom.GetIsAromatic():
            facts.append(f"aromatic({atom_term(m, atom.GetIdx())}).")
    return facts
