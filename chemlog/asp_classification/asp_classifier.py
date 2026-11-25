import logging
import os
from typing import Optional

import clingo
from rdkit import Chem

from chemlog.base_classifier import Classifier
from chemlog.preprocessing.mol_to_asp import mol_to_asp_facts, mol_id_to_term

RULES_PATH = os.path.join(os.path.dirname(__file__), "peptide_rules.pl")


class ASPPeptideClassifier:
    """
    Classifies molecules with the ASP rules in peptide_rules.pl, using clingo.

    The result for each molecule has the same format as the algorithmic classification in chemlog.cli
    (keys CHARGE, SIZE, PROTEINOGENICS, SUBSTRUCT and chebi_classes).

    timeout: maximum number of seconds clingo may search for the next (better) model, None for no limit
    batch_size: number of molecules classified in one clingo call by classify_molecules. Each call has a fixed
        cost for parsing and grounding the rules (~8 ms), batching spreads it over several molecules.
    """

    def __init__(self, rules_path: str = RULES_PATH, timeout: Optional[float] = None, batch_size: int = 50):
        with open(rules_path, "r", encoding="utf-8") as f:
            self.rules = f.read()
        self.timeout = timeout
        self.batch_size = batch_size

    def solve(self, facts: list[str]) -> Optional[list[clingo.Symbol]]:
        """Ground the rules together with `facts`, return the shown atoms of the optimal answer set."""
        messages = []
        ctl = clingo.Control(["--opt-mode=opt"], logger=lambda code, msg: messages.append(msg))
        ctl.add("base", [], self.rules)
        ctl.add("base", [], "\n".join(facts))
        ctl.ground([("base", [])])
        symbols = None
        with ctl.solve(yield_=True, async_=True) as handle:
            while True:
                handle.resume()
                if not handle.wait(self.timeout):
                    handle.cancel()
                    logging.warning(f"clingo timed out after {self.timeout}s")
                    return None
                model = handle.model()
                if model is None:
                    break
                # with --opt-mode=opt, each model improves on the previous one, the last one is optimal
                symbols = model.symbols(shown=True)
        for message in messages:
            logging.debug(f"clingo: {message}")
        return symbols

    def classify_molecules(self, mols: dict) -> dict:
        """Classify several molecules (dict mol_id -> rdkit Mol), using one clingo call per batch."""
        mol_ids = list(mols)
        results = {}
        for start in range(0, len(mol_ids), self.batch_size):
            results.update(self._classify_batch({i: mols[i] for i in mol_ids[start:start + self.batch_size]}))
        return results

    def _classify_batch(self, mols: dict) -> dict:
        facts = []
        term_to_id = {}
        for mol_id, mol in mols.items():
            facts += mol_to_asp_facts(mol, mol_id)
            term_to_id[mol_id_to_term(mol_id)] = mol_id
        symbols = self.solve(facts)
        if symbols is None:
            return {mol_id: None for mol_id in mols}
        results = {mol_id: {"chebi_classes": [], "proteinogenic_residues": [],
                            "emericellamide_match": False, "diketopiperazine_match": False}
                   for mol_id in mols}
        for symbol in symbols:
            res = results[term_to_id[str(symbol.arguments[0])]]
            if symbol.name == "peptide_size":
                res["SIZE"] = symbol.arguments[1].number
            elif symbol.name == "charge_category":
                res["CHARGE"] = symbol.arguments[1].string
            elif symbol.name == "chebi_class":
                res["chebi_classes"].append(symbol.arguments[1].string)
            elif symbol.name == "proteinogenic_residue":
                # store the letter and the atom indices of N, alpha-C and carboxy-C
                res["proteinogenic_residues"].append(
                    (symbol.arguments[1].string, *[a.arguments[1].number for a in symbol.arguments[2:]]))
            elif symbol.name in ["emericellamide_match", "diketopiperazine_match"]:
                res[symbol.name] = True
        return {mol_id: self._to_chemlog_format(res) for mol_id, res in results.items()}

    def classify(self, mol: Chem.Mol) -> dict:
        return self.classify_molecules({"m": mol})["m"]

    @staticmethod
    def _to_chemlog_format(res: dict) -> dict:
        residues = sorted(res.pop("proteinogenic_residues"), key=lambda r: r[1:])
        res["PROTEINOGENICS"] = [r[0] for r in residues]
        res["PROTEINOGENICS_additional"] = {"residues": residues}
        # like the algorithmic substructure classifier, only report the substructure that is relevant
        # for the peptide size
        res["SUBSTRUCT"] = {}
        if res["SIZE"] == 5:
            res["SUBSTRUCT"]["emericellamide"] = res["emericellamide_match"]
        if res["SIZE"] == 2:
            res["SUBSTRUCT"]["2,5-diketopiperazines"] = res["diketopiperazine_match"]
        del res["emericellamide_match"], res["diketopiperazine_match"]
        res["chebi_classes"] = sorted(res["chebi_classes"])
        return res


# The CLI calls one classifier per aspect (charge, size, ...) for each molecule. Since a single clingo call
# answers all of them, the wrappers below share the result of the last molecule.
_shared_classifier = None
_last_result = (None, None)


def _classify_shared(mol: Chem.Mol) -> dict:
    global _shared_classifier, _last_result
    if _last_result[0] is not mol:
        if _shared_classifier is None:
            _shared_classifier = ASPPeptideClassifier()
        _last_result = (mol, _shared_classifier.classify(mol))
    return _last_result[1]


class ASPChargeClassifier(Classifier):
    def classify(self, mol: Chem.Mol, *args, **kwargs):
        return _classify_shared(mol)["CHARGE"], None


class ASPPeptideSizeClassifier(Classifier):
    def classify(self, mol: Chem.Mol, *args, **kwargs):
        return _classify_shared(mol)["SIZE"], None


class ASPProteinogenicsClassifier(Classifier):
    def classify(self, mol: Chem.Mol, *args, **kwargs):
        res = _classify_shared(mol)
        return res["PROTEINOGENICS"], res["PROTEINOGENICS_additional"]


class ASPSubstructureClassifier(Classifier):
    def classify(self, mol: Chem.Mol, *args, **kwargs):
        return _classify_shared(mol)["SUBSTRUCT"], None


if __name__ == "__main__":
    classifier = ASPPeptideClassifier()
    mol = Chem.MolFromSmiles("C[C@H](N)C(=O)N[C@@H](C)C(=O)O")  # Ala-Ala
    Chem.Kekulize(mol)
    print(classifier.classify(mol))
