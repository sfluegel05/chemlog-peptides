import ast
import enum
import json
import logging
import os
import time

import click
import multiprocess as mp
import networkx as nx
import tqdm
from rdkit import Chem

from chemlog.alg_classification.charge_classifier import get_charge_category, AlgChargeClassifier
from chemlog.alg_classification.peptide_size_classifier import get_carboxy_derivatives, get_amide_bonds, \
    get_amino_groups
from chemlog.alg_classification.peptide_size_classifier import get_n_amino_acid_residues, AlgPeptideSizeClassifier
from chemlog.alg_classification.peptide_size_translated import AlgPeptideSizeClassifierTranslated
from chemlog.alg_classification.proteinogenics_classifier import get_proteinogenic_amino_acids, \
    AlgProteinogenicsClassifier
from chemlog.alg_classification.substructure_classifier import is_emericellamide, is_diketopiperazine, \
    AlgSubstructureClassifier
from chemlog.asp_classification.asp_classifier import ASPChargeClassifier, ASPPeptideSizeClassifier,     ASPProteinogenicsClassifier, ASPSubstructureClassifier
from chemlog.base_classifier import ChargeCategories
from chemlog.fol_classification.charge_verifier import ChargeVerifier
from chemlog.fol_classification.functional_groups_verifier import FunctionalGroupsVerifier
from chemlog.fol_classification.model_checking import ModelCheckerOutcome
from chemlog.fol_classification.peptide_size_verifier import PeptideSizeVerifier, FOLPeptideSizeClassifierTranslated
from chemlog.fol_classification.proteinogenics_verifier import ProteinogenicsVerifier
from chemlog.fol_classification.substruct_verifier import SubstructVerifier
from chemlog.mona_classification.peptide_size_mona import MonaPeptideSizeClassifier, MonaPeptideSizeClassifierCompiled
from chemlog.preprocessing.chebi_data import ChEBIData
from chemlog.preprocessing.mol_to_fol import mol_to_fol_atoms
from chemlog.preprocessing.pubchem_data import PubChemData
from chemlog.qbf_classification.peptide_size_qbf import QBFPeptideSizeClassifierDepQBF, QBFPeptideSizeClassifierCAQE, \
    QBFPeptideSizeClassifierDepQBFTranslated
from chemlog.timestamped_logger import TimestampedLogger


class LiteralOption(click.Option):
    def type_cast_value(self, ctx, value):
        try:
            return ast.literal_eval(value)
        except:
            raise click.BadParameter(value)


@click.group(help="CLI for classifying peptides")
def cli():
    pass


def resolve_chebi_classes(classification):
    n_amino_acid_residues = classification[ClassifierKeys.SIZE.name]
    charge_category = classification[ClassifierKeys.CHARGE.name]
    res = []
    if charge_category == ChargeCategories.SALT.name:
        res.append("24866")  # salt (there is no class peptide salt)
    elif charge_category == ChargeCategories.ANION.name:
        res.append("22563")
    elif charge_category == ChargeCategories.CATION.name:
        res.append("36961")
    elif charge_category == ChargeCategories.ZWITTERION.name:
        res.append("27369")
    if n_amino_acid_residues >= 2:
        if charge_category == ChargeCategories.ANION.name:
            # peptide anion
            res.append("60334")
        elif charge_category == ChargeCategories.CATION.name:
            # peptide cation
            res.append("60194")
        elif charge_category == ChargeCategories.ZWITTERION.name:
            res.append("60466")
            if n_amino_acid_residues == 2:
                # zwitterion, peptide zwitterion, dipeptide zwitterion
                res.append("90799")
            if n_amino_acid_residues == 3:
                res.append("155837")
        elif charge_category == ChargeCategories.NEUTRAL.name:
            res.append("16670")
            if n_amino_acid_residues == 2:
                res.append("46761")
            if n_amino_acid_residues == 3:
                res.append("47923")
            if n_amino_acid_residues == 4:
                res.append("48030")
            if n_amino_acid_residues == 5:
                res.append("48545")
            if n_amino_acid_residues >= 10:
                res.append("15841")
            else:
                # oligo
                res.append("25676")
    if ClassifierKeys.SUBSTRUCT.name in classification:
        substruct_classification = classification[ClassifierKeys.SUBSTRUCT.name]
        if "emericellamide" in substruct_classification and substruct_classification["emericellamide"]:
            res.append("64372")
        if "2,5-diketopiperazines" in substruct_classification and substruct_classification["2,5-diketopiperazines"]:
            res.append("65061")

    return res


@cli.command(help="Classify Pubchem molecules using a direct Python implementation")
@click.option('--from-batch', '-f', type=int, default=0, help='Start at this PubChem batch')
@click.option('--to-batch', '-t', type=int, default=346, help='End at this PubChem batch (exclusive)')
@click.option('--return-chebi-classes', '-c', is_flag=True, help='Return assigned ChEBI classes')
@click.option('--molecules', '-m', cls=LiteralOption, default="[]",
              help='List of PubChem IDs to classify. Default: all PubChem entries.')
def classify_pubchem(from_batch, to_batch, return_chebi_classes, molecules):
    json_logger = TimestampedLogger()
    json_logger.start_run(f"classify_pubchem", {
        "return_chebi_classes": return_chebi_classes, "from_batch": from_batch, "to_batch": to_batch})

    for batch_id in range(from_batch, to_batch):
        data_filtered = PubChemData().get_processed_batch(batch_id)

        if len(molecules) > 0:
            data_filtered = {k: v for k, v in data_filtered.items() if k in molecules}

        results = []
        logging.info(f"Starting batch {batch_id} ({len(data_filtered)} molecules)")
        with tqdm.tqdm(total=len(data_filtered), desc=f"Classifying batch {batch_id}") as pbar:
            for pubchem_id, mol in data_filtered.items():
                pbar.set_description(f"Classifying molecule {pubchem_id}")
                pbar.update()
                start_time = time.perf_counter()
                charge_category = get_charge_category(mol)
                n_amino_acid_residues, add_output = get_n_amino_acid_residues(mol)
                if n_amino_acid_residues > 1:
                    proteinogenics, _, _ = get_proteinogenic_amino_acids(mol,
                                                                         add_output["amino_residue"],
                                                                         add_output["carboxy_residue"])
                else:
                    proteinogenics = []
                results.append({
                    'pubchem_id': pubchem_id,
                    'charge_category': charge_category.name,
                    'n_amino_acid_residues': n_amino_acid_residues,
                    'proteinogenics': proteinogenics,
                    'time': f"{time.perf_counter() - start_time:.4f}"
                })

                if n_amino_acid_residues == 5:
                    emericellamide = is_emericellamide(mol)
                    results[-1]["emericellamide"] = emericellamide[0]
                if n_amino_acid_residues == 2:
                    diketopiperazine = is_diketopiperazine(mol)
                    results[-1]["2,5-diketopiperazines"] = diketopiperazine[0]

                if return_chebi_classes:
                    results[-1]['chebi_classes'] = resolve_chebi_classes(results[-1])

        json_logger.save_items(f"classify_pubchem{batch_id:03d}", results)


def strategy_call_chebi(strategy, classifier_instances, ident, row):
    logging.debug(f"Classifying CHEBI:{ident} ({row['name']})  {row['smiles']}")

    res = strategy_call(strategy, classifier_instances, row["mol"])
    res["chebi_id"] = ident
    return res

def strategy_call(strategy, classifier_instances, mol):
    res = dict()
    start_time = time.perf_counter()
    if strategy == 'fol':
        fol_structure = mol_to_fol_atoms(mol)
    for key in ClassifierKeys:
        if key in classifier_instances:
            args = []
            # FOL size and proteinogenics classifiers need functional groups, substruct classifier needs size
            if strategy == 'fol':
                if key in [ClassifierKeys.SIZE, ClassifierKeys.PROTEINOGENICS]:
                    args.append(res[ClassifierKeys.FGS.name])
                elif key in [ClassifierKeys.SUBSTRUCT]:
                    args.append(res[ClassifierKeys.SIZE.name])
                # To avoid computing the fol structure 3 times, calculate it here and pass it to the classifiers
                if key in [ClassifierKeys.FGS, ClassifierKeys.SUBSTRUCT, ClassifierKeys.PROTEINOGENICS]:
                    args.append(fol_structure)
            # Algorithmic approach produces functional groups as a side product -> reuse
            elif strategy == 'algo':
                if key in [ClassifierKeys.PROTEINOGENICS]:
                    args += [res[f"{ClassifierKeys.SIZE.name}_additional"]['amino_residue'],
                             res[f"{ClassifierKeys.SIZE.name}_additional"]['carboxy_residue']]
                elif key in [ClassifierKeys.SUBSTRUCT]:
                    args.append(res[ClassifierKeys.SIZE.name])

            classification, additional_output = classifier_instances[key].classify(mol, *args)
            logging.debug(f"Classification for {key.name}: {classification}")
            res[key.name] = classification
            if additional_output is not None:
                res[key.name + "_additional"] = additional_output

    res['time'] = f"{time.perf_counter() - start_time:.4f}"
    if ClassifierKeys.SIZE.name in res and ClassifierKeys.CHARGE.name in res:
        res['chebi_classes'] = resolve_chebi_classes(res)
    return res


class ClassifierKeys(enum.Enum):
    CHARGE = 0
    FGS = 1
    SIZE = 2
    PROTEINOGENICS = 3
    SUBSTRUCT = 4


CLASSIFIERS = {
    'mona': {
        ClassifierKeys.SIZE: MonaPeptideSizeClassifierCompiled,
    },
    'mona-from-file': {
        ClassifierKeys.SIZE: MonaPeptideSizeClassifier,
    },
    'qbf-caqe': {
        ClassifierKeys.SIZE: QBFPeptideSizeClassifierCAQE,
    },
    'qbf-depqbf': {
        ClassifierKeys.SIZE: QBFPeptideSizeClassifierDepQBF,
    },
    'qbf-depqbf-translated': {
        ClassifierKeys.SIZE: QBFPeptideSizeClassifierDepQBFTranslated,
    },
    'fol': {
        ClassifierKeys.CHARGE: ChargeVerifier,
        ClassifierKeys.SIZE: PeptideSizeVerifier,
        ClassifierKeys.FGS: FunctionalGroupsVerifier,
        ClassifierKeys.PROTEINOGENICS: ProteinogenicsVerifier,
        ClassifierKeys.SUBSTRUCT: SubstructVerifier,
    },
    'fol-translated': {
        ClassifierKeys.SIZE: FOLPeptideSizeClassifierTranslated
    },
    'algo': {
        ClassifierKeys.CHARGE: AlgChargeClassifier,
        ClassifierKeys.SIZE: AlgPeptideSizeClassifier,
        ClassifierKeys.PROTEINOGENICS: AlgProteinogenicsClassifier,
        ClassifierKeys.SUBSTRUCT: AlgSubstructureClassifier,
    },
    'algo-translated': {
        ClassifierKeys.SIZE: AlgPeptideSizeClassifierTranslated
    },
    'asp': {
        ClassifierKeys.CHARGE: ASPChargeClassifier,
        ClassifierKeys.SIZE: ASPPeptideSizeClassifier,
        ClassifierKeys.PROTEINOGENICS: ASPProteinogenicsClassifier,
        ClassifierKeys.SUBSTRUCT: ASPSubstructureClassifier,
    }
}


@cli.command(
    help="Classify ChEBI molecules")
@click.option('--chebi-version', '-v', type=int, required=True, help='ChEBI version')
@click.option('--strategy', '-s', type=click.Choice(list(CLASSIFIERS.keys()), case_sensitive=False), default='algo',
              help='Strategy to use for classification.')
@click.option('--run-name', '-n', type=str, help='Results will be stored at results/%y%m%d_%H%M_{strategy}_{run_name}/')
@click.option('--debug-mode', '-d', is_flag=True, help='Logs at debug level')
@click.option('--molecules', '-m', cls=LiteralOption, default="[]",
              help='List of ChEBI IDs to classify. Default: all ChEBI classes, sorted by SMILES length')
@click.option('--only-3star', '-3', is_flag=True, help='Only consider 3-star molecules')
@click.option('--only-peptides', '-p', is_flag=True, help='Only consider peptide molecules')
@click.option('--begin-molecule', '-b', type=int, default=0,
              help='Start at this molecule index (applied after other selectors)')
@click.option('--n-molecules', '-l', type=int, default=-1, help='End after this many molecules')
@click.option('--n-workers', '-w', type=int, default=0,
              help='Number of worker processes to use, use 0 (default) for no multiprocessing')
def classify_chebi(chebi_version, strategy, run_name, debug_mode, molecules, only_peptides, only_3star, begin_molecule,
                   n_molecules, n_workers):
    json_logger = TimestampedLogger(None, f"{strategy}_{run_name}" if run_name is not None else strategy, debug_mode)
    json_logger.start_run(f"classify_{strategy}", {"chebi_version": chebi_version, "molecules": molecules,
                                                   "run_name": run_name, "debug_mode": debug_mode,
                                                   "only_peptides": only_peptides,
                                                   "begin_molecule": begin_molecule, "n_molecules": n_molecules,
                                                   "n_workers": n_workers})

    data_filtered = _supply_chebi_data(chebi_version, molecules, only_3star, only_peptides)
    data_filtered = data_filtered[begin_molecule:]
    if n_molecules > 0:
        data_filtered = data_filtered[:n_molecules]
    logging.info(f"Classifying {len(data_filtered)} molecules")

    classifier_instances = {
        k: v() for k, v in CLASSIFIERS[strategy].items()
    }
    # no multiprocessing
    if n_workers == 0:
        logging.info("Running in single-threaded mode")
        results = []
        for i, (id, row) in tqdm.tqdm(enumerate(data_filtered.iterrows())):
            results.append(strategy_call_chebi(strategy, classifier_instances, id, row))
            if len(data_filtered) < 100 or ((i+1) % (len(data_filtered) // 100)) == 0:
                json_logger.save_items(f"classify_{strategy}", results)

        json_logger.save_items(f"classify_{strategy}", results)
        for classifier in classifier_instances.values():
            classifier.on_finish()
        return

    output_q = mp.Queue()
    input_q = mp.Queue()
    for id, row in data_filtered.iterrows():
        input_q.put((id, row))

    results = []
    i = 0

    def worker():
        while True:
            next_task = input_q.get()
            if next_task is None:
                logging.info(f"Poisoned worker {mp.current_process().name}, exiting")
                break
            id, row = next_task
            output_q.put(strategy_call_chebi(strategy, classifier_instances, id, row))

    processes = []
    n_workers = min(n_workers, input_q.qsize())
    logging.info(f"Starting {n_workers} worker processes")
    input_size = len(data_filtered)
    for _ in range(n_workers):
        input_q.put(None)

    pbar = tqdm.tqdm(total=input_size, desc=f"Classifying with {strategy.upper()}")

    for _ in range(n_workers):
        p = mp.Process(target=worker, args=())
        processes.append(p)
        p.start()

    while i < input_size:
        result = output_q.get()
        pbar.update(1)
        results.append(result)
        i += 1
        if input_size < 10 or (i % (input_size // 10)) == 0:
            json_logger.save_items(f"classify_{strategy}", results)

    pbar.close()
    for p in processes:
        p.join()

    logging.info(f"Finished classifying {i}/{input_size} molecules")
    json_logger.save_items(f"classify_{strategy}", results)
    for classifier in classifier_instances.values():
        classifier.on_finish()


def _supply_chebi_data(chebi_version, molecules, only_3star, only_peptides=False):
    data_cls = ChEBIData(chebi_version)
    data = data_cls.processed
    molecules = [str(m) for m in molecules]
    if len(molecules) > 0:
        data_filtered = data.loc[data.index.isin(molecules)]
    else:
        data_filtered = data
    if only_3star:
        data_filtered = data_filtered[data_filtered["subset"][0] == "3"]
    if only_peptides:
        trans_hierarchy = data_cls.get_trans_hierarchy()
        data_filtered = data_filtered.loc[list(set.intersection(set(nx.descendants(trans_hierarchy, "16670")),
                                                                set(data_filtered.index)))]

    # start with shortest SMILES
    data_filtered["smiles_length"] = [
        len(str(row["smiles"]) if row["smiles"] is not None else "")
        for _, row in data_filtered.iterrows()
    ]
    data_filtered.sort_values("smiles_length", inplace=True, ascending=True)
    return data_filtered

@cli.command(
    help="Classify SMILES")
@click.option('--strategy', '-s', type=click.Choice(list(CLASSIFIERS.keys()), case_sensitive=False), default='algo',
              help='Strategy to use for classification.')
@click.option('--smiles', '-s', multiple=True, help='SMILES strings to predict')
@click.option('--smiles-file', '-f', type=click.Path(exists=True), help='File containing SMILES strings (one per line)')
@click.option('--run-name', '-n', type=str, help='Results will be stored at results/%y%m%d_%H%M_{strategy}_{run_name}/')
@click.option('--debug-mode', '-d', is_flag=True, help='Logs at debug level')
def classify_smiles(strategy, smiles, smiles_file, run_name, debug_mode):
    json_logger = TimestampedLogger(None, f"{strategy}_{run_name}" if run_name is not None else strategy, debug_mode)
    # Collect SMILES strings from arguments and/or file
    smiles_list = list(smiles)
    if smiles_file:
        with open(smiles_file, 'r') as f:
            smiles_list.extend([line.strip() for line in f if line.strip()])

    if not smiles_list:
        click.echo("No SMILES strings provided. Use --smiles or --smiles-file options.")
        return
    json_logger.start_run(f"classify_{strategy}", {"strategy": strategy, "smiles": smiles_list,
                                                   "run_name": run_name, "debug_mode": debug_mode})

    logging.info(f"Classifying {len(smiles_list)} molecules")

    classifier_instances = {
        k: v() for k, v in CLASSIFIERS[strategy].items()
    }
    # no multiprocessing
    results = []

    for i, smiles in tqdm.tqdm(enumerate(smiles_list)):
        mol = _smiles_to_mol(smiles)
        if mol is None:
            results.append(None)
        else:
            results.append(strategy_call(strategy, classifier_instances, mol))
        if len(smiles_list) < 100 or ((i+1) % (len(smiles_list) // 100)) == 0:
            json_logger.save_items(f"classify_{strategy}", results)

    json_logger.save_items(f"classify_{strategy}", results)
    for classifier in classifier_instances.values():
        classifier.on_finish()

def _smiles_to_mol(smiles):
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is not None:
        # turn aromatic bond types into single/double
        try:
            Chem.Kekulize(mol)
        except Chem.KekulizeException as e:
            logging.debug(f"{Chem.MolToSmiles(mol)} - {e}")
        return mol
    return mol

@cli.command(help="Verify results from a `classify` run using first-order logic (FOL)")
@click.option('--chebi-version', '-v', type=int, required=True, help='ChEBI version')
@click.option('--results-dir', '-r', type=str, required=True, help='Directory where results.json to analyse is located')
@click.option('--debug-mode', '-d', is_flag=True, help='Returns additional states')
@click.option('--molecules', '-m', cls=LiteralOption, default="[]",
              help='List of ChEBI IDs to verify. Default: all ChEBI classes.')
@click.option('--only-3star', '-3', is_flag=True, help='Only consider 3-star molecules')
def verify(chebi_version, results_dir, debug_mode, molecules, only_3star):
    json_logger = TimestampedLogger(results_dir, debug_mode=debug_mode)
    json_logger.start_run("verify" + ("_3star" if only_3star else ""),
                          {"chebi_version": chebi_version, "results_dir": results_dir, "debug_mode": debug_mode,
                           "molecules": molecules})
    data = ChEBIData(chebi_version)
    with open(os.path.join(results_dir, "classify.json"), "r") as f:
        results = json.load(f)
    charge_verifier = ChargeVerifier()
    functional_groups_verifier = FunctionalGroupsVerifier()
    peptide_size_verifier = PeptideSizeVerifier()
    proteinogenics_verifier = ProteinogenicsVerifier()
    substruct_verifier = SubstructVerifier()
    res = []

    molecules = [str(m) for m in molecules]
    results = [r for r in results if (len(molecules) == 0 or r["chebi_id"] in molecules) and (
            not only_3star or data.processed.loc[r["chebi_id"], "subset"][0] == "3")]

    save_results_at = len(results) / 4

    for i, result in tqdm.tqdm(enumerate(results), total=len(results), desc="Verifying"):
        outcome, expected = {}, {}
        start_time = time.perf_counter()
        expected["charge"] = ChargeCategories[result["charge_category"]].name
        mol = data.processed.loc[result["chebi_id"], "mol"]
        outcome["charge"] = charge_verifier.verify_charge_category(mol, ChargeCategories[result["charge_category"]], {})
        # functional groups
        expected["functional_groups"] = {}
        if "amide_bond" in result:
            expected["functional_groups"]["amide_bond"] = result["amide_bond"]
        else:
            _, amide_c, amide_o, amide_n = get_amide_bonds(mol)
            expected["functional_groups"]["amide_bond"] = [(c, o, n) for c, o, n in zip(amide_c, amide_o, amide_n)]
        if "amino_residue" in result:
            expected["functional_groups"]["amino_residue"] = [(n,) for n in result["amino_residue"]]
        else:
            expected["functional_groups"]["amino_residue"] = [(n,) for n in
                                                              get_amino_groups(mol, [c for c, _, _ in
                                                                                     expected["functional_groups"][
                                                                                         "amide_bond"]])]
        if "carboxy_residue" in result:
            expected["functional_groups"]["carboxy_residue"] = result["carboxy_residue"]
        else:
            expected["functional_groups"]["carboxy_residue"] = list(get_carboxy_derivatives(mol))
        outcome["functional_groups"] = functional_groups_verifier.verify_functional_groups(mol, expected[
            "functional_groups"])

        # n amino acids
        expected["n_amino_acid_residues"] = result["n_amino_acid_residues"]
        if expected["n_amino_acid_residues"] > 1:
            if "longest_aa_chain" in result:
                aars = result["longest_aa_chain"]
            else:
                add_output = get_n_amino_acid_residues(mol)[1]
                aars = add_output["longest_aa_chain"]
            outcome["size"] = peptide_size_verifier.verify_n_plus_amino_acids(
                mol, expected["n_amino_acid_residues"], expected["functional_groups"],
                {f"A{i}": aar for i, aar in enumerate(aars)}
            )
        else:
            # if there are no amino acids, then there is nothing to prove
            outcome["size"] = ModelCheckerOutcome.MODEL_FOUND_INFERRED, None

        # proteinogenics
        if "proteinogenics_locations_no_carboxy" in result:
            expected["proteinogenics"] = [(code, atoms) for code, atoms in zip(result["proteinogenics"],
                                                                               result[
                                                                                   "proteinogenics_locations_no_carboxy"])]
        else:
            proteinogenics, _, proteinogenics_locations_no_carboxy = get_proteinogenic_amino_acids(
                mol, [amino[0] for amino in expected["functional_groups"]["amino_residue"]],
                expected["functional_groups"]["carboxy_residue"]
            )
            expected["proteinogenics"] = [(code, atoms) for code, atoms in
                                          zip(proteinogenics, proteinogenics_locations_no_carboxy)]
        if len(expected["proteinogenics"]) > 0:
            # only take first atoms of functional groups
            atom_level_functional_groups = {
                "amino_residue_n": [amino[0] for amino in expected["functional_groups"]["amino_residue"]],
                "carboxy_residue_c": [carboxy[0] for carboxy in expected["functional_groups"]["carboxy_residue"]]}
            outcome["proteinogenics"] = proteinogenics_verifier.verify_proteinogenics(mol, atom_level_functional_groups,
                                                                                      expected["proteinogenics"])
        else:
            outcome["proteinogenics"] = ModelCheckerOutcome.MODEL_FOUND_INFERRED, None

        # substructures
        if "emericellamide" in result and result["emericellamide"]:
            if "emericellamide_atoms" in result:
                atoms = result["emericellamide_atoms"]
            else:
                atoms = is_emericellamide(mol)[1]
            expected["emericellamide"] = True
            outcome["emericellamide"] = substruct_verifier.verify_substruct_class(mol, "emericellamide", atoms)
        if "2,5-diketopiperazines" in result and result["2,5-diketopiperazines"]:
            if "2,5-diketopiperazines_atoms" in result:
                atoms = result["2,5-diketopiperazines_atoms"]
            else:
                atoms = is_diketopiperazine(mol)
            expected["2,5-diketopiperazines"] = True
            outcome["2,5-diketopiperazines"] = substruct_verifier.verify_substruct_class(mol, "diketopiperazines",
                                                                                         atoms)

        res.append({
            "chebi_id": result["chebi_id"],
            "expected": expected,
            "outcome": {key: o[0].name for key, o in outcome.items()},
            "time": f"{time.perf_counter() - start_time:.4f}"
        })

        if debug_mode:
            res[-1]["proof_attempts"] = {key: o[1] for key, o in outcome.items() if o[1] is not None}
        if any(o[0] not in [ModelCheckerOutcome.MODEL_FOUND, ModelCheckerOutcome.MODEL_FOUND_INFERRED] for o in
               outcome.values()):
            warning_str = ""
            if outcome["charge"][0] not in [ModelCheckerOutcome.MODEL_FOUND, ModelCheckerOutcome.MODEL_FOUND_INFERRED]:
                warning_str += f"Expected charge: {expected['charge']}, got: {outcome['charge'][0].name}, tried: {outcome['charge'][1]}\n"
            if outcome["functional_groups"][0] not in [ModelCheckerOutcome.MODEL_FOUND,
                                                       ModelCheckerOutcome.MODEL_FOUND_INFERRED]:
                warning_str += f"Expected groups: {expected['functional_groups']}, got {outcome['functional_groups'][0].name}, tried: {outcome['functional_groups'][1]}\n"
            if outcome["size"][0] not in [ModelCheckerOutcome.MODEL_FOUND, ModelCheckerOutcome.MODEL_FOUND_INFERRED]:
                warning_str += f"Expected {expected['n_amino_acid_residues']} amino acids, got {outcome['size'][0].name}\n"
            if outcome["proteinogenics"][0] not in [ModelCheckerOutcome.MODEL_FOUND,
                                                    ModelCheckerOutcome.MODEL_FOUND_INFERRED]:
                warning_str += f"Expected proteinogenics: {expected['proteinogenics']}, got {outcome['proteinogenics'][0].name}, tried {outcome['proteinogenics'][1]}\n"
            if "emericellamide" in expected and outcome["emericellamide"][0] not in [ModelCheckerOutcome.MODEL_FOUND,
                                                                                     ModelCheckerOutcome.MODEL_FOUND_INFERRED]:
                warning_str += f"Expected emericellamide, got {outcome['emericellamide'][0].name}\n"
            if "2,5-diketopiperazines" in expected and outcome["2,5-diketopiperazines"][0] not in [
                ModelCheckerOutcome.MODEL_FOUND, ModelCheckerOutcome.MODEL_FOUND_INFERRED]:
                warning_str += f"Expected 2,5-diketopiperazines, got {outcome['2,5-diketopiperazines'][0].name}\n"
            logging.warning(f"Verification failed for CHEBI:{result['chebi_id']} \n{warning_str}")

        if i >= save_results_at:
            save_results_at = i + (len(results) - i) // 4 + 10
            json_logger.save_items(f"verify_{json_logger.timestamp}", res)
    json_logger.save_items(f"verify_{json_logger.timestamp}", res)


if __name__ == "__main__":
    classify_chebi(["-v", 239, "-s", "algo", "-3", "-w", "0"], standalone_mode=False)
