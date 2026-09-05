import copy
import json

import numpy as np
import pytest
from Bio.PDB import MMCIFIO, Chain, MMCIFParser, Model
from Bio.PDB.Structure import Structure
from click.testing import CliRunner

import sabr.cli as cli
from sabr import constants, renumber_structure
from sabr.corrections import apply_corrections
from sabr.numbering import number_alignment
from sabr.structure import apply_numbering, extract_chain
from tests.helpers import DATA, assert_atomic_content_equal, residue_ids

pytestmark = pytest.mark.integration


@pytest.mark.parametrize(
    "case,chain,kind,gapped,boundary",
    [("8dy0_A", "A", "HK", False, 137), ("8dy1_C", "C", "LH", True, 111)],
)
@pytest.mark.parametrize("mode", constants.MODES)
def test_real_multidomain_pipeline_and_cli(
    case, chain, kind, gapped, boundary, mode, tmp_path
):
    path = DATA / f"{case}.cif"
    source = MMCIFParser(QUIET=True).get_structure(case, path)
    original = copy.deepcopy(source)
    data = extract_chain(source, chain, None)
    if gapped:
        with pytest.raises(ValueError, match="SAbR will not run"):
            renumber_structure(source, chain, mode=mode, scfv=True)
    result = renumber_structure(
        source,
        chain,
        mode=mode,
        scfv=True,
        dangerously_allow_structural_gaps=gapped,
    )
    explicit = renumber_structure(
        source,
        chain,
        mode=mode,
        chain_type="HK,HL,KH,LH",
        dangerously_allow_structural_gaps=gapped,
    )
    assert residue_ids(result, chain) == residue_ids(explicit, chain)
    assert_atomic_content_equal(original, source)
    assert_atomic_content_equal(original, result)
    assert residue_ids(original, chain) == residue_ids(source, chain)
    expected = json.loads((DATA / "multidomain_baseline.json").read_text())[
        case
    ][mode]
    assert expected["selected"] == kind
    polymer_ids = [
        (
            list(result[0][chain])[i].id[1],
            list(result[0][chain])[i].id[2].strip(),
        )
        for i in data.residue_indices
    ]
    assert polymer_ids == [(r[2], r[3]) for r in expected["mapping"]]
    assert len(set(polymer_ids)) == len(data.sequence)
    assert all(n < 1000 for n, _ in polymer_ids[:boundary])
    assert polymer_ids[boundary] == (1001, "")
    assert all(1000 < n < 2000 for n, _ in polymer_ids[boundary:])
    # Number only the rows actually left unassigned between domain alignments.
    linker_start = (
        119
        if case == "8dy0_A" and mode == "sabr"
        else 136 if case == "8dy0_A" else boundary
    )
    assert polymer_ids[linker_start:boundary] == [
        (polymer_ids[linker_start - 1][0] + j, "")
        for j in range(1, boundary - linker_start + 1)
    ]
    output = tmp_path / f"{case}-{mode}.cif"
    args = [
        "-i",
        str(path),
        "-c",
        chain,
        "-o",
        str(output),
        "--mode",
        mode,
        "--scfv",
    ]
    if gapped:
        args.append("--dangerously-allow-structural-gaps")
    cli_result = CliRunner().invoke(cli.main, args)
    assert cli_result.exit_code == 0, cli_result.output
    if gapped:
        assert cli_result.output.startswith(
            "WARNING: Structural-gap safety check disabled"
        )
    parsed = MMCIFParser(QUIET=True).get_structure("output", output)
    assert residue_ids(parsed, chain) == residue_ids(result, chain)
    assert_atomic_content_equal(original, parsed, serialized=True)


@pytest.mark.parametrize("representation", ("KH", "HHK", "HHL"))
@pytest.mark.parametrize("linker_length", (0, 1, 3))
def test_synthetic_domain_numbering_matches_independent_single_domain_goldens(
    aligned_case, representation, linker_length, tmp_path
):
    """Combinatorial numbering fixture, not a model of protein geometry."""
    cases = [aligned_case(kind) for kind in representation]
    golden = json.loads((DATA / "numbering_baseline.json").read_text())
    size = sum(
        len(golden[kind]["schemes"]["imgt"]["numbered"])
        for kind in representation
    ) + linker_length * (len(cases) - 1)
    alignment = np.zeros((size, 128 * len(cases)), dtype=int)
    structure = Structure("synthetic")
    structure.add(Model.Model(0))
    chain = Chain.Chain("X")
    structure[0].add(chain)
    expected = []
    row = 0
    for domain, (kind, case) in enumerate(
        zip(representation, cases, strict=True)
    ):
        source = cli._read_structure(case["path"])
        if domain:
            for _ in range(linker_length):
                residue = copy.deepcopy(next(source.get_residues()))
                residue.detach_parent()
                residue.resname = "GLY"
                residue.id = (" ", 10000 + row, " ")
                chain.add(residue)
                expected.append((row, expected[-1][1] + 1, "", "G"))
                row += 1
        count = len(golden[kind]["schemes"]["imgt"]["numbered"])
        domain_alignment = case["alignment"][:count]
        alignment[row : row + count, domain * 128 : (domain + 1) * 128] = domain_alignment
        baseline = golden[kind]["schemes"]["imgt"]
        assert baseline["first_row"] == 0
        assert "".join(aa for _, _, aa in baseline["numbered"]) == case["data"].sequence[:count]
        expected.extend(
            (row + index, number + 1000 * domain, code.strip(), aa)
            for index, (number, code, aa) in enumerate(baseline["numbered"])
        )
        for original_index in case["data"].residue_indices[:count]:
            residue = copy.deepcopy(
                list(source[0][case["chain"]])[original_index]
            )
            residue.detach_parent()
            residue.id = (" ", 10000 + row, " ")
            chain.add(residue)
            row += 1
    original = copy.deepcopy(structure)
    data = extract_chain(structure, "X", None)
    for domain in range(len(cases)):
        apply_corrections(alignment[:, domain * 128 : (domain + 1) * 128])
    records = number_alignment(alignment, data.sequence, "imgt", representation)
    assert [(r, n, code.strip(), aa) for r, n, code, aa in records] == expected
    result = apply_numbering(structure, "X", data, records)
    assert residue_ids(result, "X") == [(n, code) for _, n, code, _ in expected]
    assert len(set(residue_ids(result, "X"))) == size
    assert_atomic_content_equal(original, result)
    assert residue_ids(original, "X") == residue_ids(structure, "X")
    output = tmp_path / "synthetic.cif"
    writer = MMCIFIO()
    writer.set_structure(result)
    writer.save(str(output))
    parsed = MMCIFParser(QUIET=True).get_structure("output", output)
    assert residue_ids(parsed, "X") == residue_ids(result, "X")
    assert_atomic_content_equal(original, parsed, serialized=True)
