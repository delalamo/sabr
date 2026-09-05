import hashlib
import json
from pathlib import Path

import pytest

from sabr import constants
from sabr.alignment import load_gap_penalties, load_references
from sabr.model import load_parameters
from sabr.numbering import _load_missing_imgt_positions

DATA = Path(__file__).parent / "data"
ASSETS = Path(__file__).parents[1] / "src" / "sabr" / "assets"


def test_scientific_assets_are_unchanged():
    expected = {
        "mpnn_encoder.npz": (
            "ee70ee43b4c9fdcf7528ad1c7ccfcba20"
            "f47d99de7568939db2f034f21b37ca7"
        ),
        "embeddings_noise_0.0.npz": (
            "1a738d26d0decb75f2fe0b83da1c71ac"
            "91433d80cc3d8270e669e99cf00e24bb"
        ),
        "embeddings_noise_0.2.npz": (
            "619286656b05d0921a4bf7a587bfc767"
            "d93f9e193fd2d2ff6787d7a4b7fa9ef3"
        ),
        "embeddings_noise_0.5.npz": (
            "3a4f51604dd0b2d27352ad269ed9ceb62"
            "42d08dc0a837bf5d3846af145a4444c"
        ),
        "embeddings_noise_1.0.npz": (
            "70495a54ae02a4db70ade0be2252c86b"
            "b50b4b02197b9ac33b042493dad88deb"
        ),
        "embeddings_noise_2.0.npz": (
            "b2f8e886e73655d25c68d23488ad17a"
            "2e441af9fe016da50d28b3e05928861e2"
        ),
        "softalign_encoder.npz": (
            "eb25b5b37b4aa62c0d70cfea78e7fd32"
            "8e7fbfce2b05a783f19b701145edba32"
        ),
        "softalign_embeddings.npz": (
            "640415ee69f144cb6d4a6ab0df2817bd0"
            "541dc8b0252f40c14e312f2e9448533"
        ),
        "softalign_gap.npz": (
            "ba3565786a7d9d3e69f6a41796f3b40"
            "e361ec65ed33cb4b435019549f2b5abae"
        ),
    }
    for filename, digest in expected.items():
        assert (
            hashlib.sha256((ASSETS / filename).read_bytes()).hexdigest()
            == digest
        )
    assert not (ASSETS / "embeddings.npz").exists()


def test_every_noise_asset_has_all_chain_references():
    for level in constants.NOISE_LEVELS:
        references = load_references(level)
        assert tuple(references) == constants.CHAIN_TYPES
        for embeddings, positions in references.values():
            assert embeddings.shape[1] == constants.EMBED_DIM
            assert embeddings.shape[0] == len(positions)


def test_missing_imgt_metadata_matches_every_reference_asset():
    reference_sets = [
        load_references(level) for level in constants.NOISE_LEVELS
    ]
    reference_sets.append(load_references(0.0, "softalign"))
    missing_imgt_positions = _load_missing_imgt_positions()

    all_positions = frozenset(range(1, constants.IMGT_MAX_POSITION + 1))
    for references in reference_sets:
        for chain_type, (_, positions) in references.items():
            assert all_positions.difference(positions) == (
                missing_imgt_positions[chain_type]
            )


def test_softalign_mode_loads_its_complete_parameter_set():
    references = load_references(2.0, "softalign")
    assert tuple(references) == constants.CHAIN_TYPES
    for embeddings, positions in references.values():
        assert embeddings.shape == (len(positions), constants.EMBED_DIM)

    parameters = load_parameters("softalign")
    assert load_parameters("softalign") is parameters
    assert set(parameters) == set(load_parameters("sabr"))

    gap_extend, gap_open = load_gap_penalties("softalign")
    assert gap_extend == pytest.approx(0.1942468136548996)
    assert gap_open == pytest.approx(-2.5441808700561523)


def test_model_parameters_are_cached_and_reference_arrays_are_read_only():
    load_parameters.cache_clear()
    first_parameters = load_parameters()
    assert load_parameters() is first_parameters

    load_references.cache_clear()
    first_references = load_references(0.0)
    assert load_references(0.0) is first_references
    for embeddings, positions in first_references.values():
        assert not embeddings.flags.writeable
        assert isinstance(positions, tuple)


def test_new_regression_files_match_recorded_checksums():
    manifest = json.loads((DATA / "regression_provenance.json").read_text())
    for name, checksum in {**manifest["inputs"], **manifest["outputs"]}.items():
        assert (
            hashlib.sha256((DATA / name).read_bytes()).hexdigest() == checksum
        ), name


@pytest.mark.parametrize("loader", (load_parameters, load_gap_penalties))
def test_parameter_loaders_reject_unknown_modes(loader):
    with pytest.raises(ValueError, match="mode"):
        loader("unknown")


def test_reference_loader_rejects_unknown_mode():
    with pytest.raises(ValueError, match="mode"):
        load_references(0.0, "unknown")
