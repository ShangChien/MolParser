"""Regression tests for abbreviation element identity."""

import unittest

from rdkit import Chem

from molparser.utils import substitute_markush
from molparser.utils.chem_utils import get_abbrev_smi


class ElementAbbreviationTests(unittest.TestCase):
    def test_fragments_preserve_elements(self):
        for label, expected in [("Mn", [25]), ("Mg", [12]), ("MgBr", [12, 35])]:
            with self.subTest(label=label):
                molecule = Chem.MolFromSmiles(get_abbrev_smi()[label])
                self.assertIsNotNone(molecule)
                self.assertEqual(
                    [atom.GetAtomicNum() for atom in molecule.GetAtoms()], expected
                )

    def test_substitution_preserves_elements(self):
        for label, atomic_number in [("Mn", 25), ("Mg", 12)]:
            with self.subTest(label=label):
                caption = substitute_markush(f"N*<sep><a>1:{label}</a>", {})
                self.assertIsInstance(caption, str)
                molecule = Chem.MolFromSmiles(caption.split("<sep>", 1)[0])
                self.assertIsNotNone(molecule)
                self.assertCountEqual(
                    [atom.GetAtomicNum() for atom in molecule.GetAtoms()],
                    [7, atomic_number],
                )


if __name__ == "__main__":
    unittest.main()
