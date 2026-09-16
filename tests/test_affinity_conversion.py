import json
import tempfile
import unittest
from pathlib import Path

from examples.affinity_utils import predicted_nominal_ic50_nm
from examples.boltz2_predict import extract_ic50


class PredictedNominalIC50ConversionTests(unittest.TestCase):
    def test_zero_log10_micromolar_is_1000_nm(self):
        self.assertEqual(predicted_nominal_ic50_nm(0), 1000.0)

    def test_positive_and_negative_values_keep_the_documented_sign(self):
        self.assertEqual(predicted_nominal_ic50_nm(1), 10000.0)
        self.assertEqual(predicted_nominal_ic50_nm(-1), 100.0)

    def test_public_parser_uses_the_documented_conversion(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            affinity_json = Path(tmpdir) / "affinity_example.json"
            affinity_json.write_text(
                json.dumps({"affinity_pred_value": 0}), encoding="utf-8"
            )
            mean_ic50, values = extract_ic50(tmpdir)

        self.assertEqual(mean_ic50, 1000.0)
        self.assertEqual(values, [1000.0])


if __name__ == "__main__":
    unittest.main()
