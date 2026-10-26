"""Regression checks for reviewed baselines and stale copy-number tracks."""
import argparse
import importlib.util
from pathlib import Path
import tempfile
import unittest

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def module(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / (name + '.py'))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


reviews = module('10_apply_reviewed_ploidy')
exports = module('07_export_submission_data')


class ReviewedPloidy(unittest.TestCase):
    def test_confirmed_review_overrides_stale_input_only(self):
        calls = [{'sample_id':'sample-A','baseline_ploidy':4}, {'sample_id':'sample-B','baseline_ploidy':2}]
        accepted, changes = reviews.apply_reviews(calls, [{'sample_id':'sample-A','baseline_ploidy':3,'reviewer':'investigator','evidence':'confirmed'}])
        self.assertEqual(accepted, {'sample-A':3,'sample-B':2})
        self.assertEqual(calls[0]['baseline_ploidy'], 4)
        self.assertTrue(changes[0]['changed'])
        with self.assertRaisesRegex(ValueError, 'unknown'):
            reviews.apply_reviews(calls, [{'sample_id':'missing','baseline_ploidy':3}])

    def test_review_table_validates_assignments_and_uniqueness(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'calls.tsv'
            for content in ['sample_id\tbaseline_ploidy\nA\t3.5\n', 'sample_id\tbaseline_ploidy\nA\t3\nA\t4\n']:
                path.write_text(content)
                with self.assertRaises(ValueError):
                    reviews.read_rows(path, ['sample_id','baseline_ploidy'])

    def test_submission_export_rejects_baseline_review_mismatch(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root/'ploidy.tsv').write_text('sample_id\tbaseline_ploidy\nA\t4\n')
            (root/'chr.tsv').write_text('sample_id\tChr1\nA\t4\n')
            (root/'review.tsv').write_text('sample_id\tbaseline_ploidy\treviewer\tevidence\nA\t3\tinvestigator\tconfirmed\n')
            args = argparse.Namespace(ploidy=root/'ploidy.tsv', chromosome_cn=root/'chr.tsv', manual_review=root/'review.tsv', expected_samples=1)
            with self.assertRaisesRegex(ValueError, 'accepted review'):
                exports.load_final_calls(args)
            (root/'ploidy.tsv').write_text('sample_id\tbaseline_ploidy\nA\t3\n')
            self.assertEqual(exports.load_final_calls(args)[1]['A'], 3)

    def test_submission_export_rejects_stale_or_unlabelled_tracks(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            baseline = pd.Series({'A':3})
            manual = pd.DataFrame({'baseline_ploidy':[3]},index=['A'])
            for header, value in [('\tbaseline_ploidy','\t4'), ('','')]:
                (root/'A.tsv').write_text(f'chr\tstart\tend\tCN{header}\nChr1\t1\t10000\t4{value}\n')
                with self.assertRaises(ValueError):
                    exports.export_10kb_source_data(root, root/'source.tsv.gz', baseline, pd.DataFrame(), manual)
            (root/'A.tsv').write_text('chr\tstart\tend\tCN\tbaseline_ploidy\nChr1\t1\t10000\t3\t3\n')
            self.assertEqual(exports.export_10kb_source_data(root, root/'source.tsv.gz', baseline, pd.DataFrame(), manual), (1,1))


if __name__ == '__main__':
    unittest.main()
