import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import lnctrail


class LncTrailTests(unittest.TestCase):
    def test_fasta_reads_and_converts_t_to_u(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "test.fa"
            p.write_text(">r1 description\nACGTACGT\n", encoding="utf-8")
            records, warnings = lnctrail.read_fasta(p, "test")
            self.assertEqual(len(records), 1)
            self.assertEqual(records[0].id, "r1")
            self.assertEqual(records[0].sequence, "ACGUACGU")
            self.assertTrue(records[0].original_had_t)

    def test_duplicate_ids_fail(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "dup.fa"
            p.write_text(">same one\nACGUACGU\n>same two\nACGUACGU\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "duplicated FASTA ID"):
                lnctrail.read_fasta(p, "test")

    def test_builtin_detects_exact_reverse_complement_seed(self):
        q = "AUGCCAUUGACCUAGGCUAC"
        target = "UUU" + q.translate(lnctrail.COMP)[::-1] + "AAAA"
        task = ("q1", "q1", q, lnctrail.gc_percent(q), "t1", "t1", target,
                lnctrail.gc_percent(target), 7, 8, 8, 12, 5, 1000)
        results = lnctrail.builtin_pair_task(task)
        self.assertTrue(results)
        self.assertEqual(results[0]["engine"], "LncTrail-builtin")
        self.assertEqual(results[0]["energy_kcal_mol"], "")
        self.assertGreaterEqual(results[0]["canonical_pairs"], 7)
        self.assertIn("|", results[0]["interaction_pattern"])

    def test_ranking_prefers_more_negative_energy(self):
        rows = [
            {"lncrna_id": "q", "target_id": "a", "energy_kcal_mol": -7.0, "binding_score": "", "target_start": 1},
            {"lncrna_id": "q", "target_id": "b", "energy_kcal_mol": -12.0, "binding_score": "", "target_start": 1},
        ]
        out = lnctrail.rank_rows(rows, "intarna")
        ranks = {r["target_id"]: r["target_rank"] for r in out}
        self.assertEqual(ranks["b"], 1)
        self.assertEqual(ranks["a"], 2)

    def test_infer_gene_from_refseq_description(self):
        rec = lnctrail.SequenceRecord(
            "NM_000001.2",
            "NM_000001.2 Homo sapiens transmembrane protein 173 (TMEM173), mRNA",
            "ACGUACGU"
        )
        self.assertEqual(lnctrail.infer_target_gene(rec), "TMEM173")

    def test_low_complexity_fraction_flags_homopolymer(self):
        self.assertGreater(lnctrail.low_complexity_fraction("AAAAAAAAAACCCCCCCCCC"), 0.5)
        self.assertEqual(lnctrail.low_complexity_fraction("ACGUACGUACGUACGU"), 0.0)

    def test_target_scan_keeps_one_row_per_input_pair_and_marks_no_hit(self):
        lncs = [lnctrail.SequenceRecord("q", "q TREX-like query", "ACGUACGUACGU")]
        targets = [
            lnctrail.SequenceRecord("t1", "t1 gene alpha (GENEA), mRNA", "ACGUACGUACGU"),
            lnctrail.SequenceRecord("t2", "t2 gene beta (GENEB), mRNA", "UUUUCCCCAAAA"),
        ]
        candidate = {
            "lncrna_id": "q", "lncrna_description": "q TREX-like query",
            "target_id": "t1", "target_description": targets[0].description,
            "engine": "LncTrail-builtin", "evidence_type": "sequence_heuristic",
            "lncrna_start": 1, "lncrna_end": 8, "target_start": 1, "target_end": 8,
            "integrated_priority_score": 75.0, "priority_score": 75.0,
            "target_gene": "GENEA", "binding_score": 20,
        }
        scan = lnctrail.build_target_scan([candidate], lncs, targets, "builtin")
        self.assertEqual(len(scan), len(lncs) * len(targets))
        status = {row["target_id"]: row["prediction_status"] for row in scan}
        self.assertEqual(status["t1"], "candidate_found")
        self.assertEqual(status["t2"], "no_candidate_passed_filters")
        no_hit = next(row for row in scan if row["target_id"] == "t2")
        self.assertEqual(no_hit["priority_score"], "")
        self.assertEqual(no_hit["target_gene"], "GENEB")

    def test_gene_summary_collapses_isoforms(self):
        targets = [
            lnctrail.SequenceRecord("NM_1.1", "NM_1.1 Homo sapiens gene alpha (GENEA), mRNA", "ACGUACGU"),
            lnctrail.SequenceRecord("NM_2.1", "NM_2.1 Homo sapiens gene alpha (GENEA), transcript variant 2, mRNA", "ACGUACGA"),
        ]
        rows = [
            {"lncrna_id":"q", "target_id":"NM_1.1", "target_gene":"GENEA", "integrated_priority_score":80,
             "priority_score":80, "priority_tier":"A", "engine":"LncTrail-builtin", "target_start":1, "target_end":8},
            {"lncrna_id":"q", "target_id":"NM_2.1", "target_gene":"GENEA", "integrated_priority_score":60,
             "priority_score":60, "priority_tier":"B", "engine":"LncTrail-builtin", "target_start":1, "target_end":8},
        ]
        summary = lnctrail.build_gene_summary(rows, targets)
        self.assertEqual(len(summary), 1)
        self.assertEqual(summary[0]["target_gene"], "GENEA")
        self.assertEqual(summary[0]["target_isoforms_in_input"], 2)
        self.assertEqual(summary[0]["transcripts_with_hits"], 2)
        self.assertEqual(summary[0]["best_transcript_id"], "NM_1.1")

    def test_intarna_batches_append_live_target_scan(self):
        """The main TSV must contain the previous completed batch when the next starts."""
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            query = root / "query.fa"
            targets = root / "targets.fa"
            outdir = root / "out"
            query.write_text(">q lncRNA\nACGUACGUACGUACGU\n", encoding="utf-8")
            targets.write_text(
                ">t1 gene one (GENE1), mRNA\nACGUACGUACGU\n"
                ">t2 gene two (GENE2), mRNA\nUUUUCCCCAAAA\n"
                ">t3 gene three (GENE3), mRNA\nGGGGAAAACCCC\n",
                encoding="utf-8")
            observations = {
                "calls": 0,
                "file_exists_before_first_batch": False,
                "rows_visible_before_second_batch": None,
            }

            def fake_intarna(lncs, target_batch, args, temp_dir, batch_index=None):
                observations["calls"] += 1
                live_file = outdir / "lnctrail_results.tsv"
                if observations["calls"] == 1:
                    observations["file_exists_before_first_batch"] = live_file.exists() and len(live_file.read_text(encoding="utf-8").splitlines()) == 1
                if observations["calls"] == 2:
                    observations["rows_visible_before_second_batch"] = len(live_file.read_text(encoding="utf-8").splitlines()) - 1
                return [], f"mock IntaRNA batch {batch_index}", []

            with patch.object(lnctrail, "first_available_engine", return_value="intarna"), \
                 patch.object(lnctrail, "run_intarna", side_effect=fake_intarna):
                code = lnctrail.main([
                    "--lncrna", str(query), "--targets", str(targets), "--outdir", str(outdir),
                    "--engine", "intarna", "--threads", "1", "--batch-size", "1", "--max-pairs", "10"
                ])

            self.assertEqual(code, 0)
            self.assertEqual(observations["calls"], 3)
            self.assertTrue(observations["file_exists_before_first_batch"])
            self.assertEqual(observations["rows_visible_before_second_batch"], 1)
            final_lines = (outdir / "lnctrail_results.tsv").read_text(encoding="utf-8").splitlines()
            self.assertEqual(len(final_lines) - 1, 3)
            metadata = __import__("json").loads((outdir / "run_metadata.json").read_text(encoding="utf-8"))
            self.assertEqual(metadata["status"], "completed")
            self.assertEqual(metadata["processed_targets"], 3)


if __name__ == "__main__":
    unittest.main()
