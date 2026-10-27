# LncTrail: landscape review and method selection

This review compares methods relevant to the specific use case: lncRNA FASTA + target-transcript FASTA → ranked RNA–RNA interaction candidates. It does not imply that every lncRNA mechanism can be predicted from sequence alone.

## Summary comparison

| Method/resource | Main signal | Accessibility modeled? | Scale/use case | Role in LncTrail v0.2.3 |
|---|---|---:|---|---|
| **IntaRNA 3.x** | Thermodynamic hybridization and seed constraints | Yes | High-throughput all-vs-all with native threading | Integrated thermodynamic backend |
| Built-in LncTrail engine | Exact reverse-complement seed plus contiguous extension | No | Fast exploratory screening and software smoke tests | Integrated fallback; heuristic score only |
| RNAplex (ViennaRNA) | Efficient target-site hybridization | May use accessibility profiles | Fast site search | Useful independent cross-check; not integrated |
| RNAup (ViennaRNA) | Binding energy plus site-opening cost | Yes | Focused transcript/site analysis | Orthogonal thermodynamic cross-check; not integrated |
| RIsearch2 | Seed search and energy-based extension | Not equivalent to full accessibility modelling | Large-scale screening; reusable target index | Potential alternative for large target sets; not integrated |
| RIblast | Fast seed-and-extension interaction prediction | Approximation | Transcriptome-scale screening | Potential alternative; not integrated |
| RIME | Deep RNA-sequence-model embeddings over windows | Learned model representation | Model-heavy local or web workflow | Optional shortlist rescoring through LncTrail; disabled by default |
| LncTar | Normalized interaction free energy | Joint RNA structure model | Candidate lncRNA–mRNA pairs | Historical comparator; not integrated |
| LongTarget / Triplexator | RNA–DNA triplex motifs/sites | Different modality | Genomic DNA plus lncRNA motif search | Requires a separate DNA-triplex module |
| LncTarD | Experimentally supported lncRNA–target relationships | Not a de novo predictor | Literature evidence and biological context | Potential future evidence-annotation source |
| catRAPID and related resources | RNA–protein interaction propensity | Different input/features | Protein-binding hypotheses | Requires a separate RNA–protein module |

## Why IntaRNA is supported

IntaRNA supports multiple query and target FASTA records, predicts interaction sites, reports coordinates and thermodynamic terms, models accessibility of interacting regions, accepts seed constraints, and provides native multithreading. These features align with LncTrail's input/output use case. See the [official documentation](https://backofenlab.github.io/IntaRNA/).

There is no universally optimal seed length or free-energy cutoff for eukaryotic lncRNA interactions. Test parameter settings against known interactions and appropriate sequence-matched decoys before interpreting any threshold as evidence.

## Why include a built-in engine?

The built-in engine allows installation and smoke testing without an external binary. It indexes exact reverse-complement seeds (12 nt by default in this release), extends contiguous alignments, and computes a transparent heuristic score. It returns at most one site per transcript pair by default and uses seed-specificity and low-complexity features in its ranking. It does **not** calculate thermodynamic free energy, structural accessibility, or a probability. It is for rapid triage and testing, not a replacement for thermodynamic modelling.

## Ranking and optional deep-learning evidence

LncTrail produces both site-level tables and a gene-level table that attempts to collapse RefSeq transcript isoforms using symbols inferred from FASTA descriptions. Verify that grouping against the matching genome assembly and reference annotation. `priority_score` is a heuristic ranking index, not a validated probability or p-value. It should be used to order candidates within a run and followed by independent biological assessment.

The optional RIME integration rescores a limited shortlist using a published deep-learning model that evaluates sequence windows. It needs a separate repository, model/checkpoint and environment. RIME's output, IntaRNA's energy, and LncTrail's heuristic score represent different evidence and must not be compared as if they share one scale. Check the [RIME repository](https://github.com/tartaglialabIIT/RIME) for current software requirements and license terms.

## Biological limitations and recommended follow-up

1. **Distinct mechanisms need distinct models.** Direct RNA–RNA pairing, RNA–DNA triplex formation, RNA–protein binding, chromatin-mediated regulation and miRNA competition are different biological processes.
2. **Input sequences matter.** Use correctly oriented spliced transcript sequences with relevant UTRs and appropriate isoforms. Genomic DNA and intron-containing loci are not suitable substitutes for mature RNA transcripts in RNA–RNA pairing prediction.
3. **Transcript-level results are not gene-level proof.** Isoforms may have different binding regions. Confirm transcript-to-gene mappings against an assembly-matched annotation.
4. **Predicted energy is not biological function.** A favorable duplex cannot prove expression co-occurrence, subcellular colocalization, in-vivo binding, regulation or phenotype.
5. **Prioritize with independent evidence.** Consider matched cell/tissue expression, localization, relevant interaction assays (e.g. RAP, CHIRP, CHART or CLASH as appropriate), perturbation/rescue results and carefully curated publications.
6. **Do not mix scores across methods.** IntaRNA energy, RNAplex/RNAup energies, RIME scores and LncTrail heuristic scores have different meanings.

## References checked

- IntaRNA documentation: https://backofenlab.github.io/IntaRNA/.
- IntaRNA Bioconda recipe: https://bioconda.github.io/recipes/intarna/README.html.
- RNAplex documentation: https://viennarna.readthedocs.io/en/latest/man/RNAplex.html.
- RIsearch2 software/manual: https://github.com/RTH-tools/risearch/ and https://rth.dk/resources/risearch/.
- RIblast source: https://github.com/fukunagatsu/RIblast.
- RIME paper: https://doi.org/10.1038/s41467-026-68356-w; model: https://github.com/tartaglialabIIT/RIME; web server: https://tools.tartaglialab.com/rna_rna.
- LncTar paper: https://doi.org/10.1093/bib/bbu048.
- lncRNA interaction benchmarking: https://doi.org/10.1093/bib/bby032.
- LncTarD 2.0: https://lnctard.bio-database.com/.
- LongTarget: https://github.com/LongTarget/LongTarget.

## Incremental outputs in v0.2.3

LncTrail writes progress results after each target batch rather than relying on an IntaRNA process to flush partial standard output. Batch size is configurable; a small batch gives more frequent progress updates but adds IntaRNA process startup overhead. Final global ranking is calculated once after all batches are complete.
