# bioRxiv submission sheet

**Rejected at screening 11 September 2026.** The first attempt, BIORXIV/2026/750799, was
returned without review because the affiliation field was left as "None". bioRxiv requires every
submission to carry an organizational affiliation, so that some organization can adjudicate ethical
issues or disputes. Nothing about the manuscript itself was at issue.

**Resubmitting with Guidenco Inc. as the affiliation.** This is a new submission, not an appeal, so it gets
a new manuscript ID. Start it at https://submit.biorxiv.org. An acknowledgement email goes to
ronnel@onvo.ai, and the preprint appears publicly with a DOI if it passes screening, normally
within one to two working days.

Below is what to enter.

| Field | Value |
|---|---|
| Title | Per-variant tissue assignment by a sequence-to-function model is limited by effect detectability more than by assay coverage |
| Author | Ronnel Davis, sole author and corresponding author |
| Affiliation | Guidenco Inc., Delaware, USA |
| Email | ronnel@onvo.ai (account email and manuscript correspondence address, now consistent) |
| ORCID | to be supplied (or logged in through ORCID) |
| Subject category | Genomics |
| Licence | CC-BY 4.0 |
| Manuscript file | `main.pdf` (15 pages, 155 KB) |
| Abstract | 305 words, plain text below |
| Type | New results |
| Funding | None |
| Competing interests | None declared |
| Prior submission | Not under review at a journal |
| Data and code | https://github.com/ronneldavis/alphagenome-tissue-benchmark (cited in Data availability) |
| AI disclosure | Acknowledgements state that analysis was performed with assistance from Claude Code (Anthropic), and that design, claims and interpretation are the author's |

## Abstract as plain text

Sequence-to-function models predict tissue-resolved regulatory activity from DNA, so they can be asked in which tissue an individual non-coding risk variant acts. We tested whether AlphaGenome recovers trait-tissue relationships established by orthogonal evidence, scoring fine-mapped variants from 20 GWAS of 19 complex traits across 11 tissue groups against a shared background of 318 variants. The pre-specified expected tissue ranked first for 11 of 20 trait-GWAS pairs (P = 2.7e-7), within the top three for 16, and was significantly enriched for 14. Recovery did not track how well the expected tissue is measured: track count was uncorrelated with enrichment (rho = -0.01) or rank (rho = +0.25), although the test can exclude only a strong coverage effect. Enrichment did track a tissue-agnostic effect visibility (rho = +0.81), but largely by construction, because predicted histone effects are shared across tissues (mean inter-tissue correlation 0.80); visibility did not predict the expected tissue's rank (rho = +0.31, P = 0.19), and at the variant level tissue information was modest and concentrated in the most visible variants (within-trait odds ratio 1.45). Matched against LDSC-SEG on the same 16 GWAS and tissue definitions, the aggregate method was ahead (rank-1 recovery 10 versus 8; top three 15 versus 13). Its advantage was concentrated where per-variant effects are smallest: it placed brain first for both brain-expected traits (P = 0.002 and 2.4e-4) and the islet first for type 2 diabetes, where AlphaGenome ranked them third, third and 8th, whereas AlphaGenome recovered artery for coronary artery disease where LDSC-SEG found no significant tissue. On 2 of the 4 traits defeating both, the methods chose the same wrong tissue, implicating the labels. A null per-variant tissue assignment therefore more often reflects undetected effects than absent tissue involvement, and the tissue-specific component of predicted effects is small relative to the shared one.

## Notes

- bioRxiv screens every submission and may ask the corresponding author to verify their identity, which is more common for authors without a university address. A company domain email (ronnel@onvo.ai) helps here. Expect a delay of one to two working days before the preprint appears.
- bioRxiv will offer to forward the manuscript to a journal during submission. Decline it; pick the journal deliberately later.
