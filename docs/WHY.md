# Why it's built this way

## The problem

Every sample in a genomic data release needs the same checks: contamination, sex, coverage, and whether it matches the same person on another platform. The metrics come from different tools, and inconsistent hand-made decisions are how a swapped or contaminated sample slips through.

## Design choices

**Why read tool outputs instead of computing metrics?** Because VerifyBamID, PLINK and mosdepth already compute them well, and reimplementing them would add risk without adding value. The useful part is applying the same thresholds the same way every time.

**Why separate readers from checks?** Because file formats change more often than QC rules. A reader turns a file into a dict; a check turns a dict into a list of failures. Each is tested on its own.

**Why compare genotypes as unordered allele pairs?** Because platforms write the same genotype differently (`AG`, `G/A`, `A|G`). String comparison would call these discordant and make a good sample look swapped.

**Why keep `0/0` separate from missing?** Because in VCF it means homozygous reference. Treating it as missing would quietly drop the most common genotype from the comparison.

**Why report samples with too little overlap?** Because "we could not check this sample" is different from "this sample passed". Hiding the first inside the second is how untested samples get released.

**Why PLINK's F thresholds for sex?** Because they are what PLINK prints as OK or PROBLEM, so a reviewer comparing this tool's output with PLINK's sees the same answer.

## Questions worth asking

**"Is 5% FREEMIX the right contamination threshold?"**
It is a common starting point, not a rule. The right value depends on the downstream use (somatic calling is far less tolerant than common-variant GWAS) and on how FREEMIX was estimated. It is a parameter because it should be a decision.

**"How would you tell a sample swap from a strand problem?"**
A swap affects one or a few samples: their concordance drops to what unrelated people share, while everyone else stays high. A strand or reference-allele problem affects every sample and is concentrated in particular variants. Looking at per-sample rates (kept in `per_sample_rate`) and per-variant discordance separates the two; the second is on the roadmap.

**"What about a sample with an intermediate F value?"**
Values between 0.2 and 0.8 are ambiguous. They can mean contamination, aneuploidy (for example XXY) or poor data. The check fails a sample whose F does not match its reported sex, but an intermediate F needs a person to look at it; the tool's job is to make sure someone does.

## What's next

- Per-variant discordance to separate strand problems from swaps.
- A statistical batch-effect test.
- One combined per-sample QC table for release sign-off.
