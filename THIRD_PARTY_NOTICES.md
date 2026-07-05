# Third-Party Notices

This repository is distributed as source code for OddSketch-Genome and its
benchmark workflows. Original OddSketch-Genome code is licensed under the MIT
License in `LICENSE`.

Some files, tools, and data sources are not covered by the MIT License for the
original OddSketch-Genome code. They are handled as follows.

## Bundled Source Code

### xxHash

- Files: `include/third_party/xxhash.hpp`,
  `include/third_party/xxhash_header_only.hpp`, `src/third_party/xxhash.cpp`
- Upstream: https://github.com/Cyan4973/xxHash
- License: BSD 2-Clause
- Copyright: Yann Collet and contributors

The copyright and BSD 2-Clause license notice is preserved in the vendored
source files.

### libpopcnt

- File: `include/third_party/libpopcnt.hpp`
- Upstream: https://github.com/kimwalisch/libpopcnt
- License: BSD 2-Clause
- Copyright: Kim Walisch, Wojciech Mula, and contributors

The copyright and BSD 2-Clause license notice is preserved in the vendored
source file.

## External Tools

### BinDash

BinDash source code and binaries are not included in this repository. BinDash is
used only for benchmark workflows that compare against the BinDash baseline.

- Upstream: https://github.com/zhaoxiaofei/bindash.git
- Recorded benchmark tag: `v2.6`
- Recorded benchmark commit:
  `ce2d16816beade65db992b8cd6eced00b54ca9ef`
- License reported by upstream: Apache License 2.0

If you install or redistribute BinDash separately, follow the BinDash upstream
license and citation requirements.

## External Data

### NCBI RefSeq

RefSeq genome FASTA files and the downloaded RefSeq assembly summary used in
the large benchmark are not included in this repository. The benchmark scripts
download RefSeq data from NCBI using:

```text
https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
```

RefSeq acquisition dates, URL construction, and integrity-check details are
documented in `experiments/refseq_sketch_task/README.md`.

NCBI states that it places no restrictions on use or distribution of molecular
data in its databases, but also notes that some original submitters or countries
of origin may claim intellectual-property rights and that NCBI cannot transfer
rights it does not hold. Users should follow NCBI policies and any applicable
third-party restrictions.

## Tutorial Fixtures

Files under `data/oddsketch_cli_sample/` are small synthetic tutorial fixtures
created for this repository. They are distributed with the original
OddSketch-Genome materials under the repository license, except where a file
explicitly states otherwise.
