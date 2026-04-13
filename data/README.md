# Data

Pipeline data managed by [DVC](https://dvc.org/). Run `pixi run repro` to reproduce.

## Directory structure

| Directory | Contents |
|-----------|----------|
| `metadata/` | Curated metadata CSVs (raw, augmented, cleaned) and variant definitions |
| `3di/` | 3DI (structural alphabet) sequences and trees |
| `aa/` | Amino acid sequences and trees |
| `models/` | 3DI substitution matrices (Q.3DI.AF, Q.3DI.LLM, VK23) |
| `cif/` | Downloaded mmCIF structure files (gitignored, ~2GB) |
| `foldseek/` | Intermediate foldseek database files (gitignored) |
