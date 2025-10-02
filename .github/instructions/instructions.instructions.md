---
applyTo: '**'
---
Provide project context and coding guidelines that AI and contributors should follow when generating code, answering questions, or reviewing changes.

Project: Multi-object spectroscopy with DMDs (student project)

Overview
--------
This repository contains a student research project exploring the prospects and performance of multi-object spectroscopy for the amateur-astronomy community, using Digital Micro-mirror Devices (DMDs) as a variable slit mask. The project is implemented primarily in Python and organized as Jupyter notebooks for analysis and demonstration, with supporting Python modules for reusable helper functions.

Primary goals
-------------
- Implement simulation and data-analysis pipelines that model light collection, spectrometer throughput, noise sources, and extraction for many objects simultaneously.
- Use realistic optics and detector models to estimate signal-to-noise, throughput, and limiting magnitudes.
- Produce reproducible notebooks and helper modules so others can replicate and extend the work.

Language, libraries and style
----------------------------
- Language: Python 3.11+ (where possible)
- Preferred libraries:
  - numpy for numerical work
  - astropy for astronomical unit handling, coordinates, tables, and I/O
  - pandas for tabular data where appropriate
  - matplotlib / seaborn for plotting
  - jupyter / ipykernel for interactive notebooks
  - pytest for tests
  - optionally: numba for hotspots that need JIT speed-ups
- Use type hints and small, well-documented functions. Prefer numpy vectorized operations over Python loops when possible.
- 

Repository layout (recommended)
------------------------------
- `spectrometer-calculation.ipynb` — primary analysis notebook(s)
- `notebooks/` — additional experiment or exploration notebooks
- `src/utils` — python helper functions and modules
- `data/` — raw and derived data (large raw files should be stored externally and referenced via README)
- `.github/instructions/` — this instructions file
- `environment.yml` — conda environment to reproduce the environment

Environment & reproducibility
-----------------------------
We use a conda environment to pin dependencies and reproduce results. A sample `environment.yml` is present at the repository root. To create the environment locally:

```bash
# create the environment
conda env create -f environment.yml

# activate it
conda activate s4-env

# register a Jupyter kernel so notebooks can use it
python -m ipykernel install --user --name s4-env --display-name "Python (s4-env)"
```

Notes:
- If you update `environment.yml`, update it in a single commit and include a short rationale in the commit message. To update an existing environment:

```bash
conda env update -f environment.yml --prune
```

Data handling and large files
----------------------------
- Keep small example datasets inside `data/` for examples and tests.
- For large observational or simulated datasets, provide download scripts or external links and a small script to fetch and verify checksums.
- Never commit sensitive or private data.

Notebooks: conventions and best practices
---------------------------------------
- Use notebooks for descriptions, visual analysis, and demonstrations. Put production-ready, reusable code into `src/` modules and import them from notebooks.
- Each notebook should have a short header cell describing purpose and required environment (e.g., kernel name).
- Avoid leaving long-running cells (e.g., heavy simulations) in master branches; instead include checkpoints or sample results and provide scripts to reproduce heavy runs offline.
- Prefer explicit imports at the top of each notebook cell (or at least a well-defined "environment" cell that imports everything used).
- Use `astropy.units` when working with physical quantities to avoid unit mistakes.

Coding conventions and tests
---------------------------
- Follow PEP8-style formatting. Keep functions small and single-purpose.
- Add docstrings to all public functions with argument and return descriptions.
- Code should be keept minimal, without elaborate error handling or abstractions unless necessary.

```bash
pytest -q
```

Quality gates (recommended before merging changes)
-----------------------------------------------
- Build: ensure notebooks run to their final state using the pinned environment (optionally use nbconvert or nbclient for CI).
- Lint/Typecheck: run a flake8 or ruff check and mypy if types are used.
- Unit Tests: all tests pass.

Edge cases to consider when implementing simulations
---------------------------------------------------
1. Zero or missing flux in input spectra or extreme faint sources.
3. Detector issues: hot pixels, saturation, non-linear response.
4. Telescope pointing errors and partial object coverage of DMD micro-mirrors.

Running notebooks in VS Code or JupyterLab
-----------------------------------------
- After creating and activating the conda env and installing the ipykernel, open the notebook in VS Code or JupyterLab and select the kernel named `Python (s4-env)` or `s4-env`.
- If the kernel/ interpreter doesn't appear, restart VS Code, ensure the Python extension is installed, and confirm `ipykernel` is installed in the environment.

Quick troubleshooting
---------------------
- `conda env list` — list available envs.
- `jupyter kernelspec list` — list installed kernels.
- If imports fail inside the notebook but work in the terminal, ensure the notebook kernel points to the same interpreter as `conda run -n s4-env python --version`.

Security and sharing
--------------------
- Avoid embedding API keys, credentials, or personal data in notebooks. Use environment variables or external config files not committed to the repo.

References and further reading
------------------------------
- astropy: https://www.astropy.org/
- numpy: https://numpy.org/

Contact
-------
For questions about the project or contributions, open an issue in this repository.

Acceptance criteria for AI-generated edits
-----------------------------------------
- Use Python and the preferred libraries when possible.
- Keep changes minimal and focused; follow existing module and notebook structure.
- Add or update unit tests for any new or modified numerical functionality.
- Prefer reproducible outputs and avoid running or committing heavy binary outputs.