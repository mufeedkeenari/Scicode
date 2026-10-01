# Validation report - rbsog_npt_shear_pressure

Mechanical result: PASS

Semantic review is required to assess scientific test quality and whether the final orchestrator meaningfully uses earlier steps' outputs.

## Stage: structural - PASS

| check | description | result | findings |
|---|---|---|---|
| S1 | deliverables present/readable, no placeholder text | PASS | 0 |
| S2 | >=7 step files, NN_snake_case, sequential | PASS | 0 |
| S3 | one public function per step; stub is docstring + bare return | PASS | 0 |
| S4 | exactly one _oracle_<public> per step; no _gold_ prefix | PASS | 0 |
| S5 | Oracle imports/helpers survive Studio field and driver extraction | PASS | 0 |
| S6 | every _oracle_ name referenced in tests exists | PASS | 0 |
| S7 | oracle/helper unshadowed global calls do not use public step names | PASS | 0 |
| S8 | RNG discipline: no np.random.seed / legacy draws | PASS | 0 |
| S9 | stub and oracle fully annotated with matching signatures | PASS | 0 |
| S10 | isolated step factory: >=3 cases, setup/call/gold_call strings | PASS | 0 |
| S11 | except-hygiene in test setups | PASS | 0 |
| S12 | isolated integration factory: expression schema, literal value types, public final call | PASS | 0 |
| S13 | rubric schema: count, categories, weights, Browsing 15-25% | PASS | 0 |
| S14 | golden has the <reasoning>/<final_answer> tags the block requires, finite value | PASS | 0 |
| S15 | prompt hygiene: no paper identifiers, sentence count | PASS | 0 |
| S16 | browsing sources parse with real URLs and justifications | PASS | 0 |
| S17 | no harness implementation details in model-visible text | PASS | 0 |
| S20 | no filesystem or dynamic-import access in task code | PASS | 0 |
| S21 | every public function takes at least one parameter | PASS | 0 |
| S22 | same-named functions are identical across files | PASS | 0 |
| S24 | Studio section boundaries and documentation placement | PASS | 0 |
| S25 | each test side builds its inputs and helpers from its own implementations | PASS | 0 |

## Stage: execute - PASS

| check | description | result | findings |
|---|---|---|---|
| E1 | all step cases pass twice (differential + determinism) | PASS | 0 |
| E2 | no step with missing/broken tests silently skipped | PASS | 0 |
| E3 | every result contains only supported comparable values | PASS | 0 |
| E4 | wall-clock within thresholds | PASS | 0 |
| E5 | integration expected-value/property comparisons execute and pass | PASS | 0 |
| E6 | golden <final_answer> matches the final integration case | PASS | 0 |
| E7 | golden / integration / computed triple agrees | PASS | 0 |

## Stage: precalibration - PASS

| check | description | result | findings |
|---|---|---|---|
| S23 | advisory only: numeric overlap with the golden final answer | PASS | 0 |
| P2 | cases produce distinct results (normal/boundary/edge) | PASS | 0 |
| P35 | local comparator acceptance bounds, including small signals and relative tolerance | PASS | 0 |
| P3 | invalid-input coverage (ValueError contract exercised) | PASS | 0 |
| P4 | syntactic exception evidence for mandatory LLM contract review P31 | REVIEW REQUIRED | 0 |
| P5 | explicit oracle return-site coverage only; P16 LLM coverage review still required | PASS | 0 |
| P6 | one docstring convention across step files, per surface | PASS | 0 |
| P8 | potential external-reference wording for LLM review P32 | PASS | 0 |
| P11 | possible code-call syntax for LLM review P32 | PASS | 0 |
| P12 | potential contrastive wording for LLM review P32 | PASS | 0 |
| P26 | runtime observations for mandatory LLM state/dependency review P30 | REVIEW REQUIRED | 0 |
| P27 | no step passable by a constant-return function | PASS | 0 |
| P34 | wording to assess for qualification in LLM review P32 | PASS | 0 |

## Required semantic review

REVIEW REQUIRED is a pending LLM assessment, not a mechanical failure or a semantic PASS. Evidence is retained even when only structural validation was requested.

### P16 (S12/E5/E6/E7) - REVIEW REQUIRED

Assess per-step test coverage, model/oracle dependencies and the integration tests under SKILL.md's Representative test coverage, Model and oracle dependencies and Integration tests sections. Record evidence or verification limits; mechanical agreement here does not complete that review.

Evidence collection: COLLECTED.

- `{"case": 1, "call": "compute_pressure_relative_error(positions, charges, cell, 1.6, 3.0, 10, 9.0, 6, 128, 0, 1)", "gold_call": "0.37851954139827876", "selected_final": true, "actual": "0.37851954139827876", "expected": "0.37851954139827876", "comparison_passed": true, "benchmark_extraction": "direct scalar comparison"}`
- `{"case": 2, "call": "compute_pressure_relative_error(positions, -charges, cell, 1.6, 3.0, 10, 9.0, 6, 128, 0, 1)", "gold_call": "0.37851954139827876", "selected_final": false, "actual": "0.37851954139827876", "expected": "0.37851954139827876", "comparison_passed": true}`
- `{"case": 3, "call": "compute_pressure_relative_error(shifted, charges, cell, 1.6, 3.0, 10, 9.0, 6, 128, 0, 1)", "gold_call": "0.37851954139827876", "selected_final": false, "actual": "0.3785195413982788", "expected": "0.37851954139827876", "comparison_passed": true}`
- `{"case": 4, "call": "compute_pressure_relative_error(positions, charges, cell, 2.0, 4.0, 6, 8.0, 6, 64, 1, 2)", "gold_call": "0.06056888896760986", "selected_final": false, "actual": "0.06056888896760986", "expected": "0.06056888896760986", "comparison_passed": true}`

### P30 (P26) - REVIEW REQUIRED

Assess explicit state, task dependencies and self-containment from the task source and runtime observations. Mutation of supplied arrays, RNGs, mutable objects, output/work buffers and views is permitted, as are internal caches/work state and scoped numerical settings restored on exit. Cache history must not change the required result. State carried between calls that affects the required scientific result must be supplied explicitly. Check that mutation or preservation is documented and tested when the contract relies on it, and that candidate/oracle comparisons start from separate equivalent mutable inputs. Attribute an operation to task behavior or dependency initialization before identifying a defect. Dependency initialization alone is not a task defect; no observed events does not establish contract compliance. The audit is not an execution sandbox.

Evidence collection: COLLECTED.

No observations collected; semantic review is still required.

### P31 (P4) - REVIEW REQUIRED

Determine which exceptions are tested or can escape for documented inputs and whether the public function's docstring states that contract, consistently with the scientific background. New Studio uploads put parameter, return and exception documentation in the public function's docstring. Raise ValueError for invalid data actually supplied by the task’s tests or function calls, including upstream calls reached from valid downstream inputs. Functions may rely on established preconditions. Validation is required only for invalid data passed by tests included in the task. Additional validation solely for hypothetical invalid calls is not required. Invalid-input tests are optional. Syntactic raise/except references may be caught, unused, or incidental; name matching is not a contract verdict.

Evidence collection: COLLECTED.

- `{"step": "01_compute_pressure_kernel_gaussians.py", "oracle_raise_names": ["ValueError"], "setup_except_names": ["ValueError"], "names_absent_from_visible_docstrings": []}`
- `{"step": "02_compute_narrowest_weight_factor.py", "oracle_raise_names": ["ValueError"], "setup_except_names": ["ValueError"], "names_absent_from_visible_docstrings": []}`
- `{"step": "03_compute_short_range_pressure.py", "oracle_raise_names": ["ValueError"], "setup_except_names": ["ValueError"], "names_absent_from_visible_docstrings": []}`
- `{"step": "04_compute_structure_factor_power.py", "oracle_raise_names": ["ValueError"], "setup_except_names": ["ValueError"], "names_absent_from_visible_docstrings": []}`
- `{"step": "05_compute_long_range_pressure.py", "oracle_raise_names": ["ValueError"], "setup_except_names": ["ValueError"], "names_absent_from_visible_docstrings": []}`
- `{"step": "06_compute_nonradial_normalization.py", "oracle_raise_names": ["ValueError"], "setup_except_names": ["ValueError"], "names_absent_from_visible_docstrings": []}`
- `{"step": "07_compute_nonradial_variance.py", "oracle_raise_names": ["ValueError"], "setup_except_names": ["ValueError"], "names_absent_from_visible_docstrings": []}`
- `{"step": "08_compute_pressure_relative_error.py", "oracle_raise_names": ["ValueError"], "setup_except_names": ["ValueError"], "names_absent_from_visible_docstrings": []}`


## Stage: precalibration-semantic - PASS

Evidence: a fresh Stage 3 review agent (writer pre-review, report-only) followed the reviewer skill's task-review instructions. A configuration audit, a paper-based blind solve, a code-based blind solve and a fresh rubric-compliance agent ran as isolated agents. The orchestrator verified every finding against the task, the paper (arXiv:2602.23582v1, all 29 pages) and independent recomputation. Independent results:

- Both blind solves reproduced the benchmark exactly, 0.3785195413982788.
- Two earlier blind re-implementations made during authoring agree to about 1e-15.
- An independent Ewald evaluation agrees with the SOG P_xy to 0.27 %, the decomposition error at b = 1.6. The gap vanishes as b → 1.

| check | result | finding |
| --- | --- | --- |
| P13 | PASS | The prompt and background identify the method without coined terms and state no Browsing-graded formula. All other content is instance configuration: the truncated C0 rule, the all-images rule, the i.i.d. idealization, P = 128. No graded value appears. |
| P14 | PASS | Criteria 6–8, 10 and 11 are values the prompt requests. Criteria 1–5 and 9 are method content or on-path derivations, all of which the paper solver computed. |
| P15 | PASS | All instance data is supplied. The variants that matter (truncated vs untruncated ω̃, minimum image, bound vs exact variance, radial vs non-radial proposal, correlated MH chain) are resolved by the prompt or by the method. The defensible alternatives (untruncated ω̃, dividing by the Ewald value) stay within tolerance. |
| P16 | PASS | `call` chains use public names and `gold_call` chains use oracle names. Each step has normal, boundary and edge cases representative of downstream use. Comparator allowances admit only the untruncated-ω̃ variant, which step 02 grades separately. There are four integration configurations: the literal benchmark, charge conjugation, rigid translation, and an independent 2:1 salt literal. None is differential. |
| P17 | PASS | No duplicate criteria. |
| P18 | PASS | RSE 8, variance 7, non-radial proposal 5, P_xy parts 5. Bookkeeping items weigh 2. |
| P19 | PASS | Bound-as-equality and missing-rescaling responses earn justified partial credit. Equivalent forms, other units and the row-convention wavevectors are accepted. |
| P20 | WARN | Criteria 1–3 are supported in the source (Eqs. 3.1, 3.19–3.22, 3.17/3.21), and the query is realistic. Calibration watch: criterion 3, and possibly criterion 1, may be answerable without browsing. Check against no-browsing responses at Stage 4. If criterion 3 is reclassified, the Browsing share becomes 8/48 = 16.7 %. |
| P21 | PASS | The golden has 8 numbered steps, states every graded value and carries the correct tags. The paper's Prop. 7 bound (1/36, an inequality) is handled correctly: the exact 1/64 is used. |
| P22 | PASS | 7 of 8 steps (87.5 %) need the primary paper's method. Step 04 (structure factor) is generic. |
| P23 | WARN | Step text states conventions only. Step 05's radial/non-radial split definition is the paper's (Eqs. 3.17–3.18). It is kept as the return-split convention and discloses no coefficient. Judgment call for human review. |
| P24 | PASS | Notation, units and symbols agree across prompt, golden, rubric and code. |
| P25 | PASS | No typos or placeholders. The code compiles under Python 3.12. |
| P28 | PASS | No graded value appears in reasoning-stage or code-stage text (manual numeric scan). |
| P30 | PASS | Inputs are copied and never mutated, and tests pass `.copy()`. There is no RNG, I/O or global state. Same-named helpers are identical across files. |
| P31 | PASS | All 18 invalid cases raise the documented ValueError from the intended check. Every graded exception is documented in its public docstring. |
| P32 | PASS | The first compliance pass found F1–F5 (criteria 1, 2, 4, 9, 10). All were verified and applied. Re-check 1 found three more wording issues (criteria 5, 8, 10, 11), which were applied. Re-check 2: compliant. Scientific values, tolerances and near-miss claims were reproduced independently. Records are kept with the review evidence (out_rubric, recheck, recheck2). |
| P33 | PASS | The steps cover the SOG pressure split, C0 rescaling, real space, Fourier space, the radial/non-radial split, the non-radial proposal, its normalization, the exact variance and the batch error. The MH acceptance rule (Eq. 3.25) is outside the requested ideally-mixed quantity. |

Readiness: ready to paste into Studio for calibration. At Stage 4, check the Browsing classification of criteria 3 and 1 against the no-browsing responses.
