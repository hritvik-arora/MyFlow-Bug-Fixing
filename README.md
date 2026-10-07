# Investigation of Grain Counting Inconsistencies in MYflow

A MATLAB-based investigation into grain-counting discrepancies observed during composite phase-map analysis in **MYflow**.

This repository contains the modified `grainsize.m` implementation developed during the investigation, along with the complete technical report documenting the observations, debugging process, code-level changes, and results.

---

## Overview

During testing of **MYflow**, the software was initially evaluated using the example dataset provided with the program. The example ran successfully and produced results consistent with the expected output.

However, when custom composite phase maps were introduced, an inconsistency was observed between the number of grains that could be identified visually and the number of grains reported by MYflow.

In one of the test cases, the K-feldspar phase contained **four visually distinguishable grains**, while the original implementation reported **seven grains**.

This led to a detailed investigation of the composite-mode grain-identification workflow, particularly the way RGB phase maps were processed before grain identification and connected-component analysis.

---

## Repository Contents

```text
MYflow-Grain-Counting-Investigation/
│
├── grainsize.m
├── Investigation of Grain Counting Inconsistencies in MYflow.pdf
└── README.md
