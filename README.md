# Cosmochlore

**COSM**ochlo**R**e (**Co**ntinuous **S**hape **M**easurements in **R**ust) is a fast, pure-Rust implementation of the Continuous Shape Measures (CShM) Calculation, used to quantify how closely a set of points matches an idealized reference polyhedron.

It is a from-scratch command line interface (CLI) tool that reimplements the shape-measure engine found in [`cosymlib`](https://github.com/GrupEstructuraElectronicaSimetria/cosymlib) and `SHAPE`<sup>1</sup> 2.1 in Rust. Cosmochlore accurately reproduces `SHAPE`'s 2.1 results using a pruned branch-and-bound algorithm for significantly faster performance on larger coordination numbers. One of the advantages of Rust over the old Fortran code is that Cosmochlore's error handling will always tell the user if something went wrong at run-time, the program will never silently crash or give you a number without you knowing something went wrong. 

The name of the tool comes from the mineral [Kosmochlor](https://en.wikipedia.org/wiki/Kosmochlor), a rare chromium clinopyroxene found in iron meteorites and as an accessory mineral to various other chromium-containing pyroxenes.

## What it does

Given a molecular structure (or any set of 3-dimensional points), Cosmochlore computes the **Continuous Shape Measure (CShM)** between the structure and one or more idealized reference shapes. A value of `CShM = 0` means a perfect match; larger values indicate greater distortion from the ideal geometry.

### Future Features
I'm planning to support Cosmochlore in [SymmetryMeasurements](https://github.com/Yluro/symmetry-measurements/tree/master) in the very near future. More features are planned to be added:
 - Finishing the CShM toolkit: retrieval of the rotation matrix, generalized shape coordinate, _etc_.
 - Might do some octahedral distortion parameters.
 - Might reimplement the continuous symmetry/symmetry operation measures.

_[See you in 25 years...](https://www.youtube.com/watch?v=BL57-9171pk)_

## Installation
**Build from source:**
````bash
git clone https://github.com/Yluro/cosmochlore
cd cosmochlore
cargo build --release
````
`cargo` will compile the source code taking into account your system's architecture. The easiest way to get `cargo` is to install the latest stable release of [Rust](https://rust-lang.org/) using [`rustup`](https://rustup.rs/).

The compiled binary will be at `target/release/cosmochlore`. By cloning the current main branch of the repository you will get access to experimental versions of Cosmochlore.

### Precompiled binary downloads:

Precompiled binaries for the most common operating systems are available in this repository's [releases](https://github.com/Yluro/cosmochlore/releases).

_Note: It is strongly recommended to place the cosmochlore executable in the systems `PATH`._

## Usage

Cosmochlore is a CLI tool. It can be called by the general syntax:

```bash
cosmochlore <COMMAND> <NAME> [OPTIONS]
```
### Commands

| Command | Description |
|---|---|
| `cshm` | Continuous shape measures calculation. |
| `csom` | Continuous symmetry operation measures calculation. _Unimplemented as of v0.2.0._ |
| `help` | Print the help or the help of the given subcommand(s). |

### Arguments:
| Flag | Description |
|---|---|
| `<NAME>` | **Required**. Path to the input `.xyz` file containing the structure to analyze. |

#### Optional arguments for `cshm`:

| Flag | Value | Description |
|---|---|---|
| `-n` <br> `--nc` | None | Indicates the structure does **not** include an explicit central atom (i.e. only ligand/vertex coordinates are given). If a structure contains a central atom, Cosmochlore assumes it is in the first position of the `.xyz` file.  |
|`-s` <br> `--sh` | `<SHAPES>...` | Restrict the comparison to specific built-in reference shapes by index, for the detected vertex count. If omitted, all applicable built-in shapes are used. Specified indices should be separated by whitespace. |
| `-r` <br> `--ref` | `<USER_SHAPES>...` | Path to the `reference.yaml` files that contain user-defined shapes to include in the CShM calculation. Specified files should be separated by whitespace.
|`-t` <br> `--table`| None | Write a `name_table.csv` file containing the output of the calculation. |
|`-i` <br> `--ideal`| None | Write a `name_ideal.xyz` file containing the reconstructed idealised structures for **all** the selected reference shapes. |

_Oh..., I lost my crab here. 🦀 Thanks for finding it!_

### Example usage of `cshm`
Given the following `FeCl6.xyz` file:
````xyz
7
High-spin iron(ii) complex
Fe6 4.92991 10.3899 12.9237
N004 6.20468 10.6922 14.7747
N006 5.00034 8.51382 13.9503
N008 6.80263 10.0749 11.7725
N00B 5.47238 12.2809 12.1201
N00C 3.30536 11.3249 13.862
N00D 3.74731 9.47729 11.4687
````
Running:

````cmd
cosmochlore cshm FeHS.xyz --ref ebcT-6.yaml -t -i
````

Will output:

````output
  _  ______   _____ __  __  ____   _____ _    _ _      ____  _____
 | |/ / __ \ / ____|  \/  |/ __ \ / ____| |  | | |    / __ \|  __ \
 | ' / |  | | (___ | \  / | |  | | |    | |__| | |   | |  | | |__) |
 |  <| |  | |\___ \| |\/| | |  | | |    |  __  | |   | |  | |  _  /
 | . \ |__| |____) | |  | | |__| | |____| |  | | |___| |__| | | \ \
 |_|\_\____/|_____/|_|  |_|\____/ \_____|_|  |_|______\____/|_|  \_\

Continuous Shape Measurements in Rust.
Version: 0.2.0
Authors: José Serrano Guarinos <jose.serranog@ub.edu>
Repository: https://github.com/Yluro/cosmochlore
============================================================
 Symbol   Shape                           Symmetry   CShM   
------------------------------------------------------------
 HP-6     Hexagon                         D6h        33.215 
 PPY-6    Pentagonal pyramid              C5v        23.015 
 OC-6     Octahedron                      Oh         2.109  
 TPR-6    Trigonal prism                  D3h        11.027 
 JPPY-6   Johnson pentagonal pyramid J2   C5v        27.034 
 ebcT-6   Edge-bicapped tetrahedron       D2d        14.335 
------------------------------------------------------------
Writing output table to .\tests\FeHS_table.csv...
Writing idealised polyhedra coordinates to table to .\tests\FeHS_ideal.xyz...
Program finished in 18.5727ms
````
The output of the calculation is saved in the file `FeHS_table.csv` by calling the `-t` or `--table` flag.

````csv
Symbol,Name,Symmetry,CShM
HP-6,Hexagon,D6h,33.215
PPY-6,Pentagonal pyramid,C5v,23.015
OC-6,Octahedron,Oh,2.109
TPR-6,Trigonal prism,D3h,11.027
JPPY-6,Johnson pentagonal pyramid J2,C5v,27.034
ebcT-6,Edge-bicapped tetrahedron,D2d,14.335
````

The `FeHS_ideal.xyz` file was produced by calling the `-i` or `--ideal` flag. The coordinates of the ideal octahedron placed in the correct position of the structure can be extracted from it.

<img width="411" height="355" alt="ideal reference octahedron superimposed with the problem shape" src="https://github.com/user-attachments/assets/490458e7-a4a5-4c0b-b82d-3444ace939e6" />

_Note that ebcT-6 is a non-standard reference shape included by passing the `--ref` flag. See more below._

## Reference Polyhedra
The geometries of 90 reference polyhedra are internally defined in Cosmochlore. This list was integrally derived from the `SHAPE` 2.1 list of reference polyhedra and has been discussed in numerous articles by Alemany, Llunell, Alvarez, Avnir, Cirera _et at._<sup>2</sup>

<table>
<thead>
<tr>
<th>Vertices</th>
<th>Index</th>
<th>Label</th>
<th>Shape</th>
<th>Symmetry</th>
</tr>
</thead>
<tbody>
<tr><td rowspan="3">2</td><td>0</td><td>L-2</td><td>Linear</td><td>D<sub>∞h</sub></td></tr>
<tr><td>1</td><td>vT-2</td><td>Divacant tetrahedron (V-shape, 109.47º)</td><td>C<sub>2v</sub></td></tr>
<tr><td>2</td><td>vOC-2</td><td>Tetravacant octahedron (L-shape, 90º)</td><td>C<sub>2v</sub></td></tr>

<tr><td rowspan="4">3</td><td>0</td><td>TP-3</td><td>Trigonal planar</td><td>D<sub>3h</sub></td></tr>
<tr><td>1</td><td>vT-3</td><td>Pyramid‡ (vacant tetrahedron)</td><td>C<sub>3v</sub></td></tr>
<tr><td>2</td><td>fac-vOC-3</td><td>Trivacant octahedron</td><td>C<sub>3v</sub></td></tr>
<tr><td>3</td><td>mer-vOC-3</td><td>Trivacant octahedron (T-shape)</td><td>C<sub>2v</sub></td></tr>

<tr><td rowspan="4">4</td><td>0</td><td>SP-4</td><td>Square</td><td>D<sub>4h</sub></td></tr>
<tr><td>1</td><td>T-4</td><td>Tetrahedron</td><td>T<sub>d</sub></td></tr>
<tr><td>2</td><td>SS-4</td><td>Seesaw or sawhorse‡ (cis-divacant octahedron)</td><td>C<sub>2v</sub></td></tr>
<tr><td>3</td><td>vTBPY-4</td><td>Axially vacant trigonal bipyramid</td><td>C<sub>3v</sub></td></tr>

<tr><td rowspan="5">5</td><td>0</td><td>PP-5</td><td>Pentagon</td><td>D<sub>5h</sub></td></tr>
<tr><td>1</td><td>vOC-5</td><td>Vacant octahedron‡ (Johnson square pyramid, J1)</td><td>C<sub>4v</sub></td></tr>
<tr><td>2</td><td>TBPY-5</td><td>Trigonal bipyramid</td><td>D<sub>3h</sub></td></tr>
<tr><td>3</td><td>SPY-5</td><td>Square pyramid§</td><td>C<sub>4v</sub></td></tr>
<tr><td>4</td><td>JTBPY-5</td><td>Johnson trigonal bipyramid (J12)</td><td>D<sub>3h</sub></td></tr>

<tr><td rowspan="5">6</td><td>0</td><td>HP-6</td><td>Hexagon</td><td>D<sub>6h</sub></td></tr>
<tr><td>1</td><td>PPY-6</td><td>Pentagonal pyramid</td><td>C<sub>5v</sub></td></tr>
<tr><td>2</td><td>OC-6</td><td>Octahedron</td><td>O<sub>h</sub></td></tr>
<tr><td>3</td><td>TPR-6</td><td>Trigonal prism</td><td>D<sub>3h</sub></td></tr>
<tr><td>4</td><td>JPPY-5</td><td>Johnson pentagonal pyramid (J2)</td><td>C<sub>5v</sub></td></tr>

<tr><td rowspan="7">7</td><td>0</td><td>HP-7</td><td>Heptagon</td><td>D<sub>7h</sub></td></tr>
<tr><td>1</td><td>HPY-7</td><td>Hexagonal pyramid</td><td>C<sub>6v</sub></td></tr>
<tr><td>2</td><td>PBPY-7</td><td>Pentagonal bipyramid</td><td>D<sub>5h</sub></td></tr>
<tr><td>3</td><td>COC-7</td><td>Capped octahedron*</td><td>C<sub>3v</sub></td></tr>
<tr><td>4</td><td>CTPR-7</td><td>Capped trigonal prism*</td><td>C<sub>2v</sub></td></tr>
<tr><td>5</td><td>JPBPY-7</td><td>Johnson pentagonal bipyramid (J13)</td><td>D<sub>5h</sub></td></tr>
<tr><td>6</td><td>JETPY-7</td><td>Elongated triangular pyramid (J7)</td><td>C<sub>3v</sub></td></tr>

<tr><td rowspan="13">8</td><td>0</td><td>OP-8</td><td>Octagon</td><td>D<sub>8h</sub></td></tr>
<tr><td>1</td><td>HPY-8</td><td>Heptagonal pyramid</td><td>C<sub>7v</sub></td></tr>
<tr><td>2</td><td>HBPY-8</td><td>Hexagonal bipyramid</td><td>D<sub>6h</sub></td></tr>
<tr><td>3</td><td>CU-8</td><td>Cube</td><td>O<sub>h</sub></td></tr>
<tr><td>4</td><td>SAPR-8</td><td>Square antiprism</td><td>D<sub>4d</sub></td></tr>
<tr><td>5</td><td>TDD-8</td><td>Triangular dodecahedron</td><td>D<sub>2d</sub></td></tr>
<tr><td>6</td><td>JGBF-8</td><td>Johnson-Gyrobifastigium (J26)</td><td>D<sub>2d</sub></td></tr>
<tr><td>7</td><td>JETBPY-8</td><td>Johnson-Elongated triangular bipyramid (J14)</td><td>D<sub>3h</sub></td></tr>
<tr><td>8</td><td>JBTPR-8</td><td>Johnson-Biaugmented trigonal prism (J50)</td><td>C<sub>2v</sub></td></tr>
<tr><td>9</td><td>BTPR-8</td><td>Biaugmented trigonal prism</td><td>C<sub>2v</sub></td></tr>
<tr><td>10</td><td>JSD-8</td><td>Snub disphenoid (J84)</td><td>D<sub>2d</sub></td></tr>
<tr><td>11</td><td>TT-8</td><td>Triakis tetrahedron</td><td>T<sub>d</sub></td></tr>
<tr><td>12</td><td>ETBPY-8</td><td>Elongated trigonal bipyramid (see 8)</td><td>D<sub>3h</sub></td></tr>

<tr><td rowspan="13">9</td><td>0</td><td>EP-9</td><td>Enneagon</td><td>D<sub>9h</sub></td></tr>
<tr><td>1</td><td>OPY-9</td><td>Octagonal pyramid</td><td>C<sub>8v</sub></td></tr>
<tr><td>2</td><td>HBPY-9</td><td>Heptagonal bipyramid</td><td>D<sub>7h</sub></td></tr>
<tr><td>3</td><td>JTC-9</td><td>Triangular cupola (J3) = trivacant cuboctahedron</td><td>C<sub>3v</sub></td></tr>
<tr><td>4</td><td>JCCU-9</td><td>Capped cube (Elongated square pyramid, J8)</td><td>C<sub>4v</sub></td></tr>
<tr><td>5</td><td>CCU-9</td><td>Capped cube</td><td>C<sub>4v</sub></td></tr>
<tr><td>6</td><td>JCSAPR-9</td><td>Capped sq. antiprism (Gyroelongated square pyramid J10)</td><td>C<sub>4v</sub></td></tr>
<tr><td>7</td><td>CSAPR-9</td><td>Capped square antiprism</td><td>C<sub>4v</sub></td></tr>
<tr><td>8</td><td>JTCTPR-9</td><td>Tricapped trigonal prism (J51)</td><td>D<sub>3h</sub></td></tr>
<tr><td>9</td><td>TCTPR-9</td><td>Tricapped trigonal prism</td><td>D<sub>3h</sub></td></tr>
<tr><td>10</td><td>JTDIC-9</td><td>Tridiminished icosahedron (J63)</td><td>C<sub>3v</sub></td></tr>
<tr><td>11</td><td>HH-9</td><td>Hula-hoop</td><td>C<sub>2v</sub></td></tr>
<tr><td>12</td><td>MFF-9</td><td>Muffin</td><td>C<sub>s</sub></td></tr>

<tr><td rowspan="13">10</td><td>0</td><td>DP-10</td><td>Decagon</td><td>D<sub>10h</sub></td></tr>
<tr><td>1</td><td>EPY-10</td><td>Enneagonal pyramid</td><td>C<sub>9v</sub></td></tr>
<tr><td>2</td><td>OBPY-10</td><td>Octagonal bipyramid</td><td>D<sub>8h</sub></td></tr>
<tr><td>3</td><td>PPR-10</td><td>Pentagonal prism</td><td>D<sub>5h</sub></td></tr>
<tr><td>4</td><td>PAPR-10</td><td>Pentagonal antiprism</td><td>D<sub>5d</sub></td></tr>
<tr><td>5</td><td>JBCCU-10</td><td>Bicapped cube (Elongated square bipyramid J15)</td><td>D<sub>4h</sub></td></tr>
<tr><td>6</td><td>JBCSAPR-10</td><td>Bicapped square antiprism (Gyroelongated square bipyramid J17)</td><td>D<sub>4d</sub></td></tr>
<tr><td>7</td><td>JMBIC-10</td><td>Metabidiminished icosahedron (J62)</td><td>C<sub>2v</sub></td></tr>
<tr><td>8</td><td>JATDI-10</td><td>Augmented tridiminished icosahedron (J64)</td><td>C<sub>3v</sub></td></tr>
<tr><td>9</td><td>JSPC-10</td><td>Sphenocorona (J87)</td><td>C<sub>2v</sub></td></tr>
<tr><td>10</td><td>SDD-10</td><td>Staggered dodecahedron (2:6:2)#</td><td>D<sub>2</sub></td></tr>
<tr><td>11</td><td>TD-10</td><td>Tetradecahedron (2:6:2)</td><td>C<sub>2v</sub></td></tr>
<tr><td>12</td><td>HD-10</td><td>Hexadecahedron (2:6:2, or 1:4:4:1)</td><td>D<sub>4h</sub></td></tr>

<tr><td rowspan="7">11</td><td>0</td><td>HP-11</td><td>Hendecagon</td><td>D<sub>11h</sub></td></tr>
<tr><td>1</td><td>DPY-11</td><td>Decagonal pyramid</td><td>C<sub>10v</sub></td></tr>
<tr><td>2</td><td>EBPY-11</td><td>Enneagonal bipyramid</td><td>D<sub>9h</sub></td></tr>
<tr><td>3</td><td>JCPPR-11</td><td>Capped pent. prism (Elongated pentagonal pyramid J9)</td><td>C<sub>5v</sub></td></tr>
<tr><td>4</td><td>JCPAPR-11</td><td>Capped pent. antiprism (Gyroelongated pentagonal pyramid J11)</td><td>C<sub>5v</sub></td></tr>
<tr><td>5</td><td>JAPPR-11</td><td>Augmented pentagonal prism (J52)</td><td>C<sub>2v</sub></td></tr>
<tr><td>6</td><td>JASPC-11</td><td>Augmented sphenocorona (J87)</td><td>C<sub>s</sub></td></tr>

<tr><td rowspan="13">12</td><td>0</td><td>DP-12</td><td>Dodecagon</td><td>D<sub>12h</sub></td></tr>
<tr><td>1</td><td>HPY-12</td><td>Hendecagonal pyramid</td><td>C<sub>11v</sub></td></tr>
<tr><td>2</td><td>DBPY-12</td><td>Decagonal bipyramid</td><td>D<sub>10h</sub></td></tr>
<tr><td>3</td><td>HPR-12</td><td>Hexagonal prism</td><td>D<sub>6h</sub></td></tr>
<tr><td>4</td><td>HAPR-12</td><td>Hexagonal antiprism</td><td>D<sub>6d</sub></td></tr>
<tr><td>5</td><td>TT-12</td><td>Truncated tetrahedron</td><td>T<sub>d</sub></td></tr>
<tr><td>6</td><td>COC-12</td><td>Cuboctahedron</td><td>O<sub>h</sub></td></tr>
<tr><td>7</td><td>ACOC-12</td><td>Anticuboctahedron (Triangular orthobicupola J27)</td><td>D<sub>3h</sub></td></tr>
<tr><td>8</td><td>IC-12</td><td>Icosahedron</td><td>I<sub>h</sub></td></tr>
<tr><td>9</td><td>JSC-12</td><td>Square cupola (J4)</td><td>C<sub>4v</sub></td></tr>
<tr><td>10</td><td>JEPBPY-12</td><td>Elongated pentagonal bipyramid (J16)</td><td>D<sub>6h</sub></td></tr>
<tr><td>11</td><td>JBAPPR-12</td><td>Biaugmented pentagonal prism (J53)</td><td>C<sub>2v</sub></td></tr>
<tr><td>12</td><td>JSPMC-12</td><td>Sphenomegacorona (J88)</td><td>C<sub>s</sub></td></tr>

<tr><td>20</td><td>0</td><td>DD-20</td><td>Dodecahedron†</td><td>I<sub>h</sub></td></tr>
<tr><td rowspan="2">24</td><td>0</td><td>TCU-24</td><td>Truncated cube</td><td>O<sub>h</sub></td></tr>
<tr><td>1</td><td>TOC-24</td><td>Truncated octahedron</td><td>O<sub>h</sub></td></tr>
<tr><td>48</td><td>0</td><td>TCOC-48</td><td>Truncated cuboctahedron</td><td>O<sub>h</sub></td></tr>
<tr><td>60</td><td>0</td><td>TIC-60</td><td>Truncated icosahedron (fullerene)</td><td>I<sub>h</sub></td></tr>
</tbody>
</table>

_Be noted that the index of each shape differs from `SHAPE`'s 2.1 by 1. 
I have a certain suspicion that the original numbering starting from one is due to  Fortran arrays being [silly](https://xkcd.com/163/)._

### User defined polyhedra 

The user can include custom polyhedra into the calculation by pointing to `shape.yaml` files. The shapes specified must have the same number of vertices as the problem shape. 
The syntax for writing these files is as follows:
```yaml
symbol:
  name: 
  symmetry:
  vertices:
    - [ coord_x,  coord_y, coord_z]
    - [ ...
  centre:
    - [ coord_x, coord_y, coord_z]

```
Cosmochlore's YAML parser is somewhat flexible when specifying the key names. It accepts the following synonyms:
| Key | Synonyms | Value type |
| --- | --- | --- |
| `symbol` | None | `string` |
| `name` | None | `string` |
| `symmetry` | `symm` | `string` |
| `vertices` | `ligands` | `[float, float, float]` |
| `centre` | `center`, `metal` | `[float, float, float]` |

**Example .yaml file:**
```yaml
ebcT-6:
  name: Edge-bicapped tetrahedron
  symm: D2d
  vertices:
    - [0.5000,   0.5000,  -0.6495]
    - [1.0000,   1.0000,   0.0000]
    - [1.0000,   0.0000,   1.0000]
    - [0.0000,   1.0000,   1.0000]
    - [0.5000,   0.5000,   1.6495]
    - [0.5000,   0.5000,   0.5000]
  centre:
    - [ 0.0000,   0.0000,   0.0000 ]

fvCU-6:
  name: Face-divacant cube
  symm: C2v
  vertices:
#    - [ 1.000000000000,  0.000000000000, -0.707106780000]
    - [ 0.000000000000,  1.000000000000, -0.707106780000]
#    - [-1.000000000000,  0.000000000000, -0.707106780000]
    - [-0.000000000000, -1.000000000000, -0.707106780000]
    - [ 1.000000000000,  0.000000000000,  0.707106780000]
    - [ 0.000000000000,  1.000000000000,  0.707106780000]
    - [-1.000000000000,  0.000000000000,  0.707106780000]
    - [-0.000000000000, -1.000000000000,  0.707106780000]
  center:
    - [ 0.000000000000,  0.000000000000,  0.000000000000]
```
As per the example, a single .yaml file can contain multiple reference polyhedra. 
Comments are also supported by the YAML parser, every line that starts with `#` will be ignored.

## The `cshm` algorithm

The Continuous Shape Measure of a problem shape $$Q$$ relative to a reference shape $$P$$ is defined as:

```math
S_P(Q)= \mathrm{min} \left( \frac{\sum_i^n |q_i - p_i|^2}{ \sum_i^n |p_i - p_0|^2} \right) × 100
```
where $$N$$ is the number of vertices in the structures we are comparing, $$q_i$$ and $$p_i$$ are the position vectors of the vertices of $$Q$$ and $$P$$, respectively, and $$p_0$$ the geometric centre of the problem structure $$Q$$. 
This minimization is carried out over all rotations, translations and scalings **and** over all $$N!$$ possible permutations of point pairs assignments.

The rotation/translation/scaling part is solved via an SVD-based Kabsch-style alignment. The combinatorial problem: finding the best point pair matches between $$Q$$ and $$P$$ is done by:

- **Automorphism-aware deduplication**: the reference shape's own symmetry group is precomputed, so permutations that are guaranteed to produce identical scores (automorphisms of the reference shape's point group) are never evaluated twice.
- **Branch-and-bound pruning**: partial assignments are bounded using the subadditivity property of the singular-value sum, allowing branches that provably cannot beat the current best score to be discarded early.

## Acknowledgements

This project reimplements the shape-measure methodology originally developed for the `SHAPE` program and continued in [`cosymlib`](https://github.com/GrupEstructuraElectronicaSimetria/cosymlib) by the Electronic Structure Group at the Universitat de Barcelona. Cosmochlore is an independent, from-scratch Rust implementation and is not affiliated with the original authors.

## License

**The code, binaries and sample tests are provided as is, with no warranty of any kind.** This program is licensed under the GNU General Public Licence v3.0 (GPL-3.0). See the `LICENSE.md` file or https://www.gnu.org/licenses/gpl-3.0.html for full terms.

Parts of the `cshm` module are based on `cosymlib` (`shp.f90`), Copyright (c) 2021 Pere Alemany, Efrem Bernuz, Abel Carreras and Miquel Llunell, licensed under the MIT Licence.

Parts of the `odis` module is are based on the [`OctaDist`](https://octadist.github.io/)<sup>3</sup> algorithm, Copyright (c) 2019-2026  Rangsiman Ketkaew et al., licensed under the GNU General Public Licence v3.0 (GPL-3.0).

The program relies on the [`clap`](https://crates.io/crates/clap) crate for argument parsing. Dual-licensed under Apache 2.0 or MIT licences.

The program relies on the [`nalgebra`](https://docs.rs/nalgebra/latest/nalgebra/) crate for the fast linear algebra computations. Licensed under the Apache 2.0 licence.


## References

1. Cirera, J., Ruiz, E., & Alvarez, S. (2005). Continuous Shape Measures as a Stereochemical Tool in Organometallic Chemistry. _Organometallics_, 24(7), 1556–1562. https://doi.org/10.1021/om049150z
2. Alvarez, S., Alemany, P., Casanova, D., Cirera, J., Llunell, M., & Avnir, D. (2005). Shape maps and polyhedral interconversion paths in transition metal chemistry. _Coordination Chemistry Reviews_, 249(17–18), 1693–1708. https://doi.org/10.1016/j.ccr.2005.03.031
3. Ketkaew, R., Tantirungrotechai, Y., Harding, P., Chastanet, G., Guionneau, P., Marchivie, M., & Harding, D. J. (2021). OctaDist: a tool for calculating distortion parameters in spin crossover and coordination complexes. _Dalton Transactions_, 50(3), 1086–1096. https://doi.org/10.1039/d0dt03988h
