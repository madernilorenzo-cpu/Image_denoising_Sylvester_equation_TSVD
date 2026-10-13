# Decomposizione ai valori singolari di operatori di Sylvester con applicazione all'imaging
**Singular value decomposition of Sylvester operators with application to imaging**

[Italiano](#italiano) | [English](#english)

---

## Italiano

Codice MATLAB che accompagna la tesi triennale in Analisi Numerica:

> **Decomposizione ai valori singolari di operatori di Sylvester con
> applicazione all'imaging** — Lorenzo Maderni, Corso di Laurea in
> Matematica, Alma Mater Studiorum - Università di Bologna, anno
> accademico 2022/2023. Relatore: Davide Palitta.

### Il problema

Data un'immagine digitale `F` (matrice `N x M`) degradata da rumore bianco
gaussiano, `G = F + R`, si vuole ricostruire `F` partendo da `G`.
Il problema di denoising si formalizza (Sezione 2.2 della tesi) come il
sistema lineare

```
g = H f,       H = I (x) A + B' (x) I
```

dove `f = vec(F)`, `A` e `B` derivano dalle matrici di covarianza
separabili dell'immagine (prodotto di Kronecker, eq. (2.8)). `H` è una
matrice `NM x NM` di grandi dimensioni e mal condizionata.

Due metodi a confronto (Capitolo 3):

1. **TSVD di H**: si sostituisce `H` con la sua migliore approssimazione
   di rango `k` (Teorema di Eckart-Young), costo proibitivo `O(n^6)`.
2. **Equazione di Sylvester**: si risolve

   ```
   A X + X B = C
   ```

   con l'algoritmo di Bartels-Stewart (funzione `lyap` di MATLAB,
   Sezione 2.1), costo `O(n^3)`, oppure con le versioni troncate
   `A_k1`, `B_k2` (TSVD di A e B, Sezione 3.2):

   ```
   A_k1 X + X B_k2 = C
   ```

   Troncando una sola delle due matrici l'equazione resta a soluzione
   unica: negli esperimenti della tesi si troncano gli ultimi 103
   valori singolari di `A` su 350 e nessuno di `B`.

### Struttura

```
├── src/                        funzioni numeriche
│   ├── build_sylvester_data.m  costruzione di A, B, C (Sez. 3.3.1)
│   ├── bartels_stewart.m       algoritmo di Bartels-Stewart (Sez. 2.1)
│   ├── sylvester_solve.m       wrapper: lyap se disponibile, altrimenti Bartels-Stewart
│   ├── sylvester_tsvd_solve.m  eq. di Sylvester con TSVD di A e B (Sez. 3.2)
│   ├── kron_operator.m         H = I (x) A + B' (x) I
│   ├── covariance_matrix.m     matrice di Toeplitz R_z, eq. (3.15)
│   ├── psf_matrix.m            matrici della PSF (Sez. 3.3.1)
│   └── load_example_image.m    immagine di test (con fallback sintetico)
├── experiments/
│   ├── experiment_denoising.m  ricostruzione di immagini (Sez. 3.3.3)
│   └── experiment_approxH.m    SigmaH vs ThetaA+B (Sez. 3.3.2, Fig. 3.1)
├── tests/
│   ├── run_all_tests.m
│   ├── test_build_sylvester_data.m
│   ├── test_bartels_stewart.m
│   └── test_sylvester_tsvd.m
└── output/                     figure generate (non versionata)
```

### Requisiti

- MATLAB (verificato su R2022a). Non servono toolbox: il codice usa solo
  funzioni del MATLAB base. La funzione `lyap` (Control System Toolbox),
  usata negli esperimenti originali, viene impiegata se installata,
  altrimenti si ricade sull'implementazione autonoma di Bartels-Stewart.

### Uso

```matlab
% test (dalla cartella tests/, oppure: matlab -batch "run('...\run_all_tests.m')")
run('tests\run_all_tests.m')

% esperimento principale: denoising 350 x 350 (Sez. 3.3.3)
run('experiments\experiment_denoising.m')   % figure in output\

% confronto valori singolari di H con le somme sigma_i(A) + sigma_j(B)
run('experiments\experiment_approxH.m')
```

Parametri usati negli esperimenti (valori della tesi): `rho_x = 0.4310`,
`rho_y = 0.9331`, varianza del rumore `sigma2 = 1.825e-3`, `P = I`,
`k1 = 350 - 103`, `k2 = 350`.

### Risultati

Con immagine di test 350 x 350 e rumore gaussiano a varianza 1.825e-3
(errori relativi in norma di Frobenius rispetto all'immagine originale):

| Ricostruzione | Errore relativo |
|---|---|
| Immagine degradata | 1.13e-1 |
| Equazione di Sylvester completa | 8.1e-2 |
| Equazione di Sylvester con TSVD (k1 = 247) | 7.5e-2 |

La ricostruzione troncata è visivamente indistinguibile da quella
completa, come osservato nella tesi (Figura 3.4); l'omogeneizzazione dei
colori resta al prezzo di una minore nitidezza dei bordi. Nel confronto
`SigmaH` contro `ThetaA+B` vale la disuguaglianza `SigmaH(k) <= ThetaA+B(k)`
(eq. (3.16)) e le due curve sono quasi sovrapposte.

### Bibliografia essenziale

1. R. H. Bartels, G. W. Stewart, *Solution of the Matrix Equation
   AX + XB = C*, Comm. ACM 15 (1972), 820-826.
2. D. Calvetti, L. Reichel, *Application of ADI Iterative Methods to the
   Restoration of Noisy Images*, SIAM J. Matrix Anal. Appl. 17 (1996),
   165-174.
3. C. F. Van Loan, *The ubiquitous Kronecker product*, J. Comput. Appl.
   Math. 123 (2000), 85-100.

---

## English

MATLAB code accompanying the undergraduate (bachelor) thesis in
Numerical Analysis:

> **Decomposizione ai valori singolari di operatori di Sylvester con
> applicazione all'imaging** (*Singular value decomposition of Sylvester
> operators with application to imaging*) — Lorenzo Maderni, Bachelor's
> Degree in Mathematics, Alma Mater Studiorum - University of Bologna,
> academic year 2022/2023. Advisor: Davide Palitta.

### The problem

Given a digital image `F` (an `N x M` matrix) degraded by white Gaussian
noise, `G = F + R`, we want to reconstruct `F` starting from `G`.
The denoising problem is formalized (Section 2.2 of the thesis) as the
linear system

```
g = H f,       H = I (x) A + B' (x) I
```

where `f = vec(F)`, and `A` and `B` come from the separable covariance
matrices of the image (Kronecker product, eq. (2.8)). `H` is a large,
ill-conditioned `NM x NM` matrix.

Two methods are compared (Chapter 3):

1. **TSVD of H**: `H` is replaced by its best rank-`k` approximation
   (Eckart-Young theorem), with a prohibitive `O(n^6)` cost.
2. **Sylvester equation**: one solves

   ```
   A X + X B = C
   ```

   with the Bartels-Stewart algorithm (MATLAB's `lyap` function,
   Section 2.1), with an `O(n^3)` cost, or with the truncated versions
   `A_k1`, `B_k2` (TSVD of A and B, Section 3.2):

   ```
   A_k1 X + X B_k2 = C
   ```

   Truncating only one of the two matrices keeps the equation uniquely
   solvable: in the experiments of the thesis the last 103 singular
   values of `A` are truncated out of 350, and none of `B`.

### Structure

```
├── src/                        numerical functions
│   ├── build_sylvester_data.m  construction of A, B, C (Sec. 3.3.1)
│   ├── bartels_stewart.m       Bartels-Stewart algorithm (Sec. 2.1)
│   ├── sylvester_solve.m       wrapper: lyap if available, else Bartels-Stewart
│   ├── sylvester_tsvd_solve.m  Sylvester equation with TSVD of A and B (Sec. 3.2)
│   ├── kron_operator.m         H = I (x) A + B' (x) I
│   ├── covariance_matrix.m     Toeplitz matrix R_z, eq. (3.15)
│   ├── psf_matrix.m            PSF matrices (Sec. 3.3.1)
│   └── load_example_image.m    test image (with synthetic fallback)
├── experiments/
│   ├── experiment_denoising.m  image reconstruction (Sec. 3.3.3)
│   └── experiment_approxH.m    SigmaH vs ThetaA+B (Sec. 3.3.2, Fig. 3.1)
├── tests/
│   ├── run_all_tests.m
│   ├── test_build_sylvester_data.m
│   ├── test_bartels_stewart.m
│   └── test_sylvester_tsvd.m
└── output/                     generated figures (not tracked)
```

### Requirements

- MATLAB (verified on R2022a). No toolboxes are needed: the code only
  uses core MATLAB functions. The `lyap` function (Control System
  Toolbox), used in the original experiments, is employed if installed;
  otherwise the code falls back to the standalone Bartels-Stewart
  implementation.

### Usage

```matlab
% tests (from the tests/ folder, or: matlab -batch "run('...\run_all_tests.m')")
run('tests\run_all_tests.m')

% main experiment: 350 x 350 denoising (Sec. 3.3.3)
run('experiments\experiment_denoising.m')   % figures in output\

% comparison of the singular values of H with the sums sigma_i(A) + sigma_j(B)
run('experiments\experiment_approxH.m')
```

Parameters used in the experiments (values from the thesis):
`rho_x = 0.4310`, `rho_y = 0.9331`, noise variance `sigma2 = 1.825e-3`,
`P = I`, `k1 = 350 - 103`, `k2 = 350`.

### Results

With a 350 x 350 test image and Gaussian noise with variance 1.825e-3
(relative errors in the Frobenius norm against the original image):

| Reconstruction | Relative error |
|---|---|
| Degraded image | 1.13e-1 |
| Full Sylvester equation | 8.1e-2 |
| Sylvester equation with TSVD (k1 = 247) | 7.5e-2 |

The truncated reconstruction is visually indistinguishable from the
full one, as observed in the thesis (Figure 3.4); the color
homogenization still comes at the cost of less sharp edges. In the
`SigmaH` vs `ThetaA+B` comparison the inequality
`SigmaH(k) <= ThetaA+B(k)` (eq. (3.16)) holds and the two curves are
almost superimposed.

### Key references

1. R. H. Bartels, G. W. Stewart, *Solution of the Matrix Equation
   AX + XB = C*, Comm. ACM 15 (1972), 820-826.
2. D. Calvetti, L. Reichel, *Application of ADI Iterative Methods to the
   Restoration of Noisy Images*, SIAM J. Matrix Anal. Appl. 17 (1996),
   165-174.
3. C. F. Van Loan, *The ubiquitous Kronecker product*, J. Comput. Appl.
   Math. 123 (2000), 85-100.
