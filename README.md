# Matrix-Multiply - Cosmin-Răzvan VANCEA - 333CA

## Implementare

### Varianta (C)BLAS

Am impărțit calculul în următoarele subcalcule:

* `C = A' x A`
  * Inițial pun în C rezultatul calcului de mai sus.
  * Pentru efectuarea calculului am folosit funcția `DTRMM` în forma:
  `C := alpha*op( A )*A`, cu `alpha = 1.0`, `op(A) = A**T` și `A = UPPER TRIANGULAR`

* `AB = A x B`
  * Pentru efectuarea calculului am folosit funcția `DTRMM` în forma:
  `AB := alpha*op( A )*B`, cu `alpha = 1.0`, `op(A) = A` și `A = UPPER TRIANGULAR`

* `C = [(A x B) x B'] + (A' x A) = AB x B' + C`
  * Pentru efectuarea calculului am folosit funcția `DGEMM` în forma:
  `C := alpha*op( AB )*op( B ) + beta*C`, cu `alpha = 1.0`, `op(AB) = AB`,
  `beta = 1.0`și `op(B) = B**T`

#### Rezultate pe coada ibm.nehalem:

```
Run=./tema2_blas: N=400: Time=0.057834  OK
Run=./tema2_blas: N=800: Time=0.278843  OK
Run=./tema2_blas: N=1200: Time=0.838009 OK
```

### Varianta neoptimizată

Am impărțit calculul în următoarele subcalcule:

* `AtA = A' x A`
  * Cum A este superior triunghiulară, rezultă că A' este inferior triunghiulară
  * Astfel, pot să reduc numărul de înmulțiri ignorând zerourile din ambele matrici
  * Pentru a face asta, k ia valori în intervalul: `k = [0 .. min(i, j)]`

* `BBt = B x B'`
  * Nu se știe nimic despre B, așadar se folosește algoritmul clasic de înmulțire

* `ABBt = A x (B x B') = A x BBt`
  * Se ignoră zerourile din partea inferioară a lui A (`k = [i .. N]`)

* `C = [A x (B x B')] + (A' x A) = ABBt + AtA`
  * Se însumează rezultatele temporare și se salvează în totul în C

#### Rezultate pe coada ibm.nehalem:

```
Run=./tema2_neopt: N=400: Time=1.091684   OK
Run=./tema2_neopt: N=800: Time=8.359577   OK
Run=./tema2_neopt: N=1200: Time=27.939932 OK
```

### Varianta optimizată

Am impărțit calculul în următoarele subcalcule:

* `At = A'`
  * Precalculez A' deoarece mai departe va trebui să fac calculul `A' x A`,
  iar după cum s-a observat și pe implementarea neoptimizată, accesele la matrici
  în varianta aceasta sunt nesecvențiale (`AtA[i][j] += A[k][i] * A[k][j]`).
  * Cum A este superior triunghiulară, rezultă că A' este inferior triunghiulară.
  * Formulă de calcul a transpusei este: `At[j][i] = A[i][j]`
  * Optimizări:
    1. Se sare peste zerourile de sub diagonala matricii A (`j = [0 .. i]`)
    2. Adresarea elementelor se face prin pointer arithmetic
    3. Din formula de calcul se observă că accesul la memorie este secvențial
    pentru A, dar nesecvențial pentru At (unlucky)

* `C = A' x A = A' x (A')' = At * At'`
  * Noua formulă de calcul este: `C[i][j] += At[i][k] * At[j][k]`
  * Optimizări:
    1. Se sare peste zerourile din ambele matrici
    2. Adresarea elementelor se face prin pointer arithmetic
    3. Din formula de calcul se observă că accesul la memorie este optim
    (C - constant, At - secvențial, At - secvențial)
    4. Se folosesc regiștri pentru sume parțiale
    5. Se înlocuiește apelul de funcție `min` cu codul inline.

* `BBt = B x B'`
  * Formula de calcul este: `BBt[i][j] += B[i][k] * B[j][k]`
  * Calcul este similar celui anterior, singura diferență fiind că B nu este matrice
  superior triunghiulară. Din această cauză nu se ignoră înmulțiri.
  * Optimizări:
    * (-) Ignorarea zerourilor nu se mai poate face
    * (+) În schimb, se garantează că N este multiplu de 40, deci implicit 8 și
    se face loop unrolling la cea mai din interior buclă

* `ABBt = A x (B x B') = A x BBt`
  * Formula de calcul este: `ABBt[i][j] += A[i][k] * BBt[k][j]`
  * Optimizări:
    1. Se sare peste zerourile din matricea A
    2. Adresarea elementelor se face prin pointer arithmetic
    3. Se observă din formulă că accesul la BBt este nesecvențial; se impune
    o reordonare a buclelor astfel: i - k - j
    4. După reordonare, accesul la memorie este optim (ABBt - secvențial,
    A - constant, BBt - secvențial)
    4. Se folosesc regiștri pentru constante
    5. Se face loop unrolling la cea mai din interior buclă

* `C = [A x (B x B')] + (A' x A) = ABBt + C`
  * Se adună la C(= A' x A) rezultatul calculului anterior.
  * Formula de calcul este: `C[i][j] += ABBt[i][j]`
  * Optimizări:
    1. Adresarea elementelor se face prin pointer arithmetic
    2. Accesul la memorie este optim (secvențial)
    3. Se face loop unrolling


#### Note generale:

Folosesc regiștri pentru aproape toate variabilele, singurele variabile alocate
pe stack fiind cele accesate rar (base pointerul pentru matricele temporare).
Cu ajutorul unui decompiler am confirmat că se respectă keyword-ul `register`,
singurele variabile aruncate pe stivă fiind următoarele:

```c
double *A; // [rsp+10h] [rbp-60h]
double *ABBt; // [rsp+20h] [rbp-50h]
double *BBt; // [rsp+28h] [rbp-48h]
double *At; // [rsp+30h] [rbp-40h]
double *C; // [rsp+38h] [rbp-38h]
```

#### Rezultate pe coada ibm.nehalem:

(tbf, testul N=1200 e hit or miss; alternează între ~4.05s și ~4.30s)

```
Run=./tema2_opt_m: N=400: Time=0.161066  OK
Run=./tema2_opt_m: N=800: Time=1.270034  OK
Run=./tema2_opt_m: N=1200: Time=4.057611 OK
<<< Bonus=10p >>>
```
