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
Run=./tema2_blas: N=400: Time=0.058871  OK
Run=./tema2_blas: N=600: Time=0.129748
Run=./tema2_blas: N=800: Time=0.253782  OK
Run=./tema2_blas: N=1000: Time=0.489693
Run=./tema2_blas: N=1200: Time=0.865853 OK
Run=./tema2_blas: N=1400: Time=1.322491
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
Run=./tema2_neopt: N=400: Time=1.093540   OK
Run=./tema2_neopt: N=600: Time=3.547960
Run=./tema2_neopt: N=800: Time=8.367974   OK
Run=./tema2_neopt: N=1000: Time=16.104860
Run=./tema2_neopt: N=1200: Time=28.054720 OK
Run=./tema2_neopt: N=1400: Time=45.481689
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
Run=./tema2_opt_m: N=400: Time=0.161236  OK
Run=./tema2_opt_m: N=600: Time=0.538191
Run=./tema2_opt_m: N=800: Time=1.265819  OK
Run=./tema2_opt_m: N=1000: Time=2.416722
Run=./tema2_opt_m: N=1200: Time=4.049050 OK
Run=./tema2_opt_m: N=1400: Time=6.606329
<<< Bonus=10p >>>
```

## Analiză valgrind

Nu există memory leak-uri. Rezultatele rulărilor se află în `analysis/`.

## Analiză cachegrind

Primul lucru care iese în evidență este numărul de instrucțiuni (5.7mld pentru neopt,
1.2 mld pentru opt_m și doar 200 mil pentru blas). Totuși, numărul de miss-uri în
cache-ul pentru instrucțiuni este aproape inexistent în toate cele 3 implementări.

Următorul fapt observat este accesul la memorie - opt_m accesează de 10 ori mai
puține date din memorie față de varianta neoptimizată. Acest salt uriaș poate fi
datorat optimizării constantelor și folosirii a cât mai multe registre în loc de
variabile alocate pe stivă.

Deși procentual în cazul implementării optimizate sunt mai multe missuri, acestea
sunt în valoare absolută cu mult mai puține față de varianta neoptimizată, fapt
datorat optimizărilor de reordonare a buclelor pentru a avea acces secvențial.

Ultimul lucru observat este numărul de branch-uri mai mic în cazul implementării
optimizate (loop unrolling).

În urma analizării executabilului blas se observă că acesta exploatează foarte bine
accesul la memorie, având atât procentual, cât și în valori absolute mult mai puține
miss-uri, branch-uri și mispredicts față de celelalte 2 implementări.

Rezultatele rulărilor se află în `analysis/`.

## Analiză comparativă

