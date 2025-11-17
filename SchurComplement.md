# Complemento di Schur Sequenziale

L’obiettivo di questa fase era creare e validare un prototipo **sequenziale** del metodo del **Complemento di Schur**, utile per verificare la correttezza algebrica prima di introdurre la parallelizzazione MPI.

---

## 💡 Cosa Abbiamo Implementato

### Nuova Classe: `SchurSequentialSolver`
Non abbiamo modificato le classi core (`LinearSys`, `TridiagMat`).  
Abbiamo introdotto una nuova classe manager (`SchurSequentialSolver`) che **usa** le classi esistenti tramite composizione.

Contiene:

- `std::vector<LinearSys>` per i blocchi locali \( A_{ii} \)
- `std::unique_ptr<LinearSys>` per il sistema di Schur \( S \)
- gestione dei casi limite

### Rispetto dei Contratti Core
`TridiagMat` (e quindi `ThomaSolver`) richiede matrici con \( N \ge 2 \).  
`SchurSequentialSolver` intercetta e gestisce questo vincolo.

---

# 1. Decomposizione del Dominio

Partiamo dal sistema globale:

```math
A \, \mathbf{x} = \mathbf{f}
```

Il dominio 1D è suddiviso in \( P \) blocchi. Le incognite si separano in:

- incognite interne \( \mathbf{u}_i \)  
- incognite di interfaccia \( \mathbf{u}_e \) (dimensione \( P-1 \))

Il sistema diventa:

```math
\begin{bmatrix}
A_{ii} & A_{ie} \\
A_{ei} & A_{ee}
\end{bmatrix}
\begin{bmatrix}
\mathbf{u}_i \\
\mathbf{u}_e
\end{bmatrix}
=
\begin{bmatrix}
\mathbf{f}_i \\
\mathbf{f}_e
\end{bmatrix}
```

### Dettagli Implementativi

Il costruttore calcola:

- `domain_sizes`
- `domain_starts`
- `interface_indices`

Funziona anche per suddivisioni non uniformi (es. 50/4).

---

# 2. Complemento di Schur

Prima si risolvono le incognite di interfaccia:

```math
S \, \mathbf{u}_e = \mathbf{f}^{*}
```

dove:

```math
S = A_{ee} - A_{ei}\, A_{ii}^{-1}\, A_{ie}
```

```math
\mathbf{f}^{*}
= \mathbf{f}_e - A_{ei}\, A_{ii}^{-1}\, \mathbf{f}_i
```

### Note importanti

- \( A_{ii}^{-1} \) **non viene mai calcolata esplicitamente**:  
  si risolve con `ThomaSolver()`.
- \( S \) ha dimensione \( (P-1) \times (P-1) \) e rimane **tridiagonale**.

---

# 🐛 Gestione Casi Limite

### ✔️ Caso \( P = 1 \)
Nessuna interfaccia → niente Schur.  
`SchurSequentialSolver` chiama direttamente il solver locale.

---

### ✔️ Caso \( P = 2 \)
Una sola interfaccia → \( S \) è uno scalare:

```math
s \, u = f
```

`ThomaSolver` non accetta matrici 1×1 → gestione manuale:

```math
u = \frac{f}{s}
```

---

### ✔️ Caso \( P \ge 3 \)
Caso generale.

- \( S \) ha dimensione ≥ 2  
- si costruisce in `PreProcess()`  
- si risolve in `solve()` con `ThomaSolver()`

---

# 3. Algoritmo Sequenziale

## 🔹 Pre-processing (una sola volta)

Per ogni dominio si estraggono:

- \( A_{ii} \)
- \( A_{ie} \), \( A_{ei} \)
- contributi ad \( A_{ee} \)

Per costruire \( S \) si calcola:

```math
A_{ei} \, A_{ii}^{-1} \, A_{ie}
```

risolvendo sistemi locali.

---

## 🔹 Solve (ogni timestep)

### 1. Solve Locale 1

```math
A_{ii}\, \mathbf{z}_p = \mathbf{f}_i
```

### 2. Assemblaggio RHS di Schur

```math
\mathbf{f}^{*} = \mathbf{f}_e - A_{ei}\mathbf{z}
```

### 3. Solve di Interfaccia

```math
S \, \mathbf{u}_e = \mathbf{f}^{*}
```

### 4. Solve Locale 2

```math
\mathbf{g}_p = \mathbf{f}_i - A_{ie}\mathbf{u}_e
```

```math
A_{ii} \, \mathbf{u}_i = \mathbf{g}_p
```

---

# ✅ Risultato Finale

La soluzione globale \( \mathbf{x} \) è ottenuta concatenando:

- tutti i vettori \( \mathbf{u}_i \)
- il vettore di interfaccia \( \mathbf{u}_e \)

## ✔️ Test (GTest)
Il test `SolutionAgreesWithDirectSolve` confronta:

- soluzione del solutore Schur
- solve diretto del sistema globale

Tutti i casi passano:

- \( P = 1 \)  
- \( P = 2 \)  
- \( P \ge 3 \)  
- split non uniformi (50/4, 17/5, …)
