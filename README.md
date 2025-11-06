# README — Integración por métodos numéricos (Newton–Cotes + Adaptativos)

Este repositorio reúne **scripts Octave/MATLAB** (compatibles con ambos) para integrar funciones y datos tabulados usando **Trapecios**, **Simpson 1/3**, **Simpson 3/8**, versiones **automáticas por tolerancia**, **a priori**, y **adaptativas**.  
Todos los `.m` están en **ASCII**, **una función por archivo**, sin `...` ni operadores ternarios.


---

## Tabla rápida: qué usar y cuándo

| Archivo | Propósito | Cuándo usar | Notas clave |
|---|---|---|---|
| `trapecios_xy.m` | Trapecios compuesto con datos `x,y` | Datos medidos equiespaciados, integración rápida | Recorta a `[a,b]`, valida paso casi uniforme |
| `trapecios_xy_auditoria.m` | Trapecios con auditoría | Necesitás desglose intervalo a intervalo | Devuelve estructura con aportes parciales |
| `simpson13_xy_strict.m` | Simpson 1/3 con datos `x,y` | `x,y` equiespaciados y #subintervalos en `[a,b]` **par** | Si no se cumplen condiciones → **error** |
| `simpson38_xy_strict.m` | Simpson 3/8 con datos `x,y` | `x,y` equiespaciados y #subintervalos **múltiplo de 3** | Si no se cumplen condiciones → **error** |
| `simpson13_f.m` | Simpson 1/3 con `f(x)` y `n` fijo | Querés controlar `n` (benchmark/teoría) | Devuelve `[I,h]`, requiere `n` **par** |
| `simpson13_adaptativo_f.m` | Simpson 1/3 adaptativo | Dificultad **localizada**; ahorro de evaluaciones | Test local `1/15`, malla **no uniforme** |
| `simpson38_adaptativo_f.m` | Simpson 3/8 adaptativo | Igual que arriba pero con 3/8 | Reutiliza evaluaciones; malla no uniforme |
| `simpson13_f_auto.m` | Simpson 1/3 uniforme auto (halving) | Tenés **tolerancia** y querés que el código ajuste `n` | Error a posteriori `|I_{2n}-I_n|/15` |
| `simpson13_f_auto_min.m` | 1/3 auto buscando **menor `n` par** | TP pide “primer `h` (mínimo `n`) que cumple tol” | Búsqueda fina entre `n` y `2n` |
| `simpson38_f_auto.m` | Simpson 3/8 uniforme auto | Igual que 1/3 auto pero para 3/8 | Requiere `n` **múltiplo de 3** |
| `simpson38_f_auto_min.m` | 3/8 auto con **menor `n` múltiplo de 3** | Idem anterior, versión 3/8 | Búsqueda fina en saltos de 3 |
| `simpson13_f_apriori.m` | 1/3 **a priori** con cota de `|f''''|` | Examen/Informe: “describa `h` que garantiza tol” | Usa `M4 ≥ max|f''''|`, ajusta `n` a **par** |
| `cota_M4_numerica.m` | Estima cota de `|f''''|` | `f''''` analítica engorrosa | Diferencias finitas + factor de seguridad |
| `simpson13_f_apriori_hibrido.m` | A priori + **verificación** | Camino seguro: fija `n` por cota y valida | Si no alcanza, **refina** automático |
| `simpson13_selector.m` | **Orquestador** (a priori / halving / adaptativo) | Uso “por defecto” con tol si no sabés qué conviene | Si pasás `opts.M4`, usa a priori; sino halving; fallback adaptativo |
| `benchmark_simpson13_demo.m` | Mini benchmark 1/3 (auto vs adaptativo) | Comparar costo/error en casos suave vs pico | Útil para informes y elección de estrategia |

---

## Reglas rápidas de elección

- **“Describa `h`” (malla uniforme, por teoría):**  
  1) Calculá `M4 ≥ max|f''''|` (o estimá con `cota_M4_numerica`).  
  2) Usá **`simpson13_f_apriori`** (o `simpson13_selector` con `opts.M4`).  
  3) Ajustá `n` a **par** y, si podés, verificá con una pasada a posteriori.

- **Tengo tolerancia y quiero el menor `n` uniforme:**  
  - 1/3: **`simpson13_f_auto_min`**.  
  - 3/8: **`simpson38_f_auto_min`**.

- **Tengo tolerancia pero no `M4` (malla uniforme):**  
  - 1/3: **`simpson13_f_auto`** (o `simpson13_selector` sin `M4`).  
  - Si luego querés minimizar `n`, pasá a la variante `_auto_min`.

- **Integración con datos tabulados `x,y`:**  
  - Si cumplen requisitos de Simpson → `simpson13_xy_strict` o `simpson38_xy_strict`.  
  - Si no, o querés algo robusto → `trapecios_xy` (o su versión con auditoría).

- **Función con comportamiento “raro” localizado (picos, capas límite):**  
  - 1/3: **`simpson13_adaptativo_f`**.  
  - 3/8: **`simpson38_adaptativo_f`**.  
  - Recordá: no hay un único `h`; se entrega **malla no uniforme** (`info.accepted`).

---

## Fórmulas de error (para “describir `h`”)

- **Simpson 1/3 (malla uniforme):**

$$
|E| \le \frac{(b-a)}{180}\,h^{4}\,M_4,\quad M_4 \ge \max_{[a,b]}|f^{(4)}(x)|.
$$

**Despeje:**

$$
h \le \left(\frac{180\,\text{tol}}{(b-a)\,M_4}\right)^{1/4},\quad
n=\left\lceil\frac{b-a}{h}\right\rceil\ \text{(ajustar a par)}.
$$

- **Simpson 3/8 (malla uniforme):**

$$
|E| \le \frac{(b-a)}{80}\,h^{4}\,M_4
\Rightarrow
h \le \left(\frac{80\,\text{tol}}{(b-a)\,M_4}\right)^{1/4},\quad
n\ \text{múltiplo de 3}.
$$

- **Estimación a posteriori (ambos, orden 4):**

$$
\mathrm{err\_est} \approx \frac{|I_{2n}-I_{n}|}{2^{4}-1}=\frac{|I_{2n}-I_{n}|}{15}.
$$

---

## Ejemplos mínimos

**1) 1/3 uniforme por tolerancia (auto):**
```matlab
f = @(x) sin(x);
[I,h,n,err,fe] = simpson13_f_auto(f, 0, pi, 1e-8);
```

**2) 1/3 con menor n par que cumple:**
```matlab
f = @(x) exp(-x.^2);
[I,h,n,err,fe] = simpson13_f_auto_min(f, 0, 1, 1e-6);
```

**3) 1/3 a priori (con M4) + verificación:**
```matlab
f = @(x) x.^2 .* sin(x);
a = 0; b = pi/4; tol = 1e-3;
M4 = sqrt(2) * (pi + 6 - pi^2/32); % cota analítica en [0, pi/4]
[I,h,n,bound] = simpson13_f_apriori(f, a, b, tol, M4);
```

**4) Selector (decide por vos):**
```matlab
f = @(x) sqrt(x);
[I, out] = simpson13_selector(f, 0, 1, 1e-8); % sin M4: halving; si falla, adaptativo
```

**5) Adaptativo 1/3 (malla no uniforme):**
```matlab
f = @(x) 1 ./ (1 + 100*(x-0.7).^2);
[I, info] = simpson13_adaptativo_f(f, 0, 1, 1e-8);
% info.accepted: filas [a_i, m_i, b_i] de los tramos aceptados
```

**6) Datos tabulados (Simpson 1/3 estricto):**
```matlab
% x equiespaciado; subintervalos en [a,b] deben ser pares
I = simpson13_xy_strict(x, y, a, b);
```

---

## Notas y advertencias

- **Requisitos de malla:** 1/3 → `n` **par**; 3/8 → `n` **múltiplo de 3**.  
- **`f` evaluable y finita** en todos los nodos; singularidades internas invalidan las cotas.  
- **Adaptativo:** controla error **local**, devuelve la **suma total** y la malla usada (`info`).  
- **Cuando `M4` es infinito** (p.ej., `sqrt(x)` en 0): la cota a priori no aplica; usar halving a posteriori o adaptativo y dejarlo aclarado en el informe.  
- **Conteo de evaluaciones:** los scripts reportan evaluaciones aproximadas (útil para comparar eficiencia).

---

## Sugerencia de uso “día a día”

- Caso general con tol y sin pensar mucho → **`simpson13_selector`**.  
- Te piden “describir `h`” → **`simpson13_f_apriori`** (o selector con `opts.M4`).  
- Querés **el menor `n`** uniforme → **`simpson13_f_auto_min`** (o versión 3/8).  
- Datos `x,y` → `simpson13_xy_strict` / `simpson38_xy_strict` (o `trapecios_xy` si no cumplen).
