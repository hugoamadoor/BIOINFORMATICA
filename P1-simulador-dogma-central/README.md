# Práctica 1: Simulación del dogma central de la biología molecular

Universidad de Las Palmas de Gran Canaria · Grado en Ciencia e Ingeniería de Datos  
Bioinformática (prácticas de laboratorio) · Profesora: María Dolores Afonso Suárez · Curso 2026/2027  
**Grupo 19:** Hugo Amador Hidalgo y Sara Lillo

Programa en Python que simula, paso a paso y en la terminal, el flujo completo de la
información genética a partir de una molécula de ADN:

**ADN → ADN** (replicación) · **ADN → ARN** (transcripción) · **ARN → proteína** (traducción)

En cada etapa se dibujan las hebras alineadas base a base y se explica qué molécula o
enzima actúa en ese momento.

## Requisitos e instalación

- Python 3.8 o superior. **El simulador no necesita ninguna librería externa.**
- Solo para las pruebas automáticas: `pip install -r requirements.txt` (pytest y Biopython).

## Uso

```bash
python simulador.py                     # menú interactivo: elige la molécula y el ritmo
```

También se puede lanzar directamente con opciones:

```bash
python simulador.py --ejemplo --modo paso                              # gen de ejemplo, pulsando Enter
python simulador.py --fasta datos/insulina_INS_cds.fasta --modo auto   # avanza solo
python simulador.py --secuencia ATGGCTAAATAA --modo resumen            # solo resultados
```

| Opción | Qué hace |
|---|---|
| `--ejemplo` / `--secuencia SEC` / `--fasta RUTA` | Molécula de partida (hebra codificante, 5'→3') |
| `--modo paso \| auto \| resumen` | Paso a paso con Enter, automático o solo resultados |
| `--okazaki N` | Tamaño de los fragmentos de Okazaki (por defecto 20 nt) |
| `--cebador N` | Tamaño de los cebadores de ARN (por defecto 3 nt) |
| `--retardo S` | Segundos entre pasos en el modo automático (por defecto 1) |
| `--ancho N` | Bases por línea en los dibujos (por defecto 60; bájalo si la terminal es estrecha) |
| `--sin-color` | Desactiva los colores |

Pruebas automáticas:

```bash
python -m pytest -q tests
```

## La molécula de ejemplo

Es un gen bacteriano pequeño (111 pb) diseñado para que cada etapa tenga sentido biológico:

| Posición | Parte | Para qué sirve |
|---|---|---|
| 1-6 | Caja −35 `TTGACA` | La reconoce el factor sigma de la ARN polimerasa |
| 24-29 | Caja −10 `TATAAT` | Promotor; la transcripción empieza 6 nt después |
| 36 | +1 | Primer nucleótido del ARNm |
| 36-50 | 5'UTR con Shine-Dalgarno `AGGAGG` | Sitio de unión del ribosoma |
| 51-86 | `ATG ... TAA` | Región codificante: 11 aminoácidos (MSKGEELFTGV) |
| 89-111 | Horquilla + `TTTTTTT` | Terminador intrínseco de la transcripción |

También se puede usar cualquier secuencia propia, por ejemplo la región codificante de la
insulina humana (`datos/insulina_INS_cds.fasta`). Si no tiene promotor o terminador, el
programa lo avisa y transcribe la molécula entera.

## Qué muestra cada etapa

### 1. Replicación (ADN → ADN)

Una horquilla sale del origen (extremo izquierdo) y avanza hacia la derecha por bloques.
En cada bloque se ven dos pasos:

- **a)** la **helicasa** abre la doble hélice (`^` marca su posición), las **SSB** protegen las
  hebras sencillas (`~`) y la **primasa** coloca los cebadores de ARN (en minúscula, con `u`);
- **b)** la **ADN polimerasa III** alarga los cebadores en sentido 5'→3'.

```
>> Paso 3a: la helicasa abre las posiciones 41-60 y la primasa coloca los cebadores
   Sup. parental    5' TTGACAGCTAGCTCAGTCCTAGGTATAATGCATCGACAAGGAGGTAAACC
   Rezagada (nueva) 3' AACTGTCGATCGAGTCAggaTCCATATTACGTAGCTGuuc~~~~~~~~~~
   Líder (nueva)    5' uugACAGCTAGCTCAGTCCTAGGTATAATGCATCGACAAG~~~~~~~~~~
   Inf. parental    3' AACTGTCGATCGAGTCAGGATCCATATTACGTAGCTGTTCCTCCATTTGG
```

- La **cadena líder** (verde) crece en el mismo sentido que la horquilla: es continua y
  solo necesita un cebador, en el origen.
- La **cadena rezagada** crece en sentido contrario, así que se fabrica a trozos: los
  **fragmentos de Okazaki** (en colores alternos), cada uno con su propio cebador.
- Al final, la **ADN polimerasa I** sustituye el ARN de los cebadores por ADN y la
  **ADN ligasa** sella las mellas entre fragmentos:

```
   Antes   3' AACTGTCGATCGAGTCA|GGATCCATATTACGTAGCTG|TTCCTCCATTTGG... 5'   (| = mella)
   Después 3' AACTGTCGATCGAGTCAGGATCCATATTACGTAGCTGTTCCTCCATTTGG... 5'   (cadena continua)
```

El programa comprueba que las dos moléculas hijas tienen todas las bases bien emparejadas
(A-T, G-C), que son idénticas a la original y que cada una conserva una hebra parental
(replicación semiconservativa).

### 2. Transcripción (ADN → ARN)

El factor sigma reconoce las cajas −35 y −10. La **ARN polimerasa** lee la **hebra molde**
(la inferior) de 3' a 5' y fabrica el ARNm de 5' a 3', colocando frente a cada base el
ribonucleótido complementario (A→U, T→A, G→C, C→G). Termina al llegar al terminador
(horquilla + poli-U), que suelta el ARN.

```
   Codificante      5' TTGACAGCTAGCTCAGTCCTAGGTATAATGCATCGACAAGGAGGTAAACC
   Molde (leída)    3' AACTGTCGATCGAGTCAGGATCCATATTACGTAGCTGTTCCTCCATTTGG
   ARNm             5'                                    ACAAGGAGGUAAACC
                       -35                    -10         +1
```

Se comprueba que el ARNm es igual a la hebra codificante cambiando T por U.

### 3. Traducción (ARN → proteína)

El ribosoma se une a la secuencia **Shine-Dalgarno** y empieza en el primer **AUG**
(ARNt iniciador con formil-metionina en el sitio P). En cada codón se ven los sitios
**E, P y A** del ribosoma, el **ARNt** que entra con su **anticodón**, el enlace peptídico y la
translocación. Al llegar a un **codón de paro** entra un **factor de liberación** y se suelta
la proteína.

```
>> Elongación: codón 4 = GGA -> Gly
   ARNm 5'- AUG [AGU][AAA][GGA] GAA  GAA  CUU ...-3'
                  E    P    A
   1. Entra en el sitio A el ARNt-Gly con anticodón 3'-CCU-5' (empareja con 5'-GGA-3').
   2. La peptidil transferasa une la cadena al nuevo aminoácido.
   3. El ribosoma avanza un codón; el ARNt vacío sale por el sitio E.
   Cadena: Met-Ser-Lys-Gly
```

Al final, un resumen muestra el flujo completo: ADN original → 2 moléculas hijas → ARNm → proteína.

## Casos de prueba incluidos

En `datos/casos/` hay variantes del gen de ejemplo para comprobar cada funcionalidad:

| Fichero | Qué cambia | Qué debe ocurrir |
|---|---|---|
| `mutacion_sin_sentido.fasta` | AAA→TAA en el codón 3 | Paro prematuro: Met-Ser (2 aa) |
| `mutacion_cambio_sentido.fasta` | GGA→GTA en el codón 4 | Gly→Val: MSKVEELFTGV |
| `mutacion_insercion.fasta` | Una A insertada en el codón 3 | Se desplaza la pauta: no encuentra paro y lee la 3'UTR (20 aa) |
| `inicio_GTG.fasta` | ATG→GTG en el inicio | No hay AUG: no se traduce (el simulador no usa GUG como inicio) |
| `sin_terminador.fasta` | Poli-T del terminador mutado | La transcripción llega al final de la molécula |
| `promotor_propio.fasta` | Gen mínimo de 35 pb | Promotor sin caja −35; proteína Met-Lys-Pro |

Ejemplo: `python simulador.py --fasta datos/casos/mutacion_sin_sentido.fasta --modo resumen`

## Moléculas y enzimas representadas

| Etapa | Molécula o enzima | Función en el simulador |
|---|---|---|
| Replicación | Helicasa | Abre la doble hélice bloque a bloque (`^`) |
| | Proteínas SSB | Protegen las hebras sencillas (`~`) |
| | Topoisomerasa (girasa) | Alivia el superenrollamiento por delante de la horquilla |
| | Primasa | Coloca cebadores de ARN (minúsculas, con `u`) |
| | ADN polimerasa III | Sintetiza la cadena líder (continua) y los fragmentos de Okazaki |
| | ADN polimerasa I | Sustituye los cebadores de ARN por ADN |
| | ADN ligasa | Sella las mellas entre fragmentos de Okazaki |
| Transcripción | ARN polimerasa + factor sigma | Reconoce el promotor, lee el molde 3'→5' y sintetiza el ARNm 5'→3' |
| | Terminador intrínseco | Horquilla + poli-U que suelta el ARN |
| Traducción | Ribosoma (sitios E, P, A) | Lee el ARNm codón a codón |
| | ARNt con anticodón | Lleva el aminoácido que corresponde a cada codón |
| | Peptidil transferasa | Forma el enlace peptídico |
| | Factor de liberación (RF1/RF2) | Reconoce el codón de paro y suelta la proteína |

## Simplificaciones

- Modelo bacteriano (*E. coli*): una sola horquilla que sale de un extremo. En la realidad el
  origen es bidireccional y los cromosomas bacterianos son circulares.
- Escala reducida: fragmentos de Okazaki de 20 nt y cebadores de 3 nt (en *E. coli* son de
  1000-2000 nt y unos 10 nt). Ambos tamaños se pueden cambiar con las opciones.
- Los cebadores de los extremos se rellenan igualmente; en una célula los cubriría la otra
  horquilla (o la telomerasa, en cromosomas lineales eucariotas).
- La transcripción empieza a una distancia fija de la caja −10 (6 nt) y termina en el primer
  tramo `TTTTTT`. No se modela el procesamiento eucariota (caperuza, cola poli-A, splicing).
- El anticodón se calcula por complementariedad exacta, sin tener en cuenta el balanceo
  (*wobble*) de la tercera posición.

## Estructura del código

```
simulador.py              programa principal: menú, opciones y orden de las tres etapas
dogma/
  secuencias.py           reglas de complementariedad, lectura de FASTA y gen de ejemplo
  codigo_genetico.py      tabla de los 64 codones y factores de liberación
  replicacion.py          horquilla, cebadores, cadena líder, rezagada y Okazaki
  transcripcion.py        promotor, ARN polimerasa, molde y terminador
  traduccion.py           ribosoma, ARNt, anticodones y codones de paro
  pantalla.py             colores, pausas y dibujo de hebras alineadas
datos/                    gen de ejemplo e insulina humana en FASTA
tests/test_simulador.py   pruebas (comparan con Biopython el código genético y la traducción)
```

Las tres funciones principales son `simular_replicacion`, `simular_transcripcion` y
`simular_traduccion`. Cada una dibuja su etapa y devuelve el resultado, que pasa a la siguiente.
