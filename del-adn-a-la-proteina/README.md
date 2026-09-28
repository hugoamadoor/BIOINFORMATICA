# Del ADN a la proteína: replicación, transcripción, traducción y splicing con Biopython

Universidad de Las Palmas de Gran Canaria · Grado en Ciencia e Ingeniería de Datos
Bioinformática · Profesora: María Dolores Afonso Suárez · Curso 2026/2027
**Grupo 19:** Hugo Amador Hidalgo y Sara Lillo

El informe en PDF está en [`informe/informe.pdf`](informe/informe.pdf). Este README tiene el mismo contenido, además de las instrucciones de uso del código.

## Estructura del repositorio

```
src/
  dogma.py               funciones comunes (complementaria, replicación, transcripción, ORF, traducción)
  ej1_replicacion.py     ejercicio 1
  ej2_transcripcion.py   ejercicio 2 (lee FASTA)
  ej3_traduccion.py      ejercicio 3
  ej4_splicing.py        ejercicio 4 (gen de juguete + consulta a Ensembl)
  ej5_proteinas.py       ejercicio 5 (propiedades del péptido + resumen de un fichero PDB)
  ej6_pipeline.py        ejercicio 6 (pipeline completo con registro de cada paso)
data/
  ej2_codificante.fasta, ej2_molde.fasta
  INS_ENST00000381330_cds.fasta   CDS de la insulina humana (Ensembl, GRCh38)
  P01308_INS_HUMAN.fasta          proteína de referencia (UniProt)
tests/test_dogma.py      comprueba que el código reproduce los resultados obtenidos a mano
informe/                 informe en LaTeX y PDF
```

## Instalación y uso

```bash
pip install -r requirements.txt
python src/ej1_replicacion.py
python src/ej2_transcripcion.py data/ej2_codificante.fasta
python src/ej2_transcripcion.py data/ej2_molde.fasta --hebra molde
python src/ej3_traduccion.py
python src/ej4_splicing.py
python src/ej4_splicing.py --ensembl FGFR2        # requiere conexión a Internet
python src/ej5_proteinas.py
python src/ej5_proteinas.py --pdb 4HHB.pdb        # fichero descargado de https://www.rcsb.org
python src/ej6_pipeline.py data/INS_ENST00000381330_cds.fasta \
       --referencia data/P01308_INS_HUMAN.fasta --salida resultados
python -m pytest -q tests
```

El pipeline guarda en `resultados/` las hebras nuevas, el ARNm y la proteína en FASTA, junto con `pipeline.log`.

---

## 1. Replicación del ADN

**Hebras nuevas.** La replicación es semiconservativa: se separan las dos hebras y cada una sirve de molde para sintetizar su complementaria en sentido 5'→3'. Se obtienen dos moléculas hijas idénticas a la parental, cada una con una hebra antigua y una nueva.

| | Hebra parental (molde) | Hebra nueva |
|---|---|---|
| Hija 1 | 5'-ATG CCG TTA GCT-3' | 3'-TAC GGC AAT CGA-5' (= 5'-AGCTAACGGCAT-3') |
| Hija 2 | 3'-TAC GGC AAT CGA-5' | 5'-ATG CCG TTA GCT-3' |

Si la horquilla avanza de izquierda a derecha, la hebra nueva que se forma sobre el molde inferior crece hacia la horquilla. Es la **conductora** y se sintetiza de forma continua. La que se forma sobre el molde superior crece en sentido contrario. Es la **retardada** y se sintetiza en fragmentos de Okazaki.

**Función de las enzimas.**

- **Helicasa:** rompe los puentes de hidrógeno entre las bases y abre la horquilla gastando ATP. Las proteínas SSB mantienen separadas las hebras sencillas, y la topoisomerasa alivia el superenrollamiento que se acumula por delante.
- **Primasa:** es una ARN polimerasa que sintetiza un cebador corto de ARN, necesario porque la ADN polimerasa no puede empezar una cadena desde cero. Hace falta un cebador en la hebra conductora y uno por cada fragmento de Okazaki.
- **ADN polimerasa:** añade desoxinucleótidos al extremo 3'-OH siguiendo la complementariedad A–T/G–C, siempre en sentido 5'→3'. Su actividad exonucleasa 3'→5' corrige errores. En *E. coli*, la Pol III sintetiza y la Pol I sustituye los cebadores por ADN.
- **Ligasa:** sella las mellas entre fragmentos de Okazaki formando el último enlace fosfodiéster.

**Reflexión: error no corregido.** La polimerasa por sí sola se equivoca aproximadamente una vez cada 10⁴–10⁵ nucleótidos. La corrección de pruebas y la reparación de desapareamientos (MMR) bajan esa tasa a unos 10⁻⁹–10⁻¹⁰ por nucleótido y ronda. Si un error escapa a ambos mecanismos, en la ronda siguiente la hebra errónea sirve de molde y el cambio queda fijado con su pareja correcta: 1 de las 4 moléculas lleva ya una **mutación** estable y heredable. Su efecto depende del codón afectado. Puede ser silenciosa (GCT→GCC, Ala), de cambio de sentido (TTA→TCA, Leu→Ser) o sin sentido. `ej1_replicacion.py` simula este último caso: TTA→TAA convierte el péptido `MPLA` en `MP`. Una inserción o deleción desplazaría además la pauta de lectura.

**Código.** `ej1_replicacion.py` genera la complementaria con `Seq.reverse_complement()` y confirma que las dos hebras nuevas coinciden con las obtenidas a mano.

## 2. Transcripción del ADN a ARN

**Cadena molde.** La ARN polimerasa lee el molde en sentido 3'→5' y sintetiza el ARN en sentido 5'→3'. Para que el transcrito empiece por el codón de inicio AUG, el molde debe ser la hebra inferior, 3'-TAC GGA CTT ACG-5'. La superior es la hebra **codificante**: tiene la misma secuencia que el ARN, con T en lugar de U.

**Transcrito:** 5'-AUG CCU GAA UGC-3', que se traduciría como Met–Pro–Glu–Cys. No tiene codón de paro, así que es un fragmento.

**Promotor y región codificante.** Los 12 nucleótidos del enunciado son **región codificante**: empiezan en el ATG, que marca el inicio de la *traducción*. El **promotor** no aparece en la secuencia. Estaría aguas arriba (en el lado 5' de la hebra codificante), antes del sitio de inicio de la transcripción (+1). En bacterias está formado por las cajas −10 (TATAAT) y −35 (TTGACA); en eucariotas, por la caja TATA (hacia −25/−30) y otros elementos que reconoce la ARN polimerasa II. Entre el +1 y el ATG quedaría la 5'UTR.

**Experimento de orientación** (`ej2_transcripcion.py`, que lee un FASTA). `Seq.transcribe()` solo cambia T por U, es decir, supone que se le da la hebra codificante. Si se le pasa otra cosa, el resultado es incorrecto:

| Entrada | ARN obtenido | Traducción |
|---|---|---|
| Codificante 5'→3' (correcto) | AUGCCUGAAUGC | MPEC |
| Molde 5'→3' sin invertir ni complementar | GCAUUCAGGCAU | AFRH |
| Complemento sin invertir (molde leído 3'→5') | UACGGACUUACG | YGLT |

Por eso, con la hebra molde hay que hacer `reverse_complement()` antes de transcribir. El script lo hace con la opción `--hebra molde`.

## 3. Traducción del ARNm a proteína

En 5'-AUG UAU GCU UAA-3', el codón de inicio es **AUG** (Met) y el de paro es **UAA** (ocre). El resultado es el tripéptido **Met–Tyr–Ala** (`MYA`). `Bio.Seq.translate()` devuelve `MYA*` y, con `to_stop=True`, `MYA`, igual que a mano.

**AUG→GUG.** En eucariotas, el ribosoma que recorre el ARNm pasaría de largo por el GUG y empezaría en el siguiente AUG. Aquí ese AUG está en la posición 5, fuera de fase (GUG U**AUG** CUU AA), así que se leería otra pauta (Met–Leu…) y no se obtendría el péptido original. Biopython lo refleja: con la tabla estándar y `cds=True` rechaza GUG como inicio. En bacterias, en cambio, GUG es un codón de inicio alternativo que usa una minoría de genes. El ARNt iniciador sigue colocando (f)Met, de modo que con la tabla 11 se obtiene `MYA`, aunque con menor eficiencia de inicio.

**Pérdida del codón de paro** (p. ej. UAA→CAA, Gln). El ribosoma sigue leyendo la 3'UTR hasta el siguiente paro en fase y produce una proteína con el extremo C alargado. Con una 3'UTR hipotética, el script obtiene `MYAQAG`. Si no aparece ningún paro antes de la cola poli(A), el ARNm se degrada por la vía *non-stop decay*. La proteína alargada puede plegarse mal o ser degradada.

## 4. Splicing alternativo

Se considera un gen de 5 exones con los exones 1 y 5 constitutivos y tres combinaciones alternativas además de la completa. Para probarlas, `ej4_splicing.py` usa un gen de juguete con exones de 24, 20, 22, 21 y 24 nt (el script enumera las 8 isoformas posibles):

| Isoforma | Pauta de lectura | Proteína resultante |
|---|---|---|
| 1-2-3-4-5 | correcta | completa, 36 aa (referencia) |
| 1-2-3-5 | correcta (21 nt ≡ 0 mod 3) | 29 aa: le faltan los 7 aa del exón 4 y el resto no cambia |
| 1-2-4-5 | desplazada (se salta 22 nt) | 27 aa: extremo C distinto y paro prematuro en el último exón |
| 1-3-5 | desplazada (se saltan 41 nt) | 15 aa: truncada; el paro cae junto a la última unión |

**Diferencias esperables.** Si el fragmento saltado es múltiplo de 3, se obtiene una proteína más corta a la que le falta un segmento, por ejemplo un dominio de unión, una señal de localización o un dominio transmembrana (que daría una forma soluble). Si no lo es, la pauta se desplaza, cambia todo lo que queda aguas abajo y suele aparecer un paro prematuro. Si ese paro queda a más de 50–55 nt de la última unión exón–exón, el ARNm se degrada por NMD. Si no, se produce una proteína truncada. En el gen de juguete, los paros quedan cerca de la última unión y ninguna isoforma activaría el NMD.

**¿Por qué aumenta la diversidad sin más genes?** Con 3 exones opcionales hay ya 2³ = 8 isoformas posibles. A eso se suman los exones mutuamente excluyentes, los sitios 5'/3' alternativos y la retención de intrones. En humanos hay unos 20 000 genes codificantes, y más del 90 % de los que tienen varios exones presentan splicing alternativo. Cada tejido y cada etapa del desarrollo elige isoformas mediante reguladores como las proteínas SR y las hnRNP, de modo que un mismo gen da proteínas con funciones distintas.

**FGFR2 en Ensembl.** El gen (ENSG00000066468, cromosoma 10, hebra −) tiene 59 transcritos, de los que 38 codifican proteína, con 31 longitudes distintas (de 56 a 822 aa). El canónico (MANE Select) es FGFR2-206 (ENST00000358487), con 18 exones y 821 aa. Uno de los más largos, FGFR2-215, tiene 822 aa. Sus isoformas más estudiadas se diferencian en dos exones mutuamente excluyentes que codifican la mitad del tercer dominio Ig: **IIIb**, epitelial, que une FGF7 y FGF10, y **IIIc**, mesenquimal, que une FGF2 entre otros. Un solo cambio de exón altera qué ligandos reconoce el receptor y, por tanto, su papel en el desarrollo. Las mutaciones en esta región se asocian a craneosinostosis (síndromes de Crouzon y Apert). La tabla completa de transcritos se obtiene con `python src/ej4_splicing.py --ensembl FGFR2`.

## 5. Introducción a las proteínas

En Met–Ile–Ser–Gly–Val–Lys–His, el **extremo N** es la **Met** (grupo α-amino libre) y el **extremo C** es la **His** (grupo α-carboxilo libre). Por convenio, los péptidos se escriben de N a C, que es el orden en que los sintetiza el ribosoma. `ej5_proteinas.py` clasifica los restos con la escala de Kyte–Doolittle: Met, Ile y Val son hidrofóbicos, Ser es polar, Gly es especial y Lys e His son básicos. El péptido tiene GRAVY = +0,33 y pI ≈ 8,5.

**Orden y estructura.** La secuencia contiene la información necesaria para el plegamiento (Anfinsen). El efecto hidrofóbico entierra los restos apolares en el núcleo y deja los polares expuestos. Además, el patrón de restos favorece unas estructuras u otras: Ala, Leu, Glu o Met tienden a formar hélices α; Val, Ile o Tyr, hebras β; y Gly y Pro las interrumpen. Con los mismos aminoácidos en otro orden, la proteína se pliega de otra manera.

**Mutación hidrofóbico→hidrofílico en el núcleo.** Enterrar una carga o un grupo polar sin pareja cuesta energía y rompe el empaquetamiento de van der Waals. La proteína se desestabiliza y puede plegarse mal, agregarse, ser degradada por el proteasoma o perder su función. En el script, la mutación I2K hace que el GRAVY pase de +0,33 a −0,87.

**PDB.** En la desoxihemoglobina humana (4HHB), cada subunidad adopta el plegamiento globina de 8 hélices α (A–H) y no tiene láminas β. En cambio, la GFP (1EMA) es un barril de 11 hebras β. Una sola mutación puntual en la cadena β de la hemoglobina, Glu6Val, crea en la superficie una zona hidrofóbica que hace polimerizar la desoxi-HbS (anemia falciforme). Es el caso inverso al del núcleo, pero muestra el mismo principio. `ej5_proteinas.py --pdb` resume las hélices y hebras de cualquier fichero descargado del PDB.

## 6. Actividad integradora: del ADN a la proteína

Se ha usado la **CDS de la insulina humana** (*INS*, transcrito ENST00000381330.5 de Ensembl, GRCh38): 333 nt, 64,6 % de GC, ATG…TAG y 111 codones. `ej6_pipeline.py` encadena los tres procesos y va informando de cada paso por pantalla y en `resultados/pipeline.log`:

```
0. Lectura (FASTA, validación, GC, ¿CDS completa?)
   → 1. Replicación (2 hebras nuevas, 333/333 pares Watson-Crick)
   → 2. Transcripción (molde 3'→5' → ARNm 5'→3')
   → 3. Traducción (ORF AUG…UAG → 110 aa)
   → Comparación con UniProt P01308 (110/110 idénticos)
```

**Resultados.** (1) Las dos hebras nuevas son 5'-CTAGTTGCAG…GGGCCAT-3', sintetizada sobre el molde codificante, y 5'-ATGGCCCTGT…CAACTAG-3', sintetizada sobre el molde. Las dos moléculas hijas son idénticas a la parental. (2) El ARNm es 5'-AUGGCCCUGUGG…UGCAACUAG-3'. (3) La proteína es la **preproinsulina** (110 aa, 12,0 kDa, pI 5,22) y coincide al 100 % con UniProt P01308. El dogma central termina aquí, pero la hormona activa no: tras la traducción se elimina el péptido señal (1–24) y se corta el péptido C (57–87), y las cadenas B (25–54) y A (90–110) quedan unidas por puentes disulfuro. La secuencia del gen no basta para describir la proteína funcional.

**Reflexión: ¿qué punto es más vulnerable?** Depende de si se mira la frecuencia o las consecuencias. La tasa de error es mayor en la traducción (~10⁻⁴–10⁻³ por codón) y en la transcripción (~10⁻⁵–10⁻⁴) que en la replicación (~10⁻⁹–10⁻¹⁰). Sin embargo, esos errores son transitorios: afectan a una molécula de ARN o de proteína que se sustituye. Un error de **replicación** no corregido, en cambio, se vuelve una mutación permanente, pasa a las células hijas y altera *todas* las copias de ARNm y de proteína que se produzcan a partir de él. Por eso la replicación es el punto más vulnerable en cuanto a efecto sobre la función, y también el que tiene más mecanismos de control. En la insulina, por ejemplo, hay mutaciones del gen *INS* que impiden que la proinsulina se pliegue bien y causan diabetes neonatal.

## Referencias

1. Cock, P. J. A. *et al.* (2009). Biopython. *Bioinformatics* 25, 1422–1423.
2. Watson, J. D. y Crick, F. H. C. (1953). Genetical implications of the structure of DNA. *Nature* 171, 964–967.
3. Meselson, M. y Stahl, F. W. (1958). The replication of DNA in *E. coli*. *PNAS* 44, 671–682.
4. Wang, E. T. *et al.* (2008). Alternative isoform regulation in human tissue transcriptomes. *Nature* 456, 470–476.
5. Pan, Q. *et al.* (2008). Deep surveying of alternative splicing complexity in the human transcriptome. *Nat. Genet.* 40, 1413–1415.
6. Harrison, P. W. *et al.* (2024). Ensembl 2024. *Nucleic Acids Res.* 52, D891–D899.
7. Anfinsen, C. B. (1973). Principles that govern the folding of protein chains. *Science* 181, 223–230.
8. Berman, H. M. *et al.* (2000). The Protein Data Bank. *Nucleic Acids Res.* 28, 235–242.
9. UniProt Consortium. Entrada P01308 (INS_HUMAN). https://www.uniprot.org/uniprotkb/P01308
