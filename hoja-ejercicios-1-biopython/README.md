# Hoja de Ejercicios 1: Biopython, ficheros y bases de datos

Universidad de Las Palmas de Gran Canaria · Grado en Ciencia e Ingeniería de Datos  
Bioinformática (prácticas de laboratorio) · Curso 2026/2027  
**Grupo 19:** Hugo Amador Hidalgo y Sara Lillo  
Repositorio: https://github.com/hugoamadoor/BIOINFORMATICA/tree/main/hoja-ejercicios-1-biopython

Cinco scripts en Python con Biopython para leer y escribir ficheros FASTA y GenBank y
para buscar y descargar secuencias de bases de datos biológicas (NCBI y UniProt).

## Instalación y uso

```bash
pip install -r requirements.txt
python ej1_contenido_gc.py
python ej2_sarampion.py
python ej3_covid19.py
python ej4_bacterias.py                      # o: --buscar / --acceso NC_001371
python ej5_proteinas_celebres.py
python -m pytest -q tests
```

Los ejercicios 2 a 5 necesitan conexión a Internet. El NCBI pide un correo de contacto:
cámbialo en `comun.py` o defínelo antes de ejecutar con `export NCBI_EMAIL=tu_correo@ulpgc.es`.
Los ficheros descargados se guardan en `resultados/` y no se vuelven a descargar.

## Ejercicios

| Script | Ejercicio | Qué hace |
|---|---|---|
| `ej1_contenido_gc.py` | 1. Porcentaje GC | Lee un FASTA con varias secuencias de ADN, calcula el % de GC de cada una, las traduce y guarda las proteínas en `resultados/proteinas_ej1.fasta` con los identificadores y descripciones originales. |
| `ej2_sarampion.py` | 2. El sarampión | `Entrez.esearch` con `Measles morbillivirus[Organism] AND srcdb_refseq[Properties]` (base de datos secundaria RefSeq), toma `IdList[0]`, lo descarga en FASTA con `Entrez.efetch` y muestra la descripción y los 100 primeros nucleótidos. |
| `ej3_covid19.py` | 3. La COVID-19 | Descarga el GenBank de `NC_045512`, recorre sus *features*, filtra las de tipo `gene` y muestra un listado numerado con el nombre (`qualifiers.get("gene")`) y las coordenadas de cada gen. |
| `ej4_bacterias.py` | 4. Bacterias | Descarga del NCBI el plásmido pT181 de *S. aureus* (`NZ_CP167936.1`) en GenBank, localiza todas las CDS, las traduce con la tabla bacteriana (11), comprueba que coinciden con la traducción anotada, guarda las proteínas en `resultados/` e indica cuáles son de resistencia a antibióticos. |
| `ej5_proteinas_celebres.py` | 5. Proteínas célebres | Descarga de UniProt los FASTA de la hemoglobina beta (P68871), la mioglobina (P02144) y el colágeno tipo I alfa-1 (P02452) y muestra su función, un dato curioso y propiedades calculadas (longitud, masa, pI, % de glicina). |

### Decisiones tomadas

- **Ejercicio 1, qué se traduce.** Las secuencias del fichero de ejemplo son ARNm completos
  de la insulina, con regiones no codificantes al principio y al final. Por eso, por
  defecto se traduce el **ORF más largo** (de un ATG al primer codón de parada en los tres
  marcos de lectura), que en los cuatro casos da la preproinsulina de 110 aminoácidos.
  Con la opción `--completa` se traduce la secuencia entera desde la primera base.
- **Código de acceso del ejercicio 4.** El enunciado indica `NC_001371`, pero en el NCBI ese
  registro es el plásmido ColE1 de *Escherichia coli*, no pT181 de *S. aureus* (no tiene genes
  de resistencia). Por eso el script analiza por defecto `NZ_CP167936.1`, el registro de pT181
  de *S. aureus* que devuelve la búsqueda por nombre y organismo; se deja fijo para que el
  resultado sea siempre el mismo. Con `--buscar` se repite la búsqueda en el NCBI y con
  `--acceso NC_001371` se analiza el registro del enunciado.
- **Descarga de GenBank.** Se pide con `rettype="gbwithparts"`, porque con `"gb"` algunos
  registros RefSeq llegan sin la secuencia.
- **Coordenadas (ejercicios 3 y 4).** Biopython cuenta desde 0 y no incluye el final; se
  muestran como en GenBank, desde 1 y con el final incluido.
- **Hebra (ejercicio 4).** `feature.extract()` tiene en cuenta si la CDS está en la hebra
  complementaria antes de traducirla.
- **Errores del NCBI.** El servidor del NCBI devuelve a veces errores temporales
  (`HTTP Error 500`) aunque la consulta sea correcta. Todas las llamadas a Entrez pasan por
  `consultar_ncbi()` (en `comun.py`), que espera unos segundos y reintenta hasta cuatro veces;
  si el NCBI sigue sin responder, muestra un mensaje claro en lugar del error de Python.
- **Certificados.** Si `truststore` está instalado, Python usa los certificados del
  sistema, lo que evita errores HTTPS en redes con antivirus o proxy.

## Estructura

```
comun.py                    correo del NCBI, reintentos, descargas con Entrez.efetch y carpeta de resultados
ej1_contenido_gc.py … ej5_proteinas_celebres.py
datos/insulina_humana_multiple.fasta   4 variantes de ARNm de la insulina humana (ejemplo de clase)
resultados/                 ficheros descargados y generados
tests/test_ejercicios.py    pruebas que no necesitan Internet
```
