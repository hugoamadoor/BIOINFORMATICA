"""Ejercicio 1: porcentaje de GC y traducción de varias secuencias de un FASTA.

Uso:
    python ej1_contenido_gc.py                                   # usa datos/insulina_humana_multiple.fasta
    python ej1_contenido_gc.py mi_fichero.fasta
    python ej1_contenido_gc.py mi_fichero.fasta --completa       # traduce la secuencia entera

Para cada secuencia de ADN:
  1. calcula el contenido de GC en porcentaje;
  2. obtiene la cadena de aminoácidos;
  3. guarda todas las proteínas en resultados/proteinas_ej1.fasta con el mismo
     identificador y descripción que la secuencia original.

Por defecto se traduce el marco abierto de lectura (ORF) más largo: desde un ATG hasta
el primer codón de parada. Así, en un ARNm completo se obtiene la proteína real y no se
traduce por error la región no codificante del principio (5'UTR). Con --completa se
traduce la secuencia entera desde la primera base, codón a codón.
"""
import argparse
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction

from comun import CARPETA_RESULTADOS

FICHERO_POR_DEFECTO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "datos", "insulina_humana_multiple.fasta")


def porcentaje_gc(secuencia):
    return 100 * gc_fraction(secuencia)


def orf_mas_largo(secuencia):
    """Devuelve la proteína del ORF más largo de la hebra (ATG ... parada) en los 3 marcos."""
    mejor = ""
    for marco in range(3):
        trozo = secuencia[marco:]
        trozo = trozo[: len(trozo) // 3 * 3]          # solo codones completos
        proteina = str(trozo.translate())
        # Cada fragmento que acaba en '*' (codón de parada) se lee desde su primera
        # metionina (M). El último trozo no tiene parada, así que no es un ORF completo.
        for fragmento in proteina.split("*")[:-1]:
            inicio = fragmento.find("M")
            if inicio != -1 and len(fragmento) - inicio > len(mejor):
                mejor = fragmento[inicio:]
    return mejor


def traduccion_completa(secuencia):
    """Traduce toda la secuencia desde la primera base ('*' = codón de parada)."""
    secuencia = secuencia[: len(secuencia) // 3 * 3]
    return str(secuencia.translate())


def main():
    p = argparse.ArgumentParser(description="Contenido GC y traducción de un FASTA")
    p.add_argument("fasta", nargs="?", default=FICHERO_POR_DEFECTO)
    p.add_argument("--completa", action="store_true", help="traducir la secuencia entera en vez del ORF más largo")
    a = p.parse_args()

    proteinas = []
    print(f"Fichero: {a.fasta}\n")
    print(f"{'Identificador':<18}{'Longitud':>9}{'GC (%)':>9}{'Proteína (aa)':>15}")
    for registro in SeqIO.parse(a.fasta, "fasta"):
        gc = porcentaje_gc(registro.seq)
        proteina = traduccion_completa(registro.seq) if a.completa else orf_mas_largo(registro.seq)
        print(f"{registro.id:<18}{len(registro.seq):>9}{gc:>9.2f}{len(proteina):>15}")
        proteinas.append(SeqRecord(Seq(proteina), id=registro.id, description=registro.description))

    print()
    for r in proteinas:
        print(f">{r.id}\n{r.seq[:60]}{'...' if len(r.seq) > 60 else ''}")

    salida = os.path.join(CARPETA_RESULTADOS, "proteinas_ej1.fasta")
    SeqIO.write(proteinas, salida, "fasta")
    print(f"\nSe han escrito {len(proteinas)} proteínas en {salida}")


if __name__ == "__main__":
    main()
