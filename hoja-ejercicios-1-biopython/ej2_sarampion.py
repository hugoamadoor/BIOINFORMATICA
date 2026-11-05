"""Ejercicio 2: descarga automática del genoma oficial del virus del sarampión.

Uso:  python ej2_sarampion.py

1. Entrez.esearch en la base de datos "nucleotide" con el término
   "Measles morbillivirus[Organism]", limitado a la base de datos secundaria (RefSeq,
   datos curados) con "AND srcdb_refseq[Properties]".
2. Se toma el primer identificador de la lista (IdList[0]).
3. Entrez.efetch lo descarga en formato FASTA.
4. Se muestra la descripción del registro y los 100 primeros nucleótidos.
"""
from Bio import Entrez, SeqIO

from comun import consultar_ncbi, descargar

TERMINO = "Measles morbillivirus[Organism] AND srcdb_refseq[Properties]"


def main():
    print(f"Buscando en NCBI Nucleotide: {TERMINO}")
    handle = consultar_ncbi(Entrez.esearch, db="nucleotide", term=TERMINO)
    resultado = Entrez.read(handle)
    handle.close()
    print(f"Registros encontrados: {resultado['Count']}")

    if not resultado["IdList"]:
        print("La búsqueda no ha devuelto ningún registro.")
        return
    primer_id = resultado["IdList"][0]
    print(f"Primer identificador (IdList[0]): {primer_id}\n")

    ruta = descargar("nucleotide", primer_id, "fasta", f"sarampion_{primer_id}.fasta")
    registro = SeqIO.read(ruta, "fasta")

    print(f"Descripción: {registro.description}")
    print(f"Longitud del genoma: {len(registro.seq)} nucleótidos")
    print(f"Primeros 100 nucleótidos:\n{registro.seq[:100]}")


if __name__ == "__main__":
    main()
