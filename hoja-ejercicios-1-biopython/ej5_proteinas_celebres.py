"""Ejercicio 5: información y ficheros FASTA de tres proteínas "célebres".

Uso:  python ej5_proteinas_celebres.py

Proteínas elegidas (todas humanas, de UniProtKB/Swiss-Prot, la base de datos curada):
  - Hemoglobina, subunidad beta (P68871)
  - Mioglobina (P02144)
  - Colágeno tipo I, cadena alfa-1 (P02452)

Para cada una se descarga su FASTA de UniProt (y, si UniProt no responde, del NCBI con
Entrez), se guarda en resultados/ y se muestran su información y algunas propiedades
calculadas con Biopython.
"""
import os
import urllib.request

from Bio import Entrez, SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from comun import CARPETA_RESULTADOS, consultar_ncbi

PROTEINAS = {
    "P68871": {
        "nombre": "Hemoglobina (subunidad beta)",
        "gen": "HBB",
        "funcion": "Transporta el oxígeno desde los pulmones a los tejidos dentro de los glóbulos rojos. "
                   "La hemoglobina adulta es un tetrámero de dos cadenas alfa y dos beta, cada una con un grupo hemo.",
        "curiosidad": "El cambio de un solo aminoácido (Glu6Val) produce la anemia falciforme.",
    },
    "P02144": {
        "nombre": "Mioglobina",
        "gen": "MB",
        "funcion": "Almacena oxígeno en el músculo y lo cede cuando hace falta. Es una sola cadena con un grupo hemo.",
        "curiosidad": "Fue la primera proteína cuya estructura 3D se resolvió por difracción de rayos X (Kendrew, 1958).",
    },
    "P02452": {
        "nombre": "Colágeno tipo I, cadena alfa-1",
        "gen": "COL1A1",
        "funcion": "Proteína estructural de piel, hueso y tendones. Tres cadenas se enrollan en una triple hélice.",
        "curiosidad": "Su secuencia repite el patrón Gly-X-Y, así que la glicina es casi un tercio de la región helicoidal. "
                      "Mutaciones en COL1A1 causan osteogénesis imperfecta.",
    },
}


def descargar_fasta(acceso):
    """Descarga el FASTA de UniProt; si falla, lo pide al NCBI (base de datos 'protein')."""
    ruta = os.path.join(CARPETA_RESULTADOS, f"{acceso}.fasta")
    if os.path.exists(ruta):
        return ruta
    try:
        url = f"https://rest.uniprot.org/uniprotkb/{acceso}.fasta"
        with urllib.request.urlopen(url, timeout=30) as r:
            texto = r.read().decode()
        origen = "UniProt"
    except OSError:
        handle = consultar_ncbi(Entrez.efetch, db="protein", id=acceso, rettype="fasta", retmode="text")
        texto = handle.read()
        handle.close()
        origen = "NCBI"
    with open(ruta, "w", encoding="utf-8") as f:
        f.write(texto)
    print(f"  (descargado de {origen})")
    return ruta


def main():
    for acceso, info in PROTEINAS.items():
        print("=" * 78)
        print(f"{info['nombre']}  ·  gen {info['gen']}  ·  UniProt {acceso}")
        print("=" * 78)
        registro = SeqIO.read(descargar_fasta(acceso), "fasta")
        analisis = ProteinAnalysis(str(registro.seq))
        glicina = 100 * str(registro.seq).count("G") / len(registro.seq)

        print(f"Cabecera FASTA: {registro.description}")
        print(f"Función: {info['funcion']}")
        print(f"Dato curioso: {info['curiosidad']}")
        print(f"Longitud: {len(registro.seq)} aminoácidos | masa: {analisis.molecular_weight() / 1000:.1f} kDa | "
              f"pI: {analisis.isoelectric_point():.2f} | glicina: {glicina:.1f} %")
        print(f"Secuencia (inicio): {registro.seq[:60]}...\n")

    print(f"FASTA guardados en {CARPETA_RESULTADOS}")


if __name__ == "__main__":
    main()
