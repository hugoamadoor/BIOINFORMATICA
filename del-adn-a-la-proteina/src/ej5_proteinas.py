"""Ejercicio 5. Introducción a las proteínas.

Uso:  python src/ej5_proteinas.py                    # péptido del enunciado
      python src/ej5_proteinas.py --pdb 4HHB.pdb     # resumen de hélices/láminas de un PDB descargado

El péptido Met-Ile-Ser-Gly-Val-Lys-His se escribe, por convenio, de N a C.
"""
import argparse

from Bio.SeqUtils import seq1
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.ProtParamData import kd   # escala de hidropatía de Kyte-Doolittle

PEPTIDO = "Met-Ile-Ser-Gly-Val-Lys-His"
CLASES = {**dict.fromkeys("AVILMFWC", "hidrofóbico"), "G": "especial (cadena lateral = H)",
          "P": "especial", **dict.fromkeys("STNQY", "polar sin carga"),
          **dict.fromkeys("KRH", "básico (+)"), **dict.fromkeys("DE", "ácido (-)")}


def propiedades(seq):
    pa = ProteinAnalysis(seq)
    return pa.molecular_weight(), pa.isoelectric_point(), pa.gravy()


def analizar_peptido():
    tres = PEPTIDO.split("-")
    seq = "".join(seq1(a) for a in tres)
    print(f"Péptido: {PEPTIDO}  ({seq})")
    print(f"Extremo N (amino libre):    {tres[0]} (posición 1)")
    print(f"Extremo C (carboxilo libre): {tres[-1]} (posición {len(tres)})\n")
    print(f"{'Pos':>3} {'Resto':<5} {'Kyte-Doolittle':>14}  Tipo")
    for i, (t, a) in enumerate(zip(tres, seq), 1):
        print(f"{i:>3} {t:<5} {kd[a]:>14.1f}  {CLASES[a]}")
    mw, pi, gravy = propiedades(seq)
    print(f"\nMasa: {mw:.1f} Da | pI: {pi:.2f} | GRAVY: {gravy:+.2f} (>0 hidrofóbico en conjunto)\n")

    print("Mutaciones hidrofóbico -> hidrofílico (como si el resto estuviera enterrado en el núcleo):")
    for pos, nuevo in ((2, "K"), (5, "K"), (5, "S")):
        mut = seq[:pos - 1] + nuevo + seq[pos:]
        _, pi_m, g_m = propiedades(mut)
        print(f"  {seq[pos-1]}{pos}{nuevo}: {mut}  ΔKD del resto = {kd[nuevo]-kd[seq[pos-1]]:+.1f} | "
              f"GRAVY {gravy:+.2f} -> {g_m:+.2f} | pI {pi:.2f} -> {pi_m:.2f}")


def resumen_pdb(ruta):
    """Cuenta los registros HELIX y SHEET de un fichero PDB (formato .pdb)."""
    helices, hebras, titulo = [], [], []
    with open(ruta) as f:
        for linea in f:
            if linea.startswith("TITLE"):
                titulo.append(linea[10:80].strip())
            elif linea.startswith("HELIX"):
                helices.append(int(linea[71:76]))          # longitud de la hélice
            elif linea.startswith("SHEET"):
                hebras.append(int(linea[33:37]) - int(linea[22:26]) + 1)
    print(f"{ruta}: {' '.join(titulo)}")
    print(f"  Hélices: {len(helices)} (media {sum(helices)/max(len(helices),1):.1f} restos)")
    print(f"  Hebras beta: {len(hebras)} (media {sum(hebras)/max(len(hebras),1):.1f} restos)")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pdb", help="fichero .pdb descargado de https://www.rcsb.org")
    a = p.parse_args()
    if a.pdb:
        resumen_pdb(a.pdb)
    else:
        analizar_peptido()


if __name__ == "__main__":
    main()
