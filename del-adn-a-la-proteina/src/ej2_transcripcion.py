"""Ejercicio 2. Transcripción del ADN a ARN.

Uso:  python src/ej2_transcripcion.py FASTA [--hebra codificante|molde]
      python src/ej2_transcripcion.py data/ej2_codificante.fasta
      python src/ej2_transcripcion.py data/ej2_molde.fasta --hebra molde

El FASTA contiene una única hebra escrita 5'->3'. Con --hebra se indica si es la
codificante (misma secuencia que el ARN) o la molde (complementaria al ARN).
"""
import argparse

from Bio.Seq import Seq

from dogma import (buscar_orf, duplex, leer_fasta, traducir,
                   transcribir_desde_codificante, transcribir_desde_molde)

MANUAL_ARNM = "AUGCCUGAAUGC"


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("fasta")
    p.add_argument("--hebra", choices=["codificante", "molde"], default="codificante")
    a = p.parse_args()

    ident, s = leer_fasta(a.fasta)
    print(f"Registro: {ident} ({len(s)} nt), interpretada como hebra {a.hebra}")

    if a.hebra == "codificante":
        codificante = s
        arnm = transcribir_desde_codificante(s)
    else:
        codificante = str(Seq(s).reverse_complement())
        arnm = transcribir_desde_molde(s)

    print("Dúplex (arriba: codificante; abajo: molde, leída 3'->5' por la ARN polimerasa):")
    print(duplex(codificante))
    print(f"\nARNm 5'->3': {arnm}")
    inicio, paro = buscar_orf(arnm)
    print(f"Primer AUG en la posición {inicio + 1 if inicio is not None else '-'}; "
          f"codón de paro en fase: {'no hay' if paro is None else paro + 1}")
    print(f"Traducción: {traducir(arnm)}")
    if arnm == MANUAL_ARNM:
        print("Coincide con el resultado manual (5'-AUG CCU GAA UGC-3').")

    # Experimento de orientación que pide el enunciado
    print("\nExperimento: ¿qué pasa si se usa la hebra equivocada o sin invertir?")
    casos = {
        "codificante + transcribe()          (correcto)": transcribir_desde_codificante(codificante),
        "molde 5'->3' + transcribe() directo (error)": str(Seq(str(Seq(codificante).reverse_complement())).transcribe()),
        "molde sin invertir (3'->5') + transcribe() (error)": str(Seq(codificante).complement().transcribe()),
    }
    for nombre, arn in casos.items():
        print(f"  {nombre:<52} {arn} -> {traducir(arn, hasta_paro=False)}")


if __name__ == "__main__":
    main()
