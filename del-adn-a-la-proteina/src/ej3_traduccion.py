"""Ejercicio 3. Traducción del ARNm a proteína.

Uso:  python src/ej3_traduccion.py [ARNM_5a3]
Por defecto usa el transcrito del enunciado: 5'-AUG UAU GCU UAA-3'.
"""
import sys

from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError

from dogma import buscar_orf, validar

ENUNCIADO = "AUG UAU GCU UAA"
MANUAL = "MYA"


def codones(s):
    return " ".join(s[i:i + 3] for i in range(0, len(s), 3))


def main():
    arnm = validar(sys.argv[1] if len(sys.argv) > 1 else ENUNCIADO, tipo="ARN")
    inicio, paro = buscar_orf(arnm)
    print(f"ARNm: 5'-{codones(arnm)}-3'")
    print(f"Codón de inicio: {arnm[inicio:inicio+3] if inicio is not None else '-'} (pos. {inicio + 1 if inicio is not None else '-'})")
    print(f"Codón de paro:   {arnm[paro:paro+3] if paro is not None else 'ninguno'}"
          f"{f' (pos. {paro + 1})' if paro is not None else ''}")

    seq = Seq(arnm)
    prot = str(seq.translate(to_stop=True))
    print(f"\nBio.Seq.translate()              -> {seq.translate()}   ('*' = paro)")
    print(f"Bio.Seq.translate(to_stop=True)  -> {prot}")
    tres = [CodonTable.standard_rna_table.forward_table.get(c, "Stop") for c in codones(arnm).split()]
    print("Aminoácidos (código de una letra):", " - ".join(tres))
    if arnm == validar(ENUNCIADO, "ARN"):
        print(f"Coincide con el resultado manual ({MANUAL} = Met-Tyr-Ala): {prot == MANUAL}")

    print("\nMutación 1: AUG -> GUG en el codón de inicio")
    gug = "G" + arnm[1:]
    otro = buscar_orf(gug)[0]
    if otro is None:
        print(f"  {codones(gug)}: no queda ningún AUG")
    else:
        desde = gug[otro:]
        print(f"  {codones(gug)}: el primer AUG está ahora en la pos. {otro + 1}"
              f" ({'en fase' if otro % 3 == 0 else 'fuera de fase'})")
        print(f"  Un ribosoma eucariota que escanee iniciaría ahí: {codones(desde)} -> "
              f"{Seq(desde[:len(desde)//3*3]).translate()} (otra pauta de lectura, sin paro)")
    for tabla, nombre in ((1, "estándar (eucariotas)"), (11, "bacteriana")):
        try:
            p = Seq(gug).translate(table=tabla, cds=True)
            print(f"  Tabla {tabla} {nombre}: GUG válido como inicio -> {p} (el iniciador sigue siendo Met)")
        except TranslationError as e:
            print(f"  Tabla {tabla} {nombre}: {e}")

    print("\nMutación 2: pérdida del codón de paro (UAA -> CAA)")
    sin_paro = arnm[:paro] + "C" + arnm[paro + 1:] if paro is not None else arnm
    print(f"  {codones(sin_paro)} -> {Seq(sin_paro).translate()} "
          "(el ribosoma sigue leyendo la 3'UTR hasta el siguiente paro en fase)")
    utr = "GCUGGAUGAA"   # 3'UTR hipotética para ilustrar la lectura continuada
    print(f"  con una 3'UTR hipotética ({utr}): {Seq(sin_paro + utr[:len(utr)//3*3]).translate(to_stop=True)}"
          " -> proteína alargada")


if __name__ == "__main__":
    main()
