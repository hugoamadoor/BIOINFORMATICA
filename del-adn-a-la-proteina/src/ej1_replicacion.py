"""Ejercicio 1. Replicación del ADN.

Uso:  python src/ej1_replicacion.py [SECUENCIA_5a3]
Por defecto usa la hebra superior del enunciado: 5'-ATG CCG TTA GCT-3'.
"""
import sys

from dogma import (complementaria, complemento_alineado, duplex, replicar,
                   traducir, transcribir_desde_codificante, validar)

ENUNCIADO_SUPERIOR = "ATG CCG TTA GCT"
ENUNCIADO_INFERIOR_3a5 = "TAC GGC AAT CGA"          # tal como aparece en el enunciado
MANUAL_NUEVA_SOBRE_SUPERIOR = "AGCTAACGGCAT"        # resultado a mano, 5'->3'
MANUAL_NUEVA_SOBRE_INFERIOR = "ATGCCGTTAGCT"        # resultado a mano, 5'->3'


def comprobar(nombre, obtenido, esperado):
    estado = "OK" if obtenido == esperado else "NO COINCIDE"
    print(f"  {nombre:<32} {obtenido}  (manual: {esperado})  -> {estado}")
    return obtenido == esperado


def main():
    superior = validar(sys.argv[1] if len(sys.argv) > 1 else ENUNCIADO_SUPERIOR)
    print("Molécula parental:")
    print(duplex(superior), "\n")

    if superior == validar(ENUNCIADO_SUPERIOR):
        print("Comprobación del enunciado: la hebra inferior dada es complementaria ->",
              complemento_alineado(superior) == validar(ENUNCIADO_INFERIOR_3a5), "\n")

    r = replicar(superior)
    print("Hebras nuevas (escritas 5'->3'):")
    print(f"  sobre el molde superior : 5'-{r['nueva_sobre_superior']}-3'")
    print(f"  sobre el molde inferior : 5'-{r['nueva_sobre_inferior']}-3'\n")

    print("Moléculas hijas (parental | nueva):")
    print("  Hija 1: parental superior + nueva complementaria")
    print("  " + duplex(r["hija_1"][0]).replace("\n", "\n  "))
    print("  Hija 2: nueva superior + parental inferior")
    print("  " + duplex(r["hija_2"][0]).replace("\n", "\n  "), "\n")

    print("Supuesto: la horquilla avanza de izquierda a derecha (desde el extremo 5' de la"
          " hebra superior).")
    print("  - Molde inferior (3'->5' en ese sentido): hebra CONDUCTORA, síntesis continua.")
    print("  - Molde superior (5'->3' en ese sentido): hebra RETARDADA, fragmentos de Okazaki.\n")

    if superior == validar(ENUNCIADO_SUPERIOR):
        print("Comparación con el resultado manual:")
        comprobar("nueva sobre molde superior", r["nueva_sobre_superior"],
                  MANUAL_NUEVA_SOBRE_SUPERIOR)
        comprobar("nueva sobre molde inferior", r["nueva_sobre_inferior"],
                  MANUAL_NUEVA_SOBRE_INFERIOR)
        print()

        # Reflexión: error no corregido de la polimerasa (T->A en la posición 8)
        mutada = superior[:7] + "A" + superior[8:]
        print("Simulación de un error no corregido (T8A en la hebra nueva):")
        print(f"  Ronda 1: la hebra nueva lleva {mutada} frente a un molde con T -> desapareamiento A·A")
        print("  Ronda 2: al usarse como molde, el error se copia con su pareja correcta (A-T)")
        print("           -> 1 de las 4 moléculas lleva ya una mutación estable y heredable")
        for nombre, s in (("original", superior), ("mutante ", mutada)):
            prot = traducir(transcribir_desde_codificante(s))
            print(f"  {nombre}: {s} -> proteína {prot or '(vacía)'}"
                  f"{'  (codón de paro prematuro: mutación sin sentido)' if len(prot) < len(s)//3 else ''}")
        print(f"  Complementaria de la mutante: {complementaria(mutada)}")


if __name__ == "__main__":
    main()
