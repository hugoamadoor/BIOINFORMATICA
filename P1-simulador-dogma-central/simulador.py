"""Simulador del dogma central de la biología molecular.

Uso:
    python simulador.py                         # menú interactivo
    python simulador.py --ejemplo --modo paso   # gen de ejemplo, paso a paso
    python simulador.py --fasta datos/insulina_INS_cds.fasta --modo auto
    python simulador.py --secuencia ATGGCTTAA --modo resumen

La secuencia de entrada es la hebra codificante (superior), escrita de 5' a 3'.
"""
import argparse

from dogma import pantalla
from dogma.codigo_genetico import UNA_LETRA
from dogma.replicacion import simular_replicacion
from dogma.secuencias import GEN_EJEMPLO, PARTES_GEN_EJEMPLO, leer_fasta, validar
from dogma.transcripcion import simular_transcripcion
from dogma.traduccion import simular_traduccion


def leer_argumentos():
    p = argparse.ArgumentParser(description="Simulador del dogma central (ADN -> ADN -> ARN -> proteína)")
    origen = p.add_mutually_exclusive_group()
    origen.add_argument("--ejemplo", action="store_true", help="usar el gen bacteriano de ejemplo")
    origen.add_argument("--secuencia", help="secuencia de ADN (hebra codificante, 5'->3')")
    origen.add_argument("--fasta", help="fichero FASTA con la secuencia")
    p.add_argument("--modo", choices=["paso", "auto", "resumen"],
                   help="paso: espera a Enter | auto: avanza solo | resumen: solo resultados")
    p.add_argument("--okazaki", type=int, default=20, help="tamaño de los fragmentos de Okazaki (nt)")
    p.add_argument("--cebador", type=int, default=3, help="tamaño de los cebadores de ARN (nt)")
    p.add_argument("--retardo", type=float, default=1.0, help="segundos entre pasos en el modo auto")
    p.add_argument("--ancho", type=int, default=60, help="bases por línea en los dibujos")
    p.add_argument("--sin-color", action="store_true", help="no usar colores en la terminal")
    return p.parse_args()


def elegir_secuencia(args):
    """Devuelve (nombre, secuencia) según los argumentos o preguntando al usuario."""
    if args.ejemplo:
        return "Gen bacteriano de ejemplo", GEN_EJEMPLO
    if args.secuencia:
        return "Secuencia introducida", validar(args.secuencia)
    if args.fasta:
        return leer_fasta(args.fasta)

    print("¿Qué molécula de ADN quieres usar?")
    print("  1) El gen bacteriano de ejemplo (recomendado)")
    print("  2) Escribir una secuencia")
    print("  3) Leer un fichero FASTA")
    while True:
        opcion = input("Opción [1]: ").strip() or "1"
        try:
            if opcion == "1":
                return "Gen bacteriano de ejemplo", GEN_EJEMPLO
            if opcion == "2":
                return "Secuencia introducida", validar(input("Secuencia (5'->3'): "))
            if opcion == "3":
                return leer_fasta(input("Ruta del fichero FASTA: ").strip())
        except (ValueError, OSError) as error:
            print(f"  Error: {error}")
            continue
        print("  Escribe 1, 2 o 3.")


def elegir_modo(args):
    if args.modo:
        return args.modo
    print("\n¿Cómo quieres ver la simulación?")
    print("  1) Paso a paso (pulsando Enter)")
    print("  2) Automática")
    print("  3) Solo el resumen de resultados")
    opcion = input("Opción [1]: ").strip() or "1"
    return {"1": "paso", "2": "auto", "3": "resumen"}.get(opcion, "paso")


def presentar_molecula(nombre, adn):
    pantalla.titulo("MOLÉCULA DE PARTIDA")
    print(f"   {nombre}: {len(adn)} pares de bases")
    if adn == GEN_EJEMPLO:
        posicion = 0
        for parte, secuencia in PARTES_GEN_EJEMPLO:
            print(f"   {posicion + 1:>4}-{posicion + len(secuencia):<4} {parte:<38} {secuencia}")
            posicion += len(secuencia)
    print(f"   Hebra codificante: 5'-{pantalla.abreviar(adn)}-3'")
    pantalla.pausa()


def resumen_final(adn, replicacion, transcripcion, traduccion):
    pantalla.titulo("RESUMEN: FLUJO DE LA INFORMACIÓN GENÉTICA")
    proteina = traduccion["proteina"]
    print(f"   ADN original      5'-{pantalla.abreviar(adn)}-3'")
    fragmentos = len(replicacion["fragmentos"])
    print(f"   ADN -> ADN        2 moléculas hijas idénticas a la original "
          f"({fragmentos} fragmento{'s' if fragmentos != 1 else ''} de Okazaki en la rezagada)")
    print(f"   ADN -> ARN        5'-{pantalla.abreviar(transcripcion['arnm'])}-3'")
    if proteina:
        cadena = "-".join(proteina) if len(proteina) <= 15 else "-".join(proteina[:15]) + "-..."
        print(f"   ARN -> proteína   {cadena}")
        print(f"                     {''.join(UNA_LETRA[aa] for aa in proteina)}  "
              f"({len(proteina)} aminoácidos)")
    else:
        print("   ARN -> proteína   (no se ha formado ninguna proteína)")


def main():
    args = leer_argumentos()
    pantalla.activar_colores(not args.sin_color)
    pantalla.RETARDO = args.retardo
    pantalla.ANCHO = args.ancho

    nombre, adn = elegir_secuencia(args)
    pantalla.MODO = elegir_modo(args)

    presentar_molecula(nombre, adn)
    replicacion = simular_replicacion(adn, args.okazaki, args.cebador)
    # La transcripción se hace sobre la molécula original (las hijas son idénticas a ella).
    transcripcion = simular_transcripcion(adn)
    traduccion = simular_traduccion(transcripcion["arnm"])
    resumen_final(adn, replicacion, transcripcion, traduccion)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError) as error:
        raise SystemExit(f"Error: {error}")
    except KeyboardInterrupt:
        raise SystemExit("\nSimulación interrumpida.")
