"""Todo lo que tiene que ver con mostrar la simulación en la terminal:
colores, títulos, pausas y el dibujo de hebras alineadas."""
import os
import sys
import time

# Ritmo de la simulación:
#   "paso"    -> se detiene y espera a que el usuario pulse Enter
#   "auto"    -> avanza sola, con una pequeña espera entre pasos
#   "resumen" -> no dibuja los pasos intermedios, solo los resultados
MODO = "paso"
RETARDO = 1.0      # segundos entre pasos en el modo "auto"
ANCHO = 60         # bases por línea al dibujar hebras largas

USAR_COLOR = True
COLORES = {
    "rojo": "\033[91m", "verde": "\033[92m", "amarillo": "\033[93m", "azul": "\033[94m",
    "magenta": "\033[95m", "cian": "\033[96m", "gris": "\033[90m", "negrita": "\033[1m",
}
FIN_COLOR = "\033[0m"


def activar_colores(activar):
    """Activa los colores solo si se piden y la salida es una terminal de verdad."""
    global USAR_COLOR
    USAR_COLOR = activar and sys.stdout.isatty()
    if USAR_COLOR and os.name == "nt":
        os.system("")   # en Windows, esto activa los códigos de color de la consola


def color(texto, nombre):
    if not USAR_COLOR or nombre is None:
        return texto
    return COLORES[nombre] + texto + FIN_COLOR


def titulo(texto):
    print()
    print(color("=" * 72, "cian"))
    print(color("  " + texto, "negrita"))
    print(color("=" * 72, "cian"))


def subtitulo(texto):
    print()
    print(color(">> " + texto, "amarillo"))


def abreviar(secuencia, maximo=70):
    """Acorta una secuencia larga para los resúmenes: principio...final."""
    if len(secuencia) <= maximo:
        return secuencia
    mitad = (maximo - 3) // 2
    return f"{secuencia[:mitad]}...{secuencia[-mitad:]} ({len(secuencia)} nt)"


def mostrar_detalle():
    """¿Hay que dibujar los pasos intermedios?"""
    return MODO != "resumen"


def pausa():
    if MODO == "paso":
        input(color("   [Enter] para continuar...", "gris"))
    elif MODO == "auto":
        time.sleep(RETARDO)


def regla(inicio, fin):
    """Línea con la posición (empezando en 1) cada 10 bases."""
    casillas = [" "] * (fin - inicio)
    for i in range(inicio, fin):
        numero = str(i + 1)
        if i % 10 == 0 and i - inicio + len(numero) <= len(casillas):
            for j, cifra in enumerate(numero):
                casillas[i - inicio + j] = cifra
    return "".join(casillas)


def dibujar(filas, longitud):
    """Dibuja varias hebras alineadas, una debajo de otra.

    Cada fila es (etiqueta, extremo_izquierdo, casillas, extremo_derecho), donde
    'casillas' es una lista con un texto por posición (puede llevar color).
    Si la molécula es larga, se dibuja en trozos de ANCHO posiciones.
    """
    for inicio in range(0, longitud, ANCHO):
        fin = min(inicio + ANCHO, longitud)
        print(color(f"   {'posición':<17}   {regla(inicio, fin)}", "gris"))
        for etiqueta, izquierda, casillas, derecha in filas:
            izq = izquierda if inicio == 0 else "  "
            der = derecha if fin == longitud else ""
            print(f"   {etiqueta:<17}{izq:>2} {''.join(casillas[inicio:fin])} {der}")
        print()
