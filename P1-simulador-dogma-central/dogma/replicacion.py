"""Replicación del ADN (ADN -> ADN).

Modelo simplificado de una horquilla de replicación bacteriana (E. coli) que sale del
origen, en el extremo izquierdo de la molécula, y avanza hacia la derecha:

  - La HELICASA abre la doble hélice por bloques. Las proteínas SSB cubren las hebras
    sencillas y la TOPOISOMERASA (girasa) alivia la tensión por delante de la horquilla.
  - La PRIMASA coloca cebadores cortos de ARN (se dibujan en minúscula y con u).
  - La ADN POLIMERASA III solo sintetiza en sentido 5'->3':
      * Sobre la hebra inferior (3'->5' de izquierda a derecha) la hebra nueva crece
        en el mismo sentido en que avanza la horquilla: es la CADENA LÍDER, continua,
        y le basta un cebador.
      * Sobre la hebra superior (5'->3') la hebra nueva tiene que crecer hacia atrás:
        es la CADENA REZAGADA, que se fabrica a trozos (FRAGMENTOS DE OKAZAKI), cada
        uno con su propio cebador.
  - La ADN POLIMERASA I sustituye el ARN de los cebadores por ADN.
  - La LIGASA sella las mellas que quedan entre fragmentos.
"""
from . import pantalla
from .secuencias import COMPLEMENTO_ADN, COMPLEMENTO_ARN, complementaria

# Colores con los que se dibujan los fragmentos de Okazaki (se van alternando).
COLORES_FRAGMENTOS = ["amarillo", "azul", "magenta"]


def pintar(hebra, fragmento_de, color_normal):
    """Convierte una hebra en casillas de colores para dibujarla.

    ' ' = aún no se ha abierto | '~' = hebra sencilla cubierta por SSB |
    minúscula = cebador de ARN | mayúscula = ADN nuevo
    """
    casillas = []
    for i, letra in enumerate(hebra):
        if letra == "~":
            casillas.append(pantalla.color("~", "gris"))
        elif letra.islower():
            casillas.append(pantalla.color(letra, "rojo"))
        elif fragmento_de is not None and fragmento_de[i] is not None:
            casillas.append(pantalla.color(letra, COLORES_FRAGMENTOS[fragmento_de[i] % 3]))
        else:
            casillas.append(pantalla.color(letra, color_normal))
    return casillas


def dibujar_horquilla(superior, inferior, lider, rezagada, fragmento_de, abierto_hasta):
    marcas = [" "] * len(superior)
    if abierto_hasta < len(superior):
        marcas[abierto_hasta] = pantalla.color("^", "cian")
    pantalla.dibujar([
        ("Sup. parental", "5'", list(superior), "3'"),
        ("Rezagada (nueva)", "3'", pintar(rezagada, fragmento_de, None), "5'"),
        ("Líder (nueva)", "5'", pintar(lider, None, "verde"), "3'"),
        ("Inf. parental", "3'", list(inferior), "5'"),
        ("", "", marcas, ""),
    ], len(superior))
    if abierto_hasta < len(superior):
        print(f"   ^ helicasa en la posición {abierto_hasta + 1}; a su derecha la doble hélice "
              "sigue cerrada\n     y la topoisomerasa alivia la tensión por delante.")
    else:
        print("   La horquilla ha llegado al final de la molécula.")


def a_adn(letra):
    """Pol I cambia un ribonucleótido del cebador por el desoxirribonucleótido equivalente."""
    return letra.upper().replace("U", "T")


def simular_replicacion(superior, tam_okazaki=20, tam_cebador=3):
    """Simula una ronda de replicación y devuelve las hebras y los fragmentos de Okazaki."""
    pantalla.titulo("1. REPLICACIÓN DEL ADN  (ADN -> ADN)")
    inferior = complementaria(superior)
    n = len(superior)

    lider = [" "] * n             # hebra nueva que se forma sobre la hebra inferior
    rezagada = [" "] * n          # hebra nueva que se forma sobre la hebra superior
    fragmento_de = [None] * n     # nº de fragmento de Okazaki de cada base de la rezagada
    fragmentos = []               # (inicio, inicio del cebador, fin) de cada fragmento

    if pantalla.mostrar_detalle():
        print("   Leyenda:  ~ hebra sencilla protegida por SSB   minúscula = cebador de ARN"
              "\n             verde = cadena líder   amarillo/azul/magenta = fragmentos de Okazaki")
        print(f"   Escala: fragmentos de Okazaki de {tam_okazaki} nt y cebadores de {tam_cebador} nt"
              "\n   (en E. coli son de 1000-2000 nt y unos 10 nt).")
        pantalla.subtitulo("Molécula de partida: doble hélice cerrada")
        dibujar_horquilla(superior, inferior, lider, rezagada, fragmento_de, 0)
        pantalla.pausa()

    for numero, inicio in enumerate(range(0, n, tam_okazaki)):
        fin = min(inicio + tam_okazaki, n)

        # 1) La helicasa rompe los puentes de hidrógeno de este tramo y las SSB lo cubren.
        for i in range(inicio, fin):
            lider[i] = "~"
            rezagada[i] = "~"

        # 2) La primasa coloca los cebadores de ARN.
        if numero == 0:
            # Un único cebador para la cadena líder, justo en el origen.
            for i in range(0, min(tam_cebador, fin)):
                lider[i] = COMPLEMENTO_ARN[inferior[i]].lower()
        # Un cebador para cada fragmento de la rezagada, en el lado de la horquilla.
        inicio_cebador = max(inicio, fin - tam_cebador)
        for i in range(inicio_cebador, fin):
            rezagada[i] = COMPLEMENTO_ARN[superior[i]].lower()

        if pantalla.mostrar_detalle():
            pantalla.subtitulo(f"Paso {numero + 1}a: la helicasa abre las posiciones {inicio + 1}-{fin}"
                               " y la primasa coloca los cebadores")
            dibujar_horquilla(superior, inferior, lider, rezagada, fragmento_de, fin)
            pantalla.pausa()

        # 3) La ADN polimerasa III alarga los cebadores, siempre en sentido 5'->3'.
        #    Cadena líder: de izquierda a derecha, siguiendo a la horquilla (continua).
        for i in range(inicio, fin):
            if lider[i] == "~":
                lider[i] = COMPLEMENTO_ADN[inferior[i]]
        #    Cadena rezagada: de derecha a izquierda, alejándose de la horquilla,
        #    hasta chocar con el fragmento anterior -> nuevo fragmento de Okazaki.
        for i in range(inicio_cebador - 1, inicio - 1, -1):
            rezagada[i] = COMPLEMENTO_ADN[superior[i]]
            fragmento_de[i] = numero
        fragmentos.append((inicio, inicio_cebador, fin))

        if pantalla.mostrar_detalle():
            pantalla.subtitulo(f"Paso {numero + 1}b: la ADN polimerasa III alarga la líder de forma "
                               f"continua y fabrica el fragmento de Okazaki nº {numero + 1}")
            dibujar_horquilla(superior, inferior, lider, rezagada, fragmento_de, fin)
            pantalla.pausa()

    # 4) La ADN polimerasa I quita los cebadores de ARN y rellena el hueco con ADN.
    cebadores = sum(1 for letra in lider + rezagada if letra.islower())
    lider = [a_adn(letra) for letra in lider]
    for numero, (inicio, inicio_cebador, fin) in enumerate(fragmentos):
        for i in range(inicio_cebador, fin):
            rezagada[i] = a_adn(rezagada[i])
            # El hueco lo rellena Pol I alargando el fragmento siguiente (el de la derecha).
            fragmento_de[i] = numero + 1 if numero + 1 < len(fragmentos) else numero

    # Mellas: entre el ADN de cada fragmento y el que puso Pol I al rellenar su cebador.
    mellas = [inicio_cebador for (_, inicio_cebador, _) in fragmentos[:-1]]

    if pantalla.mostrar_detalle():
        pantalla.subtitulo(f"La ADN polimerasa I sustituye los {cebadores} nucleótidos de ARN de los "
                           "cebadores por ADN")
        dibujar_horquilla(superior, inferior, lider, rezagada, fragmento_de, n)
        print("   Nota: en un cromosoma real los cebadores de los extremos los cubre la otra horquilla"
              "\n   (el origen es bidireccional); aquí se rellenan igualmente para simplificar.")
        pantalla.pausa()

        pantalla.subtitulo(f"La ADN ligasa sella las {len(mellas)} mellas entre fragmentos de Okazaki")
        if n <= 120:
            con_mellas = "".join(("|" if i in mellas else "") + letra for i, letra in enumerate(rezagada))
            print("   Antes   3' " + con_mellas + " 5'   (| = mella)")
            print("   Después 3' " + "".join(rezagada) + " 5'   (cadena continua)")
        else:
            print(f"   Mellas en las posiciones: {', '.join(str(m + 1) for m in mellas)}")
        pantalla.pausa()

    lider, rezagada = "".join(lider), "".join(rezagada)

    # 5) Resultado: replicación semiconservativa.
    pantalla.subtitulo("Resultado: dos moléculas hijas, cada una con una hebra parental y una nueva")
    hija_1_ok = rezagada == complementaria(superior)
    hija_2_ok = lider == complementaria(inferior)
    pantalla.dibujar([
        ("Hija 1: parental", "5'", list(superior), "3'"),
        ("        nueva", "3'", pintar(rezagada, None, "amarillo"), "5'"),
        ("Hija 2: nueva", "5'", pintar(lider, None, "verde"), "3'"),
        ("        parental", "3'", list(inferior), "5'"),
    ], n)
    print(f"   Fragmentos de Okazaki: {len(fragmentos)} | cebadores retirados: {cebadores} | "
          f"mellas selladas: {len(mellas)}")
    print(f"   Todas las bases bien emparejadas (A-T, G-C) en la hija 1: {hija_1_ok}; "
          f"en la hija 2: {hija_2_ok}")
    print(f"   Las dos hijas son idénticas a la molécula original: {hija_1_ok and hija_2_ok}")
    pantalla.pausa()

    return {"superior": superior, "inferior": inferior, "lider": lider, "rezagada": rezagada,
            "fragmentos": fragmentos, "cebadores": cebadores, "mellas": mellas}
