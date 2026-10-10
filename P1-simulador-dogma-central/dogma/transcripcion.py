"""Transcripción (ADN -> ARN).

Modelo bacteriano simplificado:
  - El factor sigma de la ARN POLIMERASA reconoce el promotor: las cajas -35 (TTGACA)
    y -10 (TATAAT). La transcripción empieza unos nucleótidos después de la caja -10,
    en la posición +1.
  - La ARN polimerasa abre la burbuja de transcripción y lee la HEBRA MOLDE (la inferior)
    de 3' a 5'. Frente a cada base del molde coloca el ribonucleótido complementario
    (A->U, T->A, G->C, C->G), así que el ARN crece de 5' a 3'.
  - El ARN resultante tiene la misma secuencia que la HEBRA CODIFICANTE (la superior),
    pero con U en lugar de T.
  - Termina en un terminador intrínseco: una horquilla rica en G y C seguida de una
    serie de T (poli-U en el ARN). Al llegar ahí, el ARN se suelta.
"""
from . import pantalla
from .secuencias import COMPLEMENTO_ARN, complementaria

CAJA_35 = "TTGACA"
CAJA_10 = "TATAAT"
DISTANCIA_CAJA_10_A_INICIO = 6   # nt entre la caja -10 y el +1 (en E. coli suelen ser 5-9)
SENAL_TERMINACION = "TTTTTT"      # tramo de T del terminador (poli-U en el ARN)
NT_POR_PASO = 15                  # nucleótidos que avanza la ARN polimerasa en cada dibujo


def buscar_promotor(codificante):
    """Devuelve (posición caja -35, posición caja -10, posición +1). -1 si no aparece."""
    caja_10 = codificante.find(CAJA_10)
    if caja_10 == -1:
        return -1, -1, -1
    caja_35 = codificante.rfind(CAJA_35, 0, caja_10)
    inicio = caja_10 + len(CAJA_10) + DISTANCIA_CAJA_10_A_INICIO
    return caja_35, caja_10, inicio


def buscar_terminador(codificante, desde):
    """Posición donde acaba el ARN (justo después del tramo de T), o el final de la molécula."""
    posicion = codificante.find(SENAL_TERMINACION, desde)
    if posicion == -1:
        return -1
    return posicion + len(SENAL_TERMINACION)


def marcas(longitud, etiquetas):
    """Línea de marcas: pone cada texto en su posición (por ejemplo '-35' o '+1')."""
    casillas = [" "] * longitud
    for posicion, texto in etiquetas:
        for j, letra in enumerate(texto):
            if 0 <= posicion + j < longitud:
                casillas[posicion + j] = pantalla.color(letra, "cian")
    return casillas


def dibujar(codificante, molde, arn, etiquetas):
    n = len(codificante)
    pantalla.dibujar([
        ("Codificante", "5'", list(codificante), "3'"),
        ("Molde (leída)", "3'", [pantalla.color(b, "azul") for b in molde], "5'"),
        ("ARNm", "5'", [pantalla.color(b, "verde") for b in arn], "3'"),
        ("", "", marcas(n, etiquetas), ""),
    ], n)


def simular_transcripcion(codificante):
    """Simula la transcripción de la hebra superior y devuelve el ARNm."""
    pantalla.titulo("2. TRANSCRIPCIÓN  (ADN -> ARN)")
    molde = complementaria(codificante)
    n = len(codificante)

    # Iniciación: el factor sigma busca el promotor.
    caja_35, caja_10, inicio = buscar_promotor(codificante)
    if caja_10 == -1 or inicio >= n:
        print("   No se encuentra la caja -10 (TATAAT): se transcribe desde el principio.")
        inicio = 0
    else:
        print(f"   El factor sigma reconoce el promotor: caja -10 (TATAAT) en la posición {caja_10 + 1}"
              + (f" y caja -35 (TTGACA) en la {caja_35 + 1}." if caja_35 != -1 else "."))
        print(f"   La transcripción empieza en la posición {inicio + 1} (+1).")

    fin = buscar_terminador(codificante, inicio)
    if fin == -1:
        print("   No hay terminador (TTTTTT): se transcribe hasta el final de la molécula.")
        fin = n

    etiquetas = [(inicio, "+1")]
    if caja_10 != -1 and inicio > 0:
        etiquetas.append((caja_10, "-10"))
    if caja_35 != -1 and inicio > 0:
        etiquetas.append((caja_35, "-35"))

    arn = [" "] * n
    if pantalla.mostrar_detalle():
        print("   Regla de complementariedad usada (molde -> ARN): A->U  T->A  G->C  C->G")
        pantalla.subtitulo("Iniciación: la ARN polimerasa se une al promotor y abre la burbuja")
        dibujar(codificante, molde, arn, etiquetas)
        pantalla.pausa()

    # Elongación: la ARN polimerasa recorre el molde de 3' a 5' (de izquierda a derecha).
    for tramo in range(inicio, fin, NT_POR_PASO):
        hasta = min(tramo + NT_POR_PASO, fin)
        for i in range(tramo, hasta):
            arn[i] = COMPLEMENTO_ARN[molde[i]]
        if pantalla.mostrar_detalle():
            pantalla.subtitulo(f"Elongación: la ARN polimerasa lee el molde hasta la posición {hasta}")
            dibujar(codificante, molde, arn, etiquetas + [(hasta - 1, "^")])
            print(f"   ARN sintetizado hasta ahora: 5'-{pantalla.abreviar(''.join(arn[inicio:hasta]))}-3'")
            pantalla.pausa()

    arnm = "".join(arn[inicio:fin])

    pantalla.subtitulo("Terminación y resultado")
    if fin < n or codificante.endswith(SENAL_TERMINACION):
        print("   La horquilla del terminador se forma en el ARN y la cola de U lo despega del molde:")
        print("   la ARN polimerasa se suelta y libera el ARNm.")
    print(f"   ARNm ({len(arnm)} nt): 5'-{pantalla.abreviar(arnm)}-3'")
    igual = arnm == codificante[inicio:fin].replace("T", "U")
    print(f"   Comprobación: el ARNm es igual a la hebra codificante cambiando T por U: {igual}")
    pantalla.pausa()

    return {"arnm": arnm, "inicio": inicio, "fin": fin, "caja_10": caja_10, "caja_35": caja_35}
