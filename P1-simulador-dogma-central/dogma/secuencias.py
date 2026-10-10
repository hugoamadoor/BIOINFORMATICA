"""Herramientas básicas para trabajar con secuencias de ADN y ARN.

Convenio de todo el programa: la secuencia que da el usuario es la hebra SUPERIOR
(la codificante), escrita de 5' a 3' y de izquierda a derecha. Las demás hebras se
dibujan alineadas debajo de ella, base frente a base.
"""

# Reglas de Watson y Crick para el ADN: A con T y G con C.
COMPLEMENTO_ADN = {"A": "T", "T": "A", "G": "C", "C": "G"}

# Base del ADN molde -> ribonucleótido que se coloca enfrente al fabricar ARN.
# En el ARN no hay timina: frente a una A del molde se pone uracilo (U).
COMPLEMENTO_ARN = {"A": "U", "T": "A", "G": "C", "C": "G"}


def complementaria(adn):
    """Hebra complementaria base a base (alineada, sin darle la vuelta)."""
    return "".join(COMPLEMENTO_ADN[base] for base in adn)


def adn_a_arn(adn):
    """La misma secuencia escrita como ARN: se cambia T por U."""
    return adn.replace("T", "U")


def validar(texto):
    """Limpia la secuencia (mayúsculas, sin espacios ni números) y comprueba las bases."""
    secuencia = "".join(letra for letra in texto.upper() if letra.isalpha())
    if not secuencia:
        raise ValueError("La secuencia está vacía.")
    no_validas = set(secuencia) - set("ACGT")
    if no_validas:
        raise ValueError(f"La secuencia tiene letras que no son bases de ADN: {sorted(no_validas)}")
    return secuencia


def leer_fasta(ruta):
    """Lee el primer registro de un fichero FASTA y devuelve (nombre, secuencia)."""
    nombre, lineas = "", []
    with open(ruta, encoding="utf-8") as fichero:
        for linea in fichero:
            linea = linea.strip()
            if linea.startswith(">"):
                if lineas:          # ya teníamos un registro completo: nos quedamos con él
                    break
                nombre = linea[1:]
            elif linea:
                lineas.append(linea)
    return nombre or ruta, validar("".join(lineas))


# Gen de ejemplo con la estructura de un gen bacteriano (como en E. coli).
# Cada parte tiene una función concreta en la transcripción o en la traducción.
PARTES_GEN_EJEMPLO = [
    ("Caja -35 del promotor", "TTGACA"),
    ("Separador (17 nt)", "GCTAGCTCAGTCCTAGG"),
    ("Caja -10 del promotor (Pribnow)", "TATAAT"),
    ("Separador hasta el +1", "GCATCG"),
    ("5'UTR con secuencia Shine-Dalgarno", "ACAAGGAGGTAAACC"),
    ("Región codificante (ATG ... TAA)", "ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTTAA"),
    ("3'UTR", "GC"),
    ("Terminador (horquilla + poli-T)", "CCCGCCGAAAGGCGGGTTTTTTT"),
]
GEN_EJEMPLO = "".join(secuencia for _, secuencia in PARTES_GEN_EJEMPLO)
