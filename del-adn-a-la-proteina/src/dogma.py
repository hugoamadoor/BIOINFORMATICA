"""Funciones comunes del dogma central: replicación, transcripción y traducción.

Convención: todas las secuencias se manejan como cadenas en sentido 5'->3',
salvo cuando se indica explícitamente lo contrario para su visualización.
"""
from Bio import SeqIO
from Bio.Seq import Seq

BASES_ADN = set("ACGT")
BASES_ARN = set("ACGU")


def limpiar(secuencia: str) -> str:
    """Quita espacios, guiones y marcas 5'/3' y pasa a mayúsculas."""
    s = secuencia.upper()
    for marca in ("5'", "3'", "5’", "3’"):
        s = s.replace(marca, "")
    return "".join(c for c in s if c.isalpha())


def validar(secuencia: str, tipo: str = "ADN") -> str:
    """Devuelve la secuencia limpia o lanza ValueError si tiene bases no válidas."""
    s = limpiar(secuencia)
    permitidas = BASES_ADN if tipo == "ADN" else BASES_ARN
    invalidas = set(s) - permitidas
    if not s:
        raise ValueError("La secuencia está vacía.")
    if invalidas:
        raise ValueError(f"Bases no válidas para {tipo}: {sorted(invalidas)}")
    return s


def leer_fasta(ruta: str):
    """Lee el primer registro de un FASTA y devuelve (identificador, secuencia)."""
    registro = next(SeqIO.parse(ruta, "fasta"))
    return registro.id, validar(str(registro.seq))


def complementaria(hebra_5a3: str) -> str:
    """Hebra complementaria antiparalela, escrita también en sentido 5'->3'."""
    return str(Seq(hebra_5a3).reverse_complement())


def complemento_alineado(hebra_5a3: str) -> str:
    """Complemento base a base (sin invertir): la hebra de abajo leída 3'->5'."""
    return str(Seq(hebra_5a3).complement())


def duplex(superior_5a3: str, sep: int = 3) -> str:
    """Representa un dúplex antiparalelo en dos líneas, agrupando en codones."""
    def agrupar(s):
        return " ".join(s[i:i + sep] for i in range(0, len(s), sep))
    return (f"5'-{agrupar(superior_5a3)}-3'\n"
            f"3'-{agrupar(complemento_alineado(superior_5a3))}-5'")


def replicar(superior_5a3: str) -> dict:
    """Una ronda de replicación semiconservativa.

    Devuelve las dos hebras nuevas (5'->3') y las dos moléculas hijas, cada una
    formada por una hebra parental y una nueva.
    """
    inferior_5a3 = complementaria(superior_5a3)
    nueva_sobre_superior = complementaria(superior_5a3)   # molde: hebra superior
    nueva_sobre_inferior = complementaria(inferior_5a3)   # molde: hebra inferior
    return {
        "parental_superior": superior_5a3,
        "parental_inferior": inferior_5a3,
        "nueva_sobre_superior": nueva_sobre_superior,
        "nueva_sobre_inferior": nueva_sobre_inferior,
        "hija_1": (superior_5a3, nueva_sobre_superior),
        "hija_2": (nueva_sobre_inferior, inferior_5a3),
    }


def transcribir_desde_codificante(codificante_5a3: str) -> str:
    """ARNm = hebra codificante con U en lugar de T (Bio.Seq.transcribe)."""
    return str(Seq(codificante_5a3).transcribe())


def transcribir_desde_molde(molde_5a3: str) -> str:
    """La ARN polimerasa lee el molde 3'->5' y sintetiza el ARN 5'->3'."""
    return str(Seq(molde_5a3).reverse_complement().transcribe())


def buscar_orf(arnm: str):
    """Primer AUG y primer codón de paro en fase. Devuelve (inicio, paro) en nt (0-based)."""
    inicio = arnm.find("AUG")
    if inicio == -1:
        return None, None
    for i in range(inicio, len(arnm) - 2, 3):
        if arnm[i:i + 3] in ("UAA", "UAG", "UGA"):
            return inicio, i
    return inicio, None


def traducir(arnm: str, tabla: int = 1, hasta_paro: bool = True) -> str:
    """Traduce un ARNm (se recorta al último codón completo)."""
    s = arnm[: len(arnm) - len(arnm) % 3]
    return str(Seq(s).translate(table=tabla, to_stop=hasta_paro))


def contenido_gc(secuencia: str) -> float:
    return 100 * sum(secuencia.count(b) for b in "GC") / len(secuencia)
