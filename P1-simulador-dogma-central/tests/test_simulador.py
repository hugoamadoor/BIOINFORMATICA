"""Pruebas del simulador (se ejecutan con: python -m pytest -q).

Biopython solo se usa aquí, para comprobar que nuestras reglas dan el mismo
resultado que una librería de bioinformática de referencia.
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from dogma import pantalla  # noqa: E402
from dogma.codigo_genetico import CODIGO_GENETICO  # noqa: E402
from dogma.replicacion import simular_replicacion  # noqa: E402
from dogma.secuencias import GEN_EJEMPLO, complementaria, leer_fasta, validar  # noqa: E402
from dogma.transcripcion import simular_transcripcion  # noqa: E402
from dogma.traduccion import anticodon, simular_traduccion  # noqa: E402

DATOS = os.path.join(os.path.dirname(__file__), "..", "datos")
pantalla.MODO = "resumen"        # sin pausas ni dibujos intermedios
pantalla.activar_colores(False)


def test_complementariedad():
    assert complementaria("ATGC") == "TACG"
    assert anticodon("AUG") == "UAC"


def test_validar_rechaza_letras_raras():
    assert validar("atg c\n") == "ATGC"
    with pytest.raises(ValueError):
        validar("ATGX")


def test_replicacion_semiconservativa():
    r = simular_replicacion(GEN_EJEMPLO, tam_okazaki=20, tam_cebador=3)
    assert r["lider"] == GEN_EJEMPLO                       # la líder copia la hebra superior
    assert r["rezagada"] == complementaria(GEN_EJEMPLO)    # la rezagada copia la inferior
    assert not any(b.islower() for b in r["lider"] + r["rezagada"])   # no quedan cebadores de ARN


def test_fragmentos_de_okazaki_cubren_toda_la_hebra():
    r = simular_replicacion(GEN_EJEMPLO, tam_okazaki=20, tam_cebador=3)
    fragmentos = r["fragmentos"]
    assert len(fragmentos) == 6                            # 111 nt en bloques de 20
    assert fragmentos[0][0] == 0 and fragmentos[-1][2] == len(GEN_EJEMPLO)
    for (_, _, fin), (inicio_siguiente, _, _) in zip(fragmentos, fragmentos[1:]):
        assert fin == inicio_siguiente                     # sin huecos ni solapes
    assert len(r["mellas"]) == len(fragmentos) - 1         # la ligasa sella una mella entre cada par


def test_transcripcion_del_gen_de_ejemplo():
    t = simular_transcripcion(GEN_EJEMPLO)
    assert t["inicio"] == 35                               # +1 = 6 nt después de la caja -10
    assert t["arnm"].startswith("ACAAGGAGG")
    assert t["arnm"].endswith("UUUUUU")                    # acaba en el terminador
    assert "T" not in t["arnm"]


def test_traduccion_del_gen_de_ejemplo():
    t = simular_transcripcion(GEN_EJEMPLO)
    p = simular_traduccion(t["arnm"])
    assert "-".join(p["proteina"]) == "Met-Ser-Lys-Gly-Glu-Glu-Leu-Phe-Thr-Gly-Val"
    assert p["paro"] == "UAA"


def test_sin_codon_de_inicio():
    assert simular_traduccion("GGGCCCUUU")["proteina"] == []


# --- Comparación con Biopython ---------------------------------------------------

def test_codigo_genetico_igual_que_biopython():
    tabla = pytest.importorskip("Bio.Data.CodonTable").standard_rna_table
    for codon, aminoacido in CODIGO_GENETICO.items():
        if aminoacido == "STOP":
            assert codon in tabla.stop_codons
        else:
            from Bio.SeqUtils import seq3
            assert seq3(tabla.forward_table[codon]) == aminoacido


def test_insulina_igual_que_biopython():
    Seq = pytest.importorskip("Bio.Seq").Seq
    from Bio.SeqUtils import seq1
    _, adn = leer_fasta(os.path.join(DATOS, "insulina_INS_cds.fasta"))
    t = simular_transcripcion(adn)
    p = simular_traduccion(t["arnm"])
    assert t["arnm"] == str(Seq(adn).transcribe())
    assert seq1("".join(p["proteina"])) == str(Seq(adn).translate(to_stop=True))
    assert len(p["proteina"]) == 110
