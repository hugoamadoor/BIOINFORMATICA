"""Comprueba que el código reproduce los resultados obtenidos a mano (pytest)."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from dogma import (complementaria, complemento_alineado, replicar, traducir,  # noqa: E402
                   transcribir_desde_codificante, transcribir_desde_molde, validar)
from ej4_splicing import analizar  # noqa: E402

DATA = os.path.join(os.path.dirname(__file__), "..", "data")


def test_ej1_hebra_inferior_del_enunciado():
    assert complemento_alineado("ATGCCGTTAGCT") == "TACGGCAATCGA"


def test_ej1_replicacion_manual():
    r = replicar("ATGCCGTTAGCT")
    assert r["nueva_sobre_superior"] == "AGCTAACGGCAT"
    assert r["nueva_sobre_inferior"] == "ATGCCGTTAGCT"
    assert r["hija_1"] == r["hija_2"]


def test_ej2_transcripcion_desde_ambas_hebras():
    assert transcribir_desde_codificante("ATGCCTGAATGC") == "AUGCCUGAAUGC"
    assert transcribir_desde_molde(complementaria("ATGCCTGAATGC")) == "AUGCCUGAAUGC"


def test_ej3_traduccion_manual():
    assert traducir("AUGUAUGCUUAA") == "MYA"


def test_ej4_isoformas():
    assert not analizar((1, 2, 3, 4, 5))["paro_prematuro"]
    assert analizar((1, 2, 3, 5))["en_fase"]
    assert not analizar((1, 3, 5))["en_fase"]


def test_ej6_insulina_coincide_con_uniprot():
    from Bio import SeqIO
    cds = str(SeqIO.read(os.path.join(DATA, "INS_ENST00000381330_cds.fasta"), "fasta").seq)
    ref = str(SeqIO.read(os.path.join(DATA, "P01308_INS_HUMAN.fasta"), "fasta").seq)
    assert traducir(transcribir_desde_codificante(validar(cds))) == ref
