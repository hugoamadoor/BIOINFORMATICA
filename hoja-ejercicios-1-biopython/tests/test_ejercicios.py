"""Pruebas sin conexión a Internet (python -m pytest -q)."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from Bio.Seq import Seq  # noqa: E402
from Bio.SeqFeature import FeatureLocation, SeqFeature  # noqa: E402
from Bio.SeqRecord import SeqRecord  # noqa: E402

from ej1_contenido_gc import orf_mas_largo, porcentaje_gc, traduccion_completa  # noqa: E402
from ej3_covid19 import genes_del_registro  # noqa: E402
from ej4_bacterias import traducir_cds  # noqa: E402


def test_porcentaje_gc():
    assert porcentaje_gc(Seq("GGCCAATT")) == 50


def test_orf_mas_largo_ignora_la_utr():
    # 5'UTR con un falso inicio en otro marco, y después el ORF real ATG AAA TTT GGG TAA
    adn = Seq("CCATGTTAG" + "ATGAAATTTGGGTAA" + "CC")
    assert orf_mas_largo(adn) == "MKFG"


def test_traduccion_completa_recorta_codones_incompletos():
    assert traduccion_completa(Seq("ATGAAATAGC")) == "MK*"


def test_genes_en_coordenadas_genbank():
    r = SeqRecord(Seq("A" * 100))
    r.features = [SeqFeature(FeatureLocation(9, 40, strand=1), type="gene", qualifiers={"gene": ["X"]}),
                  SeqFeature(FeatureLocation(9, 40, strand=1), type="CDS", qualifiers={"gene": ["X"]})]
    assert genes_del_registro(r) == [("X", 10, 40)]


def test_cds_en_la_hebra_complementaria():
    # ATG AAA TGA en la hebra -, escrito al revés y complementado en la hebra +
    cds = str(Seq("ATGAAATGA").reverse_complement())
    r = SeqRecord(Seq("GG" + cds + "GG"))
    f = SeqFeature(FeatureLocation(2, 11, strand=-1), type="CDS", qualifiers={})
    assert str(traducir_cds(f, r.seq)) == "MK"


def test_consultar_ncbi_reintenta_si_el_servidor_falla():
    from urllib.error import HTTPError

    from comun import consultar_ncbi

    llamadas = []

    def falla_dos_veces(**parametros):
        llamadas.append(parametros)
        if len(llamadas) < 3:
            raise HTTPError("https://eutils.ncbi.nlm.nih.gov", 500, "Internal Server Error", None, None)
        return "ok"

    assert consultar_ncbi(falla_dos_veces, espera=0, term="x") == "ok"
    assert len(llamadas) == 3
